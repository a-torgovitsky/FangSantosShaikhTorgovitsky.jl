@with_kw struct Design
    dgp = []
    samplesize = 1000
    bsreps = 500
    mcreps = 500
    seedinit = 123456
    testptlist = [] # empty means test lower bound, see below
    cone_lambdalist = [.5]
    gmm = false
    skiptszero = false
    compute_bounds = true
end
export Design

DESIGNINFO_FN = "designinfo.csv"
RESULTS_FN = "results.csv"

function rundesigns(resultsdir,
                    designlist::Array{Design, 1};
                    debug::Bool = false)
    # Function to add time stamp to log output
    timestamp_logger(logger) = TransformerLogger(logger) do log
        merge(log, (; message =
                    "$(Dates.format(now(), "mm-dd HH:MM:SS")) $(log.message)"))
    end
    # Log everything to both console and log.out
    loglevel = (debug ? Logging.Debug : Logging.Info)
    TeeLogger(
        MinLevelLogger(FileLogger(joinpath(resultsdir, "log.out")),
                       loglevel) |> timestamp_logger,
        ConsoleLogger(stdout, loglevel) |> timestamp_logger
    ) |> global_logger

    fn_results = joinpath(resultsdir, RESULTS_FN)
    if isfile(fn_results)
        @info "Results file already exists; trying to resume"
        designlist, dsum = recover_designs(designlist, resultsdir)
        results = CSV.read(fn_results)
    else
        dsum = summarize_designs(designlist, resultsdir)
        results = DataFrame()
    end

    for (i, d) in enumerate(designlist)
        _, designname = design_to_id(d)
        @info '='^length(designname)
        @info "Starting $(designname)"

        if isempty(d.testptlist)
            testptlist = [dsum[i, :lb]]
        else
            testptlist = d.testptlist
        end

        # Generate procedures to run
        procs = design_to_proclist(d, testptlist)

        # Run them
        r = runmc(procs, design_to_draw(d), design_to_mcp(d))

        # Organize results
        r = hcat(DataFrame(design = fill(i, nrow(r))), r)
        if isempty(results)
            results = similar(r, 0)
        else
            harmonize!(results, r)
        end
        append!(results, r)
        if isfile(fn_results)
            rm(fn_results)
        end
        CSV.write(fn_results, results, writeheader = true)
        @info "Done ($(i)/$(length(designlist)))"
    end
    @info "Done with everything!"
    @info "Results are in $resultsdir"
    return fn_results
end
export rundesigns

function design_to_id(d::Design)
    dimv = size(d.dgp.vsupp, 1)
    dimw = size(d.dgp.wsupp, 1)
    id = Dict{Symbol, Any}(:dimvsupp => dimv,
                           :d => dimv,
                           :dimwsupp => dimw,
                           :p => dimw + 2,
                           :samplesize => d.samplesize,
                           :targetparam => d.dgp.tp.printname)
    name = "mixedlogit_dv=$(dimv)_dw=$(dimw)_tp=$(d.dgp.tp.printname)"*
           "_n=$(d.samplesize)"
    return id, name
end

function design_to_proclist(d::Design, testptlist::Array{<:Real, 1})
    lpm = dgp_to_lpm(d.dgp)

    # Only do GMM in "nearly" point identified cases
    do_gmm = (size(lpm.A_obs)[1] >= size(lpm.A_obs)[2]) * d.gmm
    proclen = 1 + # for cone
              2 * Int(do_gmm)
    proclist = Array{Procedure, 1}(undef, proclen)

    # Cone
    cp = ConeParams(bsreps = d.bsreps, lambda = d.cone_lambdalist)
    ct = ((data, baseid) -> conetest(data, testptlist, lpm, p = cp,
                                     baseid = baseid))
    name = "cone"
    proclist[1] = Procedure(f = ct, name = name)

    # GMM
    if do_gmm
        pgmm1 = GMMParams(bsreps = d.bsreps, recenter = false)
        cgmm = ((data, baseid) -> wald_bootstrap(data, testptlist, lpm,
                                             p = pgmm1, baseid = baseid))
        name = "gmm bs"
        proclist[2] = Procedure(f = cgmm, name = name)

        pgmm2 = GMMParams(bsreps = d.bsreps, recenter = true)
        cgmm = ((data, baseid) -> wald_bootstrap(data, testptlist, lpm,
                                             p = pgmm2, baseid = baseid))
        name = "gmm bs rc"
        proclist[3] = Procedure(f = cgmm, name = name)
    end

    return proclist
end

function design_to_draw(d::Design)
    (seed) -> drawsample_safe(d.dgp, d.samplesize, seed,
                              scramble = d.mcreps)
end

function design_to_mcp(d::Design)
    MonteCarloParams(mcreps = d.mcreps, seeddgp = d.seedinit)
end

function summarize_designs(designlist::Array{Design, 1}, resultsdir::String)
    dsum = DataFrame(design = 1:length(designlist))
    dsum.p = [(size(d.dgp.wsupp)[1] + 2) for d in designlist]
    dsum.d = [size(d.dgp.vsupp)[1] for d in designlist]
    dsum.samplesize = [d.samplesize for d in designlist]
    dsum.bsreps = [d.bsreps for d in designlist]
    dsum.mcreps = [d.mcreps for d in designlist]
    id = hcat([[compute_true_bounds(d.dgp, d.dgp.tp)...]
               for d in designlist]...)
    dsum.lb = round.(id[1,:], digits = 3)
    dsum.ub = round.(id[2,:], digits = 3)
    dsum.prbaddraws = round.([simulate_baddraws(design_to_draw(d),
                                                design_to_mcp(d))
                              for d in designlist], digits = 3)
    dsum.minpr = round.(minimum.([compute_true_choice_prob(d.dgp)
                                  for d in designlist]), digits = 3)
    dsum.maxpr = round.(maximum.([compute_true_choice_prob(d.dgp)
                                  for d in designlist]), digits = 3)
    # More stuff, saved but not printed
    dsum.gmm_fullrank = [gmm_rank_check(dgp_to_lpm(d.dgp))[1]
                                for d in designlist]
    dsum.optimal_lambda1 = [optimal_lambda1(dsum.samplesize[i], dsum.p[i])
                            for i in 1:nrow(dsum)]
    dsum.optimal_lambda2 = [optimal_lambda2(dsum.samplesize[i], dsum.p[i])
                            for i in 1:nrow(dsum)]
    dsum.targetparam = [d.dgp.tp.printname for d in designlist]

    print_designs(dsum)
    fn = joinpath(resultsdir, DESIGNINFO_FN)
    CSV.write(fn, dsum)
    return dsum
end

function recover_designs(designlist::Array{Design, 1}, resultsdir::String)
    @info "Attemping to recover designs."

    fnd = joinpath(resultsdir, DESIGNINFO_FN)
    fnr = joinpath(resultsdir, RESULTS_FN)
    if isfile(fnd) & isfile(fnr)
        @info "Found both info and results file."
        dsum = CSV.read(fnd)
        r = CSV.read(fnr)
    else
        @info "Didn't find info or results file, starting over"
        return designlist, summarize_designs(designlist, resultsdir)
    end

    creps = fill(NaN, nrow(dsum))
    for i in 1:nrow(dsum)
        creps[i] = nrow(r[(r[!, :p] .== dsum[i, :p]) .&
                          (r[!, :d] .== dsum[i, :d]) .&
                          (r[!, :samplesize] .== dsum[i, :samplesize]), :])
    end
    cidx = (creps .== dsum[!, :mcreps]) # indices of completed simulations

    print_designs(dsum[cidx, :], title = "Designs already COMPLETED.")
    print_designs(dsum[.!cidx, :], title = "Designs remaining to be run.")

    return (designlist[.!cidx], dsum[.!cidx, :])
end

function print_designs(dsum; title = "Summary of the designs to be run")
    pretty_table(dsum[:, 2:11];
                 header = ["p", "d", "n", "#BS", "#MC", "LB",
                           "UB", "Pr Bad", "Min Pr", "Max Pr"],
                 title = title,
                 show_row_number = true,
                 row_number_column_title = "Design #",
                 crop = :none)
end
