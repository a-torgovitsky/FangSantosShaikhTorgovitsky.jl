@with_kw struct MonteCarloParams
    mcreps::Integer = 10
    seeddgp::Integer = 123456
    parallel::Bool = true
    maxfails::Integer = 20
end
export MonteCarloParams

@with_kw struct Procedure
    f::Function = (x -> println("you need to fill this!"))
    name::String = "unnamed"
end
export Procedure

# Many procedures, single dgp
function runmc(procs::Array{Procedure, 1}, draw::Function,
               mcp::MonteCarloParams = MonteCarloParams())

    results = DataFrame()
    if length(procs) == 0
        println("Nothing to do: length of procs is 0.")
        return results
    end

    tlist = [(procs, draw, mcp.seeddgp, mcp.maxfails, m)
             for m in 1:mcp.mcreps]
    if (mcp.parallel) & (nprocs() > 1)
        @info "Running in parallel with $(nprocs()) processes."
        resultslist = pmap(onedraw_manyprocs, tlist)
    else
        resultslist = map(onedraw_manyprocs, tlist)
    end
    results = vcat(resultslist...)

    return results
end

function onedraw_manyprocs(t::Tuple{Array{Procedure, 1},
                                    Function,
                                    Integer,
                                    Integer,
                                    Integer})
    onedraw_manyprocs(t[1], t[2], t[3], t[4], t[5])
end

function onedraw_manyprocs(procs::Array{Procedure, 1},
                           draw::Function,
                           seedbase::Integer,
                           maxfails::Integer,
                           repnum::Integer)

    results = DataFrame()
    success = false
    procfails = 0
    while !success
        s = seedbase + repnum + procfails*Int(1e6)
        id = "MC rep $repnum"
        if procfails > 0
            id = id*", retry $procfails"
        end
        @debug id
        data, drawfails = draw(s)

        try
            for (i,p) in enumerate(procs)
                r = p.f(data, id)
                r = hcat(DataFrame(rep = fill(repnum, nrow(r)),
                                   proc = fill(p.name, nrow(r))), r)
                r[:, :r_f_drawfails] = fill(drawfails, nrow(r))
                r[:, :r_f_procfails] = fill(procfails, nrow(r))
                if isempty(results)
                    results = similar(r, 0)
                else
                    harmonize!(results, r)
                end
                append!(results, r)
                success = true
            end
        catch err
            str = "repetition $repnum, attempt $procfails."
            if isa(err, ConeError) | isa(err, GMMError)
                @warn "Handled error in "*str
                procfails += 1
            else
                @error "Unexpected error in "*str
                throw(err)
            end
        end
        if procfails > maxfails
            @error "Too many failures $procfails, allowed $maxfails."*
                   " Something is probably wrong..."
        end
    end
    return results
end

function simulate_baddraws(draw::Function, mcp::MonteCarloParams)
    getdrawfails(m) = draw(mcp.seeddgp + m)[2]
    countlist = fill(0, mcp.mcreps)
    if (mcp.parallel) & (nprocs() > 1)
        countlist = pmap(getdrawfails, 1:mcp.mcreps)
    else
        countlist = map(getdrawfails, 1:mcp.mcreps)
    end
    return sum(countlist) ./ mcp.mcreps
end
