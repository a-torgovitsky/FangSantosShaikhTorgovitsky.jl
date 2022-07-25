function menu(savelocation::String = ".",
              choice::Int = -1;
              quicktest::Bool = false,
              autoname::Bool = false)
    ConsoleLogger(stdout, Logging.Info) |> global_logger

    printdiv('=')
    print(read(joinpath(@__DIR__, "banner.txt"), String))
    printdiv('=')
    println("Inference for Large-Scale Systems of Linear Inequalities")
    println("Fang, Santos, Shaikh, and Torgovitsky (2021)")
    printdiv('=')
    if choice <= 0
        println("Which simulation would you like to run?")
        println("\t 1. Identified sets (Figure 3)")
        println("\t 2. Comparison with Wald in nearly point identified cases (Figure 4)")
        println("\t 3. Size of FSST in many designs with λᵇ (Table 1)")
        println("\t 4. Size of FSST in many designs with λʳ (Table 2)")
        println("\t 5. Power (Figure 5)")
        print("Enter choice: ")
        choice = readline()
        choice = parse(Int64, choice)
    end
    @info "Running simulation $choice."
    if quicktest
        @info "TESTING MODE: SMALL NUMBER OF REPLICATIONS"
        mcreps = 2
        mcreps_small = 2
        bsreps = 5
    else
        mcreps = MCREPS
        mcreps_small = MCREPS_SMALL
        bsreps = BSREPS
    end

    fn = ""
    if choice == 1
        resultsdir = setup(savelocation;
                           stub = "mixedlogit_dfidset", autoname = autoname)
        fn = idsets(resultsdir);
        calloutput(dirname(fn), "plot_idset();")
    elseif choice == 2
        resultsdir = setup(savelocation;
                           stub = "mixedlogit_nearid", autoname = autoname)
        fn = simnearid(resultsdir, bsreps = bsreps, mcreps = mcreps);
        calloutput(dirname(fn), "plot_nearid();")
        calloutput(dirname(fn), "plot_errors();")
    elseif choice == 3
        resultsdir = setup(savelocation;
                           stub = "mixedlogit_size_auto", autoname = autoname)
        fn = simsize(resultsdir, bsreps = bsreps, mcreps = mcreps,
                     mcreps_small = mcreps_small);
        calloutput(dirname(fn), "tables_size();")
        calloutput(dirname(fn), "plot_errors();")
    elseif choice == 4
        resultsdir = setup(savelocation;
                           stub = "mixedlogit_size_mult", autoname = autoname)
        fn = simsize(resultsdir, cone_lambda = ["mult"], bsreps = bsreps,
                     mcreps = mcreps, mcreps_small = mcreps_small);
        calloutput(dirname(fn), "tables_size();")
        calloutput(dirname(fn), "plot_errors();")
    elseif choice == 5
        resultsdir = setup(savelocation;
                           stub = "mixedlogit_power", autoname = autoname)
        fn = simpower(resultsdir, bsreps = bsreps, mcreps = mcreps);
        calloutput(dirname(fn), "plot_power();")
        calloutput(dirname(fn), "plot_errors();")
    else
        error("Bad entry. Goodbye.")
    end
    return fn
end
export menu

function calloutput(dir::String, functioncall::String)
    srcdir = @__DIR__
    plotcommand = `Rscript -e "setwd('$srcdir');
                   source('output.R');
                   setwd('$(abspath(dir))');
                   $(functioncall)"`;
    run(plotcommand)
    println("Hurry up and check out $(dir)!!!")
end
export calloutput
