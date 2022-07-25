function setup(savelocation::String;
               stub::String="unnamed",
               autoname = false)

    @info "Saving in: $savelocation"
    tag = ""

    # Build directory name
    if !autoname
        print("Enter a directory tag "*
              "(or hit Enter to generate from the date): ")
        tag = readline()
    end
    if isempty(tag)
        tag = Libc.strftime("%d-%b-%y-at-%H-%M-%S", time())
    end
    if isempty(stub)
        dirname = tag
    else
        dirname = stub*"-"*tag
    end
    resultsdir = joinpath(savelocation, dirname)

    # Create directories
    if isdir(resultsdir)
        @warn "Simulation directory already exists."
    else
        mkdir(resultsdir)
    end

    @info "Saving in $resultsdir"
    return resultsdir
end
