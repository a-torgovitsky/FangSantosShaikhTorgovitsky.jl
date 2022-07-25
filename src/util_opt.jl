function checkandresolve(m::Model; id::String = "",
                         okcodes::Array{MOI.TerminationStatusCode, 1} =
                            [MOI.OPTIMAL],
                         defaultparams::Dict = Dict())
    ok = checksolveresult(m, id = id, okcodes = okcodes)
    if ok
        return ok
    else
        return resolve(m, id = id,
                       okcodes = okcodes, defaultparams = defaultparams)
    end
end

function checksolveresult(m::Model; id::String = "",
                          okcodes::Array{MOI.TerminationStatusCode, 1} =
                          [MOI.OPTIMAL])
    status = termination_status(m)
    ok = (status in okcodes)
    @debug "$id termination status: $status."
    if !ok
        @debug "$id termination status was $status."
    end
    return ok
end

function resolve(m::Model; id::String = "",
                 okcodes::Array{MOI.TerminationStatusCode, 1} = [MOI.OPTIMAL],
                 defaultparams::Dict = Dict())

    # These loop sets are hard-coded based on Gurobi options -- see manual
    methods = Dict("primal simplex" => 0,
                   "dual simplex" => 1,
                   "barrier" => 2)
    @debug "Trying to resolve $id."
    for sf in -1:3 # ScaleFlag
        for nf in 0:3 # NumericFocus
            for (name, number) in methods
                @debug "Trying: \n"*
                       "\tMethod = $name\n"*
                       "\tNumericFocus = $nf\n"*
                       "\tScaleFlag = $sf"
                set_parameter(m, "Method", number)
                set_parameter(m, "NumericFocus", nf)
                set_parameter(m, "ScaleFlag", sf)
                optimize!(m)
                newstatus = termination_status(m)
                setparams(m, defaultparams) # restore defaults
                if newstatus in okcodes
                    @debug "Problem solved, resuming normal operation."
                    return true
                else
                    @debug "Solve result was $(string(newstatus))."
                end
            end
        end
    end

    @debug "Changing methods did not fix things in $id."
    return false
end

function setparams(m::Model, dp::Dict)
    for (key, value) in dp
        set_parameter(m, key, value)
    end
end
