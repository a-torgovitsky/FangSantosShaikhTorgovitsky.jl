@with_kw struct LPModel
    A_obs::AbstractArray{<:Real, 2}
    beta_obs::Union{Function, Nothing}
    A_shp::AbstractArray{<:Real, 2}
    beta_shp::AbstractArray{<:Real, 1}
    A_tgt::AbstractArray{<:Real, 2}

    function LPModel(aobs, betaobs, ashp, betashp, atgt)
        if ~(size(aobs)[2] == size(ashp)[2] == size(atgt)[2])
            error("A_obs, A_shp, and A_tgt should have the "*
                  "same number of columns. But they have sizes "*
                  "Aobs = $(size(aobs)), "*
                  "A_shp = $(size(ashp)), "*
                  "A_tgt = $(size(atgt))")
        end
        if ~(size(ashp)[1] == length(betashp))
            error("Number of rows in A_shp is $(size(ashp[1])), "*
                  "but length of beta_shp is $(size(betashp)).")
        end

        new(aobs, betaobs, ashp, betashp, atgt)
    end
end
export LPModel

function add_shape(lpm::LPModel, m::Model, varname::Symbol)
    x = m[varname]
    for i = 1:size(lpm.A_shp)[1]
        @constraint(m, sum(lpm.A_shp[i,j]*x[j] for j in 1:size(lpm.A_shp)[2])
                        == lpm.beta_shp[i])
    end
    return m
end

function add_target(lpm::LPModel, m::Model, varname::Symbol)
    if size(lpm.A_tgt)[1] > 1
        error("Not coded yet for more than 1-dimensional target parameters.")
    end
    @variable(m, testpoint == 0)
    testpointsymbol = :testpoint
    x = m[varname]
    @constraint(m, sum(lpm.A_tgt[1,j]*x[j] for j in 1:size(lpm.A_shp)[2])
                   == testpoint)
    return m, testpointsymbol
end

function feas_check(lpm::LPModel, t::Real)
    m = Model(with_optimizer(Gurobi.Optimizer, GLOBAL_GENV[],
                             OutputFlag = 0))

    dm, dx = size(lpm.A_obs) # dm = # moments, dx = # parameters
    @variable(m, 0. <= x[1:dx])
    m = add_shape(lpm, m, :x)
    m, testpointsymbol = add_target(lpm, m, :x)
    fix(m[testpointsymbol], t)

    optimize!(m) # JuMP default is to check feasibility

    ok = checksolveresult(m, id = "Feasibility check",
                          okcodes = [MOI.OPTIMAL, MOI.INFEASIBLE])
    if !ok
        @error "Feasibility check failed for t = $t."
    end
    if termination_status(m) == MOI.OPTIMAL
        return true
    else
        return false
    end
end

function feas_check(lpm::LPModel, tlist::AbstractArray{<:Real,1})
    feas = fill(true, length(tlist))
    for (i,t) in enumerate(tlist)
        feas[i] = feas_check(lpm, t)
    end
    return feas
end

function compute_bounds(lpm::LPModel, betaobs::AbstractArray{<:Real, 1})
    m = Model(with_optimizer(Gurobi.Optimizer, GLOBAL_GENV[],
                             OutputFlag = 0))

    dx = size(lpm.A_obs, 2)
    @variable(m, x[1:dx] >= 0)

    @assert size(lpm.A_obs, 1) == length(betaobs)
    @constraint(m, matchdata[j=1:size(lpm.A_obs, 1)],
                sum(lpm.A_obs[j,i]*x[i] for i in 1:dx) == betaobs[j])
    m = add_shape(lpm, m, :x)

    bounds = NaN*ones(2)
    opt = NaN*ones(dx, 2)
    for (i,s) in enumerate([MOI.MIN_SENSE, MOI.MAX_SENSE])
        @objective(m, s, sum(lpm.A_tgt[i]*x[i] for i in 1:dx))
        optimize!(m)
        status = checksolveresult(m, id = "compute_bounds",
                                  okcodes = [MOI.OPTIMAL, MOI.INFEASIBLE])
        opt[:,i] = value.(x)
        if status == MOI.INFEASIBLE # max will also be infeasible
            bounds = [+Inf, -Inf]
            break
        else
            bounds[i] = objective_value(m)
        end
    end

    return Tuple(bounds), opt
end
