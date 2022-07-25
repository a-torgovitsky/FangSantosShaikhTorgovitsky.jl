@with_kw struct GMMParams
    bsreps::Int = 500
    bsextra::Real = .2 # as a fraction of bsreps
    seedstart::Int = 1
    recenter::Bool = false
end

struct GMMError <: Exception
    var::String
end
Base.showerror(io::IO, e::GMMError) = print(io, "GMMError: ", e.var)

function constrained_gmm(betahatobs::Array{<:Real, 1},
                         vc::Array{<:Real, 2},
                         lpm::LPModel,
                         recenter::Array{<:Real, 1} = [NaN],
                         id::String = "")

    ############################################################################
    # Set up the problem
    ############################################################################
    A = vcat(lpm.A_obs, lpm.A_shp, lpm.A_tgt)
    dimbeta, dimx = size(A)
    dimbetaobs = size(lpm.A_obs)[1]
    if any(isnan.(recenter))
        recenter = zeros(dimbetaobs)
    end
    betahatobs = betahatobs - recenter

    function create_beta(b::AbstractArray{<:Real, 1}, t::Real)
        return vcat(b, lpm.beta_shp, t)
    end

    # in the paper notation
    # p = dimbeta
    # d = dimx
    if (dimbeta < dimx)
        throw(GMMError("This is intended as a method only for the p > d case."))
    end

    ############################################################################
    # Constrained GMM estimator
    ############################################################################
    m = Model(with_optimizer(Gurobi.Optimizer,
                             GLOBAL_GENV[],
                             NumericFocus = 1,
                             OutputFlag = 0))
    @variable(m, x[1:dimx] >= 0)
    m = add_shape(lpm, m, :x)

    # Quadratic objective
    # x'Bx + c'x + d
    #
    # (Aobs * x - betahatobs)' * Xi * (Aobs * x - betahatobs)
    # =
    # x' * (Aobs' * Xi * Aobs) * x
    # +
    # x'(-2 * Aobs' * Xi * betahatobs)
    # +
    # betahatobs' * Xi * betahatobs
    Xi = inv(vc)
    B = lpm.A_obs' * Xi * lpm.A_obs
    c = -2 * lpm.A_obs' * Xi * betahatobs
    d = betahatobs' * Xi * betahatobs
    @objective(m, Min,
                 sum(x[i]*x[j]*B[i,j] for i = 1:dimx, j = 1:dimx)
               + sum(x[i]*c[i] for i = 1:dimx) + d)

    optimize!(m)
    ok = checkandresolve(m, id = id)
    xcstar = value.(x)

    return xcstar, ok
end

function wald_bootstrap(data::DataFrame,
                        testptlist::Array{<:Real, 1},
                        lpm::LPModel;
                        p::GMMParams = GMMParams(),
                        baseid::String = "")

    betahatobs, vc = lpm.beta_obs(data)
    id = baseid * " constrained GMM, sample."
    xcstar, f_sample = constrained_gmm(betahatobs, vc, lpm)
    if !f_sample
        throw(GMMError("Constrained GMM failed in sample."))
    end
    betahat_tgt = lpm.A_tgt * xcstar

    if p.recenter
        rc = lpm.A_obs * xcstar - betahatobs
    else
        rc = zeros(size(lpm.A_obs, 1))
    end

    bsmax = Int(round(p.bsreps * (1 + p.bsextra)))
    xcstar_bs = fill(NaN, (length(xcstar), p.bsreps)) # not needed to save
    tgt_bs = fill(NaN, p.bsreps)
    b = 1
    btry = 1
    f_redraw = 0
    f_bsfail = 0
    while b <= p.bsreps
        data_bs = resampledata(data, nrow(data), p.seedstart + btry)
        betahatobs_bs, _ = lpm.beta_obs(data_bs)
        if !checkresample(betahatobs_bs, betahatobs)
            @warn "Redraw $b did not pass checkresample."
            f_redraw += 1
        else
            id = baseid * "constrained GMM, b = $b, btry = $btry"
            x, ok = constrained_gmm(betahatobs_bs, vc, lpm, rc)
            if ok
                xcstar_bs[:, b] = x
                tgt_bs[b] = lpm.A_tgt * x
                b += 1
            else
                f_bsfail += 1
            end
        end
        if btry > bsmax
            throw(GMMError("Too many failures."))
        end
        btry += 1
    end

    results = wald_rule(betahat_tgt, var(tgt_bs), testptlist)
    results[!, :f_redraw] = fill(f_redraw, nrow(results))
    results[!, :f_bsfail] = fill(f_bsfail, nrow(results))
    return results
end

function wald_rule(tgt, var, testptlist)
    results = DataFrame(
        t = testptlist,
        r_wald = fill(NaN, length(testptlist)),
        r_pval = fill(NaN, length(testptlist)),
        r_cv10 = fill(quantile(Chisq(1),.90), length(testptlist)),
        r_cv05 = fill(quantile(Chisq(1),.95), length(testptlist)),
        r_cv01 = fill(quantile(Chisq(1),.99), length(testptlist))
    )
    for (i,t) in enumerate(testptlist)
        results[i, :r_wald] = (var^-1) * (tgt - t)^2
        results[i, :r_pval] = ccdf(Chisq(1), results[i, :r_wald])
    end
    return results
end

function gmm_rank_check(lpm::LPModel)
    rk = rank(lpm.A_obs)
    if rk < size(lpm.A_obs, 2)
        return false, rk
    else
        return true, rk
    end
end
