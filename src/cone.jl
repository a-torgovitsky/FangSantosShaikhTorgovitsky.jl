@with_kw struct ConeParams
    lambda::Array{Any, 1} = ["auto"]
    rho::Real = 1e-4
    bsreps::Int = 100
    bswiggle::Real = .2 # as a fraction of bsreps
    seedstart::Int = 1
    imgtol::Real = 1e-5

    function ConeParams(lambda, rho, bsreps, bswiggle, seedstart, imgtol)

        if ~(0 <= rho)
            throw(DomainError(rho))
        end
        if ~(bsreps > 0)
            throw(DomainError(bsreps))
        end
        if ~(bswiggle >= 0)
            throw(DomainError(bswiggle))
        end
        if ~(seedstart > 0)
            throw(DomainError(seedstart))
        end
        if ~(imgtol > 0)
            throw(DomainError(imgtol))
        end

        new(lambda, rho, bsreps, bswiggle, seedstart, imgtol)
    end
end

struct ConeError <: Exception
    var::String
end
Base.showerror(io::IO, e::ConeError) = print(io, "ConeError: ", e.var)

function conetest(data::DataFrame,
                  testptlist::Array{<:Real, 1},
                  lpm::LPModel;
                  p::ConeParams = ConeParams(),
                  baseid::String = "")

    r = fill(DataFrame(), length(testptlist))
    for (i,t) in enumerate(testptlist)
        r[i] = conetest(data, t, lpm, p = p, baseid = baseid)
    end
    results = vcat(r...)
    return results
end

function conetest(data::DataFrame,
                  testpt::Real,
                  lpm::LPModel;
                  p::ConeParams = ConeParams(),
                  baseid::String = "")

    ############################################################################
    # Feasibility check
    ############################################################################
    feas = feas_check(lpm, testpt)
    if ~feas # No point in continuing since will always reject
        @warn "Test point was infeasible. Returning empty results."
        return create_results_df(p.lambda)
    end

    ############################################################################
    # Set up the problem
    ############################################################################
    n = nrow(data)
    rootn = sqrt(n)

    # in the paper notation
    # dimbeta = p
    # dimx = d
    A = vcat(lpm.A_obs, lpm.A_shp, lpm.A_tgt)
    dimbeta, dimx = size(A)
    dimbetaobs = size(lpm.A_obs)[1]

    if (dimx >= dimbeta) # A has more columns than rows
        flagrange = false
    else
        flagrange = true
    end

    # This function is just for convenience
    function create_beta(b::AbstractArray{<:Real, 1})
        return vcat(b, lpm.beta_shp, testpt)
    end

    @debug "Sample size is $n."
    @debug "Dimension of x (columns of A) is $dimx."
    @debug "Dimension of beta (rows of A) is $dimbeta."
    @debug "Flagrange (more rows than columns) is $flagrange."

    ############################################################################
    # Sample and bootstrapped betahatobs, stored for later
    ############################################################################
    # Note: betahatobs_vc should be the asymptotic variance (not finite sample)
    betahatobs, betahatobs_vc = lpm.beta_obs(data)
    betahat = create_beta(betahatobs)

    # Draw bsmax > p.bsreps in case we have problems with certain draws
    bsmax = Int(round(p.bsreps * (1 + p.bswiggle)))
    betahatobs_bs = fill(NaN, (length(betahatobs), bsmax))
    f_redraw = 0
    b = 1
    btry = 1
    while btry <= bsmax
        data_bs = resampledata(data, n, p.seedstart + b)
        bobs, _ = lpm.beta_obs(data_bs)
        if checkresample(bobs, betahatobs)
            betahatobs_bs[:, btry] = bobs
            btry += 1
        else
            @debug "Redraw $b did not pass checkresample."
            f_redraw += 1
        end
        b += 1
    end

    ############################################################################
    # Estimation and bootstrapping of betahatstar
    ############################################################################
    betahatobs_star = fill(NaN, (length(betahatobs)))
    betahatobs_star_bs = fill(NaN, (length(betahatobs), p.bsreps))
    f_betahatobs_star_bs = 0

    if !flagrange # dimx >= dimbeta, the easy case
        betahatobs_star = betahatobs
        betahatobs_star_bs = betahatobs_bs[:, 1:p.bsreps]
    else # dimx < dimbeta
        m_betahatstar = Model(with_optimizer(Gurobi.Optimizer, GLOBAL_GENV[]))
        dp_betastar = Dict("OutputFlag" => 0,
                           "Method" => 0, # Primal simplex seems best
                           "NumericFocus" => 0,
                           "ScaleFlag" => -1)
        setparams(m_betahatstar, dp_betastar)
        @variable(m_betahatstar, x[1:dimx])
        m_betahatstar = add_shape(lpm, m_betahatstar, :x)
        m_betahatstar, testpointsymbol = add_target(lpm, m_betahatstar, :x)
        fix(m_betahatstar[testpointsymbol], testpt)
        Xi = inv(betahatobs_vc) # Weighting matrix
        ########################################################################
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
        ########################################################################
        B = lpm.A_obs' * Xi * lpm.A_obs
        cpre = -2 * lpm.A_obs' * Xi
        function solve_betahatstar(m::Model, beta::Array{Float64, 1}, id)
            c = cpre * beta
            d = beta' * Xi * beta
            @objective(m, Min,
                         sum(x[i]*x[j]*B[i,j] for i = 1:dimx, j = 1:dimx)
                       + sum(x[i]*c[i] for i = 1:dimx) + d)
            optimize!(m)
            ok = checkandresolve(m, id = id, defaultparams = dp_betastar)
            xhatstar = value.(x)
            return xhatstar, ok
        end

        id = baseid * ": betahatstar, t = $testpt, sample"
        xhatstar, ok = solve_betahatstar(m_betahatstar, betahatobs, id)
        if !ok
            throw(ConeError("$id did not solve correctly."))
        end
        betahatobs_star = lpm.A_obs * xhatstar

        b = 1
        btry = 1
        while b <= p.bsreps
            id = baseid * ": betahatstar, t = $testpt, bootstrap $b, try $btry"
            xhatstar, ok = solve_betahatstar(m_betahatstar,
                                             betahatobs_bs[:, btry],
                                             id)
            if ok
                betahatobs_star_bs[:, b] = lpm.A_obs * xhatstar
                b += 1
            else
                f_betahatobs_star_bs += 1
            end
            btry += 1
            if btry > bsmax
                throw(ConeError("Too many failures in"*
                                " Betahatstar t = $testpt."))
            end
        end
    end
    betahat_star = create_beta(betahatobs_star)
    betahat_star_bs = fill(NaN, (length(betahat), p.bsreps))
    for b = 1:p.bsreps
        betahat_star_bs[:,b] = create_beta(betahatobs_star_bs[:, b])
    end

    ############################################################################
    # Studentization
    #   the p.rho terms gets used as a ridge regularization to improve
    #   numerical stability
    ############################################################################
    Eye = Matrix{Float64}(I,dimbeta,dimbeta)
    OmegaHatI = zeros(dimbeta, dimbeta)

    if flagrange
        # Variance of betahatstar and then OmegaHatI
        diffstar = betahatobs_star_bs .- betahatobs_star
        diffstar = vcat(diffstar, zeros(length(lpm.beta_shp) + 1, p.bsreps))

        OmegaHatI = diffstar * diffstar'
        OmegaHatI *= (n/p.bsreps)
        rhobar = norm(OmegaHatI) * p.rho
        OmegaHatI = safe_sqrt(OmegaHatI + rhobar * Eye, p.imgtol)
    else
        dimrem = dimbeta - dimbetaobs
        betahatstar_vc = [betahatobs_vc             zeros(dimbetaobs, dimrem);
                          zeros(dimrem, dimbetaobs) zeros(dimrem, dimrem)      ]
        rhobar = norm(betahatstar_vc) * p.rho
        OmegaHatI = safe_sqrt(betahatstar_vc + rhobar * Eye, p.imgtol)
    end

    ############################################################################
    # Compute the restricted estimator
    ############################################################################
    m_restricted = Model(with_optimizer(Gurobi.Optimizer,
                                        GLOBAL_GENV[],
                                        OutputFlag = 0))

    @variable(m_restricted, b[1:dimbeta])
    @variable(m_restricted, x[1:dimx] >= 0)
    @variable(m_restricted, phi0)
    @variable(m_restricted, phi1[1:2] >= 0)
    @variable(m_restricted, phid[1:2, 1:dimx] >= 0)
    @variable(m_restricted, phip[1:2, 1:dimbeta])

    @constraint(m_restricted, phiboundneg[k = 1:2, j = 1:dimbeta],
                phi1[k] >= -1*phip[k, j])
    @constraint(m_restricted, phiboundpos[k = 1:2, j = 1:dimbeta],
                phi1[k] >= phip[k, j])
    @constraint(m_restricted, phi0bound[k = 1:2], phi0 >= phi1[k])

    # Ax = b
    @constraint(m_restricted, ax_eq_b[j=1:dimbeta],
                sum(A[j,jj]*x[jj] for jj = 1:dimx) == b[j])

    # the shape and target components of b are constrained
    @constraint(m_restricted, bshp[j=1:length(lpm.beta_shp)],
                b[dimbetaobs + j] == lpm.beta_shp[j])
    # Not coded for A_tgt to have multiple components
    # Although there's nothing conceptual preventing this, just laziness
    @assert (size(lpm.A_tgt)[1] == 1)
    @constraint(m_restricted, b[dimbeta] == testpt)

    if flagrange
        M1 = -1 * A' * OmegaHatI
        M2 = A' * A
        for i = 1:dimx
            @constraint(m_restricted, [k = 1:2],
                           sum(M1[i,j] * phip[k, j] for j = 1:dimbeta)
                        +  sum(M2[i,ii] * phid[k, ii] for ii = 1:dimx)
                        +  rootn * sum(A'[i,j]*b[j] for j in 1:dimbeta)
                        == rootn * dot(A'[i,:], betahat_star))
        end
    else
        for j = 1:dimbeta
            @constraint(m_restricted, [k = 1:2],
                           sum(-1*OmegaHatI[j,jj]*phip[k, jj]
                               for jj = 1:dimbeta)
                        +  sum(A[j,i]*phid[k, i] for i = 1:dimx)
                        +  rootn * b[j]
                        == rootn * betahat[j])
        end
    end
    @objective(m_restricted, Min, phi0)

    optimize!(m_restricted)
    id = baseid * ": betahat restricted, t = $testpt, sample"
    if !checkandresolve(m_restricted, id = id)
        throw(ConeError("$id did not solve correctly."))
    end
    betahatr = value.(b)

    ############################################################################
    # Define cone problem
    ############################################################################
    m_cone = Model(with_optimizer(Gurobi.Optimizer, GLOBAL_GENV[]))
    dp_cone = Dict("OutputFlag" => 0,
                   "Method" => -1, # Gurobi chooses automatically
                   "NumericFocus" => 0,
                   "ScaleFlag" => -1)
    setparams(m_cone, dp_cone)

    @variable(m_cone, sigmapos[1:dimbeta] >= 0)
    @variable(m_cone, sigmaneg[1:dimbeta] >= 0)
    @constraint(m_cone, sum(  sigmapos[j]
                            + sigmaneg[j] for j in 1:dimbeta) <= 1)
    @variable(m_cone, s[1:dimbeta])
    @constraint(m_cone, absvalue[j=1:dimbeta],
                sigmapos[j] - sigmaneg[j] -
                sum(OmegaHatI[j, jj]*s[jj] for jj = 1:dimbeta) == 0)
    @constraint(m_cone, negative[j=1:dimx],
                sum(A'[j,jj]*s[jj] for jj in 1:dimbeta) <= 0)

    if flagrange
        @variable(m_cone, x[1:dimx])
        @constraint(m_cone, rangea[j=1:dimbeta],
                    sum(A[j, i]*x[i] for i = 1:dimx) == s[j])
    end

    ############################################################################
    # Generic solving function for cone problem
    ############################################################################
    function solve_cone(c::Array{Float64, 1}, id)
        @assert length(c) == dimbeta
        s = m_cone[:s] # make sure s corresponds to the right model
        @objective(m_cone, Max, sum(s[i]*c[i] for i in 1:dimbeta))
        optimize!(m_cone)
        ok = checkandresolve(m_cone, id = id, okcodes = [MOI.OPTIMAL],
                             defaultparams = dp_cone)
        if ok
            solution = objective_value(m_cone)
        else
            # Gurobi sometimes throws an error from objective_value is the
            # solve was a failure
            solution = NaN
        end
        return solution, ok
    end

    ############################################################################
    # Test statistic
    ############################################################################
    if flagrange
        invroot_betahatobs_vc = inv(sqrt(betahatobs_vc))
        V = rootn * invroot_betahatobs_vc * (betahatobs - betahatobs_star)
        ts_range = maximum(abs.(V))
    else
        ts_range = 0
    end

    id = baseid * ": cone, t = $testpt, sample (test statistic)"
    ts_cone, ok = solve_cone(betahat_star, id)
    if !ok
        throw(ConeError("$id did not solve correctly."))
    end
    ts_cone = rootn * ts_cone

    ############################################################################
    # Set up lambdas, including solving for automatic one
    ############################################################################
    lambdalist = fill(NaN, length(p.lambda))
    f_autolambda = 0
    for (i, l) in enumerate(p.lambda)
        if l == "auto"
            q_lambda = fill(NaN, p.bsreps)
            for b in 1:p.bsreps
                id = baseid * ": auto lambda, t = $testpt, bootstrap $b"
                c_cone = betahat_star_bs[:, b] - betahat_star
                r, ok = solve_cone(c_cone, id)
                if ok
                    q_lambda[b] = rootn * r
                    b += 1
                else
                    f_autolambda += 1
                end
            end
            q_lambda = filter(!isnan, q_lambda) # only keep the ok solves
            if length(q_lambda) < Int(round(p.bsreps * (1 - p.bswiggle)))
                @error "Too many failures in Auto lambda t = $testpt."
            end
            alphan = 1. / sqrt(log(log(n)))
            lambdalist[i] = min(1 ./ quantileinf(q_lambda, 1 - alphan), 1)
        elseif l == "mult"
            lambdalist[i] = optimal_lambda1(n, dimbeta)
        elseif l == "add"
            lambdalist[i] = optimal_lambda2(n, dimbeta)
        else
            lambdalist[i] = l
        end
    end

    ############################################################################
    # Critical value
    ############################################################################
    # Cone components, depends on lambda
    ts_bs_cone = [fill(NaN, p.bsreps) for _ in 1:length(lambdalist)]
    f_ts_bs_cone = fill(0, length(lambdalist))
    for (k,l) in enumerate(lambdalist)
        for b in 1:p.bsreps
            id = baseid * ": cone, t = $testpt,"*
                          " lambda = $(p.lambda[k]), bootstrap $b"
            c_cone = betahat_star_bs[:, b] - betahat_star + l * betahatr
            r, ok = solve_cone(c_cone, id)
            if ok
                ts_bs_cone[k][b] = rootn * r
            else
                f_ts_bs_cone[k] += 1
            end
        end
    end

    # Range component, does not depend on lambda
    ts_bs_range = fill(NaN, p.bsreps)
    if flagrange
        for b in 1:p.bsreps
            betaobs_diff = (betahatobs_bs[:, b] - betahatobs)
            betastar_diff = betahatobs_star_bs[:, b] - betahatobs_star
            betadiff = betaobs_diff - betastar_diff

            V = rootn * invroot_betahatobs_vc * betadiff
            ts_bs_range[b] = maximum(abs.(V))
        end
    else
        ts_bs_range .= -Inf
    end

    results = create_results_df(p.lambda)
    results.t .= testpt
    results.lambda = p.lambda
    results.r_ts_cone .= ts_cone
    results.r_ts_range .= ts_range
    ts = max(ts_range, ts_cone)
    results.r_ts .= ts
    for k = 1:length(lambdalist)
        results[k, :r_lambda] = lambdalist[k]

        # Only use the good bootstrap draws
        okidx = findall(!isnan, ts_bs_cone[k])
        ts_bs = max.(ts_bs_range[okidx], ts_bs_cone[k][okidx])
        if length(okidx) < Int(round(p.bsreps * (1 - p.bswiggle)))
            @error "Too many failures in Cone t = $testpt,"*
                   " lambda = $(p.lambda[k])."
        end

        results[k, :r_pval] = findpvalue(ts, ts_bs)
        results[k, :r_cv10] = quantileinf(ts_bs, .90)
        results[k, :r_cv05] = quantileinf(ts_bs, .95)
        results[k, :r_cv01] = quantileinf(ts_bs, .99)
        results[k, :r_f_ts_bs_cone] = f_ts_bs_cone[k]
    end
    results.r_f_redraw .= f_redraw
    results.r_f_betahatobs_star_bs .= f_betahatobs_star_bs
    results.r_f_autolambda .= f_autolambda

    return results
end

function create_results_df(lambda_in)
    nrow = length(lambda_in)
    results = DataFrame(t = fill(NaN, nrow),
                        lambda = lambda_in,
                        r_pval = fill(NaN, nrow),
                        r_lambda = fill(NaN, nrow),
                        r_ts = fill(NaN, nrow),
                        r_ts_cone = fill(NaN, nrow),
                        r_ts_range = fill(NaN, nrow),
                        r_cv10 = fill(NaN, nrow),
                        r_cv05 = fill(NaN, nrow),
                        r_cv01 = fill(NaN, nrow),
                        r_f_redraw = fill(NaN, nrow),
                        r_f_betahatobs_star_bs = fill(NaN, nrow),
                        r_f_autolambda = fill(NaN, nrow),
                        r_f_ts_bs_cone = fill(NaN, nrow))
    return results
end

function safe_sqrt(M::Matrix{Float64}, imgtol::Float64)
    # declare to be symmetric -- can make a difference
    M = real(Symmetric(M))
    M = sqrt(M)
    maxim = maximum(imag(M))
    if maxim > imgtol
        @warn "Studentizing matrix has a non-trivial imaginary"*
              " component of size $(maxim)."
    end
    M = real(M) # get rid of small imaginary parts
    return M
end

function optimal_lambda1(n::Integer, p::Integer)
    d1 = log(max(exp(1), p))
    d2 = log(max(exp(1), log(max(exp(1), n))))
    d = sqrt(d1*d2)
    lambda = 1/d
    @assert lambda <= 1
    @assert lambda > 0
    return(lambda)
end

function optimal_lambda2(n::Integer, p::Integer)
    d1 = sqrt(log(max(exp(1), p)))
    d2 = sqrt(2*log(2*log(max(exp(1), n))))
    lambda = 1/(d1 + d2)
    @assert lambda <= 1
    @assert lambda > 0
    return(lambda)
end
