@with_kw struct MixedLogitTP
    weval::AbstractArray{<:Real, 1} = [1.0, 1.0]
    eeval::Real = -1.0
    printname::String = ""

    function MixedLogitTP(weval, eeval, printname)
        printname = "df_elast=$eeval"
        new(weval, eeval, printname)
    end
end

struct MixedLogitDGP
    vsupp::AbstractArray{<:Real, 2}
    vpr::AbstractArray{<:Real, 1}
    wsupp::AbstractArray{<:Real, 2}
    wpr::AbstractArray{<:Real, 1}
    tp::MixedLogitTP

    function MixedLogitDGP(vsupp, vpr, wsupp, wpr, tp)
        @assert size(vsupp, 1) == length(vpr)
        @assert size(wsupp, 1) == length(wpr)
        @assert size(vsupp, 2) == size(wsupp, 2)
        @assert isapprox(1, sum(vpr))
        @assert isapprox(1, sum(wpr))
        @assert length(tp.weval) == size(wsupp, 2)

        new(vsupp, vpr, wsupp, wpr, tp)
    end
end

# Keep everything of the DGP but change the target parameter
function MixedLogitDGP(dgp::MixedLogitDGP, tp::MixedLogitTP)
    return MixedLogitDGP(dgp.vsupp, dgp.vpr, dgp.wsupp, dgp.wpr, tp)
end

function MixedLogitDGP(L::Integer, R::Integer;
                       wbds = (0.0, 3.0),
                       vbds = [(0.0, 0.5), (-3.0, 0.0)],
                       tp = MixedLogitTP())

    vsupp0 = generate_sobol_grid(Int.(ceil(sqrt(R))), vbds[1])
    vsupp1 = generate_sobol_grid(Int.(ceil(sqrt(R))), vbds[2])
    vsupp = create_cartesian(vsupp0, vsupp1)
    vpr = fill(1 ./ size(vsupp, 1), size(vsupp, 1))

    wsupp1 = collect(range(wbds[1], stop = wbds[2], length = L))
    wsupp = hcat(ones(L), wsupp1)
    wpr = fill(1 ./ size(wsupp, 1), size(wsupp, 1)) # pr w is uniform

    return MixedLogitDGP(vsupp, vpr, wsupp, wpr, tp)
end

function generate_sobol_grid(N::Integer, bounds::Tuple{<:Real, <:Real})
    grid = fill(NaN, N)
    s = SobolSeq(1)
    for i = 1:N
        grid[i] = next!(s)[1]
    end
    grid = bounds[1] .+ grid .* (bounds[2] - bounds[1])
    return grid
end

function create_cartesian(w1::AbstractArray{<:Real, 1},
                          w2::AbstractArray{<:Real, 1})
    return create_cartesian(w1[:,:], w2)
end

function create_cartesian(w1::AbstractArray{<:Real, 2},
                          w2::AbstractArray{<:Real, 1})

    w = NaN*ones(size(w1, 1)*length(w2), size(w1, 2) + 1)
    counter = 1
    for r1 in eachrow(w1)
        for r2 in w2
            w[counter, :] = vcat(r1, r2)
            counter += 1
        end
    end

    return w
end

# y = 1(w'v >= u)
# with u ~ logistic
function drawsample(dgp::MixedLogitDGP, n::Integer = 500, seed::Integer = 1)
    Random.seed!(seed)

    dv = Categorical(dgp.vpr)
    iv = rand(dv, n)
    v = dgp.vsupp[iv,:]

    dw = Categorical(dgp.wpr)
    iw = rand(dw, n)
    w = dgp.wsupp[iw, :]

    du = Logistic()
    u = rand(du, n)

    y = -1*ones(n)
    for i = 1:n
        y[i] = compute_choice(w[i,:], v[i,:], u[i])
    end

    names = append!([:y], [Symbol("w$k") for k in 1:size(w,2)])
    data = DataFrame(hcat(y, w), names)

    return data
end

function drawsample_safe(dgp::MixedLogitDGP, n::Integer = 500,
                         seed::Integer = 1; maxfails::Integer = 10,
                         scramble = 4529023023)
    fails = 0
    while true
        data = drawsample(dgp, n, seed + fails * scramble)
        bobs, vc = beta_obs_np(data)
        # only the size of "bold" matters for the comparison, so fill with 1's
        if checkresample(bobs, fill(1.0, size(dgp.wsupp)[1]), vc)
            return data, fails
        else
            fails += 1
        end
        if (fails >= maxfails)
            error("Draw fails exceeded maxfails. Something is wrong!")
        end
    end
end

function checkresample(bnew::Array{<:Real, 1}, bold::Array{<:Real, 1},
                       vcnew::Union{Array{<:Real, 2}, Nothing} = nothing)
    length(bnew) == length(bold) ? ok = true : ok = false
    if ~isnothing(vcnew)
        if (sum(diag(vcnew) .> 0) == length(bold))
            ok = true * ok
        else
            ok = false
        end
    end
    return ok
end

function compute_choice(w::AbstractArray{<:Real, 1},
                        v::AbstractArray{<:Real, 1},
                        u::Real)
    w = compute_index(w, v)
    if w >= u
        return 1
    else
        return 0
    end
end

function compute_true_choice_prob(dgp::MixedLogitDGP)
    return [compute_choice_prob(w, dgp.vsupp, dgp.vpr)
                    for w in eachrow(dgp.wsupp)]
end

function compute_true_bounds(dgp::MixedLogitDGP,
                             tp::MixedLogitTP)
    truechoicepr = compute_true_choice_prob(dgp)
    return compute_bounds(MixedLogitDGP(dgp, tp), truechoicepr)
end

function compute_true_bounds(dgp::MixedLogitDGP,
                             tplist::Array{MixedLogitTP, 1})
    r = DataFrame(target = fill("", length(tplist)),
                  lb = fill(NaN, length(tplist)),
                  ub = fill(NaN, length(tplist)))
    for (i, tp) in enumerate(tplist)
        bd, _ = compute_true_bounds(dgp, tp)
        r[i,:target] = tp.printname
        r[i,:lb] = bd[1]
        r[i,:ub] = bd[2]
    end
    return r
end

function compute_true_bounds(dgp::MixedLogitDGP)
    compute_true_bounds(dgp, dgp.tp)
end

function compute_bounds(dgp::MixedLogitDGP, choicepr::AbstractArray{<:Real, 1})
    lpm = dgp_to_lpm(dgp)
    bounds, _ = compute_bounds(lpm, choicepr)
    return bounds
end

function true_tp_value(dgp::MixedLogitDGP)
    lpm = dgp_to_lpm(dgp)
    return lpm.A_tgt*dgp.vpr
end

function true_tp_value(dgp::MixedLogitDGP, tp::MixedLogitTP)
    return true_tp_value(MixedLogitDGP(dgp, tp))
end

function dgp_to_lpm(dgp::MixedLogitDGP; gentrue::Bool = false)

    if gentrue
        bobs = (d -> beta_true(dgp))
    else
        bobs = (d -> beta_obs_np(d))
    end

    lpm = LPModel(  A_obs = A_obs(dgp.wsupp, dgp.vsupp),
                    beta_obs = bobs,
                    A_shp = A_shp(dgp.vsupp),
                    beta_shp = beta_shp(),
                    A_tgt = A_tgt(dgp.vsupp, dgp.tp))

    return lpm
end

function A_obs(wsupp::AbstractArray{<:Real, 2},
               vsupp::AbstractArray{<:Real, 2})
    A_obs = zeros(size(wsupp, 1), size(vsupp, 1))

    for i in 1:size(wsupp,1)
        for j in 1:size(vsupp, 1)
            A_obs[i,j] = compute_type_choice_prob(wsupp[i,:], vsupp[j,:])
        end
    end

    return A_obs
end

function beta_obs_np(data::DataFrame)
    df = groupby(data, [:w1, :w2])
    df = combine(df, (mean = :y => mean),
                     (avar = :y => (x -> var(x, corrected = false))),
                     (nk = :y => length))
    df = sort(df, [:w1, :w2]) # always return in the same order
    m = df.mean
    v = diagm(0 => nrow(data) .* (df.avar ./ df.nk)) # asymptotic variance
    return m, v
end

function beta_true(dgp::MixedLogitDGP)
    btrue = A_obs(dgp.wsupp, dgp.vsupp) * dgp.vpr
    bvar = Matrix{Float64}(I, length(btrue), length(btrue))
    return btrue, bvar
end

function A_shp(vsupp::AbstractArray{<:Real, 2})
    return ones(1, size(vsupp, 1))
end

function beta_shp()
    return ones(1)
end

function A_tgt(vsupp::AbstractArray{<:Real, 2}, tp::MixedLogitTP)
    elasts = [compute_deriv_type_elasticity(tp.weval, v)
              for v in eachrow(vsupp)]
    A = [Int(e <= tp.eeval) for e in elasts]
    return A'
end

function compute_choice_prob(w::AbstractArray{<:Real, 1},
                             vsupp::AbstractArray{<:Real, 2},
                             vpr::AbstractArray{<:Real, 1})
    pgv = [compute_type_choice_prob(w, v) for v in eachrow(vsupp)]
    return dot(pgv, vpr)
end

function compute_index(w::AbstractArray{<:Real, 1},
                       v::AbstractArray{<:Real, 1})
    @assert length(w) == length(v)
    return dot(w, v)
end

function compute_type_choice_prob(w::AbstractArray{<:Real, 1},
                                  v::AbstractArray{<:Real, 1})
    i = compute_index(w, v)
    return logistic(i)
end
logistic(z::Real) = inv(exp(-z) + one(z))

function compute_deriv_type_choice_prob(w::AbstractArray{<:Real, 1},
                                        v::AbstractArray{<:Real, 1})
    i = compute_index(w, v)
    dl = dlogistic(i)
    return dl*v[2]
end

function dlogistic(z::Real)
    l = logistic(z)
    return l*(1-l) # see pg. 58 of Train 2009
end

function compute_deriv_type_elasticity(w::AbstractArray{<:Real, 1},
                                       v::AbstractArray{<:Real, 1})
    pr = compute_type_choice_prob(w,v)
    dpr = compute_deriv_type_choice_prob(w,v)
    return (dpr * (w[2]/pr))
end

function compute_change_type_choice_prob(wfrom::AbstractArray{<:Real, 1},
                                         wto::AbstractArray{<:Real, 1},
                                         v::AbstractArray{<:Real, 1})
    prfrom = compute_type_choice_prob(wfrom, v)
    prto = compute_type_choice_prob(wto, v)
    return (prto - prfrom)
end
