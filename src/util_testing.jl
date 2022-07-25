# Draws a sample of size M, with replacement
function resampledata(data::DataFrame, m::Integer, seed::Integer;
                      replace::Bool=true)

    Random.seed!(seed)
    n = nrow(data)
    if (m > n) & !replace
        println("""WARNING: Resample size is larger than original size,"
                but replace is false. Setting replace to true to proceed.""")
        replace = true
    end
    idx = sample(1:1:n, m, replace=replace)
    datars = data[idx,:]

    return datars
end

# Compute quantiles as Q(t) = inf\{ q : F(q) \geq t \}
function quantileinf(v::Array{<:Real, 1}, q::Real)
    qq = Float64[q]
    return quantileinf(v, qq)[1]
end

function quantileinf(v::Array{<:Real, 1}, q::Array{<:Real, 1})
    # sort the data vector
    sort!(v)
    v = filter(x -> !isnan(x), v)
    n = length(v)

    r = fill(NaN, length(q)) # default is NaN

    idxinf = findall(x -> (x == 0.), q) # 0 quantile is -Inf
    r[idxinf] .= -Inf

    idx = findall(x -> (x .> 0)&(x .<= 1), q) # quantiles in (0,1] well-defined
    if length(idx) > 0
        for i in idx
            pos = Int(ceil(q[i] * n))
            r[i] = v[pos]
        end
    end

    return r
end

# findpvalue:
# Sort the draws
# Find the largest draw that TS is greater than or equal to
# The empirical CDF at this draw is idx/length(draws)
# Thus, it is the
#     idx/length(draws) quantile
# of the bootstrap distribution
# So, TS is larger than the idx/length(draws) quantile
# But it is smaller than the (idx + 1)/length(draws) quantile
# So, we would reject if
#     1 - alpha = idx/length(draws)
# but we would fail to reject if
#     1 - alpha = idx/length(draws) + epsilon
# This implies
#     alpha = 1 - idx/length(draws)
# is the smallest value of alpha for which we would reject
# and thus is the p-value
function findpvalue(ts::Real, draws::Array{<:Real, 1})
    draws = sort(draws)
    idx = findlast(x -> (ts >= x), draws)
    if isnothing(idx) # ts < minimum(x), so never reject
        pvalue = 1
    else
        pvalue = 1 - idx/length(draws)
    end
end
