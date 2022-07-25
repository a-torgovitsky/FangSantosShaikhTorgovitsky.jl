function printdiv(c::Char, repeat = 1)
    div = c^80
    for i = 1:repeat
        println(div)
    end
end
export printdiv

function harmonize!(df1::DataFrame, df2::DataFrame)
    addcols!(df1, df2)
    addcols!(df2, df1)

    # put columns of df2 in order of df1
    perm = Array{Int64, 1}(undef, length(names(df2)))
    for (j,n) in enumerate(names(df1))
        perm[j] = findfirst(names(df2) .== n)
    end
    select!(df2, perm)
end

# Add any columns in dfadd but not in df to df
function addcols!(df::DataFrame, dfadd::DataFrame)
    allowmissing!(df)

    for n in names(dfadd)
        if !(n in names(df))
            Tvec = typeof.(dfadd[:, n])
            Tvec = unique(Tvec)
            T = Union{Tvec..., Missing} # make sure type will be ok
            df[:, n] = Array{T, 1}(missing, nrow(df))
        end
    end
end
