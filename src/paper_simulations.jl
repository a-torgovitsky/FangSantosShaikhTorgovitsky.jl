# Figure 3
function idsets(resultsdir)
    dgplist = [MixedLogitDGP(l, r, wbds = WBDS, vbds = VBDS, tp = TP)
               for l in [4, 16], r in [16, 100, 400, 1600]]
    weval = WEVAL[2]
    elist = 0:-.02:-3

    longlist = Iterators.product(dgplist, elist)
    dl = length(longlist)
    df = DataFrame(dimx = fill(NaN, dl),
                   dimbeta = fill(NaN, dl),
                   weval = fill(NaN, dl),
                   eeval = fill(NaN, dl),
                   eevalstar = fill(dgplist[1].tp.eeval, dl),
                   lb = fill(NaN, dl),
                   ub = fill(NaN, dl),
                   truth = fill(NaN, dl))

    for (i, l) in enumerate(longlist)
        dgp = l[1]
        tp = MixedLogitTP(weval = [1.0, weval], eeval = l[2])
        bds = compute_true_bounds(dgp, tp)

        df[i, :dimx] = size(dgp.vsupp)[1]
        df[i, :dimbeta] = size(dgp.wsupp)[1] + 2
        df[i, :weval] = weval
        df[i, :eeval] = l[2]
        df[i, :lb] = bds[1]
        df[i, :ub] = bds[2]
        df[i, :truth] = true_tp_value(dgp, tp)
    end
    fn_results = joinpath(resultsdir, "results.csv")
    CSV.write(fn_results, df, writeheader = true)
    return fn_results
end
export idsets

# Figure 4
function simnearid(resultsdir; bsreps = BSREPS, mcreps = MCREPS)
    dgplist = [MixedLogitDGP(l, r, wbds = WBDS, vbds = VBDS, tp = TP)
               for l in [4, 16], r in [4, 16] if (l >= r)]

    designlist = [Design(dgp = d,
                         samplesize = n,
                         bsreps = bsreps,
                         mcreps = mcreps,
                         testptlist = [],
                         cone_lambdalist = ["auto", "mult"],
                         gmm = true)
                  for d in dgplist, n in SAMPLESIZELIST[2:3]][:]

    fn_results = rundesigns(resultsdir, designlist)
    return fn_results
end
export simnearid

# Table 1
function simsize(resultsdir; bsreps = BSREPS, mcreps = MCREPS,
                 mcreps_small = MCREPS_SMALL, cone_lambda = ["auto"])
    llists = [[4, 16],
              [4, 16, 6^2],
              [4, 16, 6^2, 7^2],
              [4, 16, 6^2, 7^2, 9^2]] # for different sample sizes
    rlist = [100, 400, 1600, 4900, 10_000, 225^2, 317^2]

    dgplist = [MixedLogitDGP(l, r, wbds = WBDS, vbds = VBDS, tp = TP)
               for l in sort(unique(vcat(llists...))), r in rlist]

    designlist = [Design(dgp = d,
                         samplesize = SAMPLESIZELIST[i],
                         bsreps = bsreps,
                         mcreps = (size(d.vsupp)[1] <= 10_000 ?
                                   mcreps : mcreps_small),
                         testptlist = [],
                         cone_lambdalist = cone_lambda,
                         gmm = false)
                  for i in 1:length(SAMPLESIZELIST),
                      d in dgplist if (size(d.wsupp)[1] in llists[i])][:]
    # do the really time-consuming designs first
    designlist = sort(designlist,
                      by = (d -> size(d.dgp.wsupp)[1] * size(d.dgp.vsupp)[1]),
                      rev = true)
    fn_results = rundesigns(resultsdir, designlist)
    return fn_results
end
export simsize

# Figures 5 and 6
function simpower(resultsdir; bsreps = BSREPS, mcreps = MCREPS, step = .02)
    dgplist = [MixedLogitDGP(4, r, wbds = WBDS, vbds = VBDS, tp = TP)
               for r in [16, 1600]]
    designlist = [Design(dgp = d,
                         samplesize = n,
                         bsreps = bsreps,
                         mcreps = mcreps,
                         testptlist = collect(range(0, 1, step = step)),
                         cone_lambdalist = [0, "auto", "mult"],
                         gmm = false)
                  for d in dgplist, n in SAMPLESIZELIST[2:3]][:]
    fn_results = rundesigns(resultsdir, designlist)
    return fn_results
end
export simpower
