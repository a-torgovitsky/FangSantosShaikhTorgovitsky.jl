# Default design properties
const SAMPLESIZELIST = [1000, 2000, 4000, 8000]
const BSREPS = 250
const MCREPS = 5000
const MCREPS_SMALL = 1000

# Default DGP parameters
const WEVAL = [1.0, 1.0]
const EEVAL = -1.0
const WBDS = (0.0, 2.0)
const VBDS = [(0.5, 1.0), (-3.0, 0.0)]
const TP = MixedLogitTP(weval = WEVAL, eeval = EEVAL)
