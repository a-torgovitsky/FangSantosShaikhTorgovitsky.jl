load_data <- function(fn = "results.csv") {
    r <- read_csv(fn, col_types = cols())
}

# compute rejection probabilities based on p-values
compute_rejectpr <- function(pvals, levels) {
    rejectprs <- rep(NA, length(levels))
    for (i in seq_along(rejectprs)) {
        rejectprs[i] <- mean(pvals <= levels[i])
    }
    tibble(level = levels, rejectpr = rejectprs)
}

# "loc" codes t as lb/ub of identified set (or NA if neither)
create_loc <- function(r, digits = 3, tol = 1e-3) {
    r$t <- round(r$t, digits = digits)
    loc <- rep(NA, length(r$t))
    loc[(r$t > r$lb - tol) & (r$t < r$lb + tol)] <- "lb"
    loc[(r$t > r$ub - tol) & (r$t < r$ub + tol)] <- "ub"
    # conditional to avoid warning on power graphs where lb/ub not included
    if (nlevels(loc) > 0) fct_relevel(loc, "lb", "ub")
    r$loc <- loc
    return(r)
}

merge_designinfo <- function(r) {
    dgpi <- read_csv("designinfo.csv", col_types = cols()) %>%
        mutate(pointid = (lb == ub))
    r <- left_join(r, dgpi, by = "design")
    r <- create_loc(r)
    r <- r %>% arrange(d, p) %>% mutate(dgp = create_dgp(d,p))
    return(r)
}

rename_proc <- function(r) {
    r <- r %>%
        mutate(lambda = replace_na(lambda, "none")) %>%
        mutate(proc = if_else(lambda == "mult", "cone-mult", proc)) %>%
        mutate(proc =
               recode(proc,
                      "cone" = "FSST",
                      "cone-mult" = "FSST (RoT)",
                      "gmm bs" = "BS Wald",
                      "gmm bs rc" = "BS Wald (RC)"))
}

create_dgp <- function(d, p) {
    str_c("$d = ", d, ", p = ", p, "$")
}

# Turn replication-by-replication results into
#      rejection rate-by-simulation results
clean_data <- function(levels = c(.05, .01)) {
    r <- load_data() %>%
        select(design, rep, proc, t, lambda, r_pval, r_lambda) %>%
        group_by_at(c("design", "proc", "t", "lambda")) %>%
        summarize(compute_rejectpr(r_pval, levels), .groups = "keep") %>%
        ungroup()
    r <- merge_designinfo(r)
    r <- rename_proc(r)
    return(r)
}
