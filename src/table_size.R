# Formats large numbers as squares to keep the column widths roughly even
square_string <- function(s, large = 10000) {
    s <- as.integer(s)
    snew <- rep(NA_character_, length(s))
    for (i in seq_along(s)) {
        ss <- as.integer(s[i])
        if (ss >= large) {
            snew[i] <- paste0("$", as.character(round(sqrt(ss))), "^{2}$")
        }
        else {
            snew[i] <- s[i]
        }
    }
    return(snew)
}

tables_size <- function() {
    r <- clean_data(levels = c(.01, .05, .10))
    r %>%
        group_by(loc, level) %>%
        group_split %>%
        lapply(table_size_one)
    return(invisible(NULL))
}

table_size_one <- function(df) {
    loc <- df$loc[1]
    level <- df$level[1]

    df <- df %>%
        mutate(rejectpr = numform::f_num(rejectpr, digits = 3)) %>%
        select(samplesize, p, d, rejectpr) %>%
        arrange(samplesize, p, d) %>%
        spread(d, rejectpr, fill = "--")

    fn <- paste0("size_table_", loc, "_", as.integer(100*level), ".tex")
    cols <- ncol(df)
    start_table(fn, paste("c", "c", rep("c", cols - 2), collapse = ''))
    insert_toprule(fn)
    r <- c("", "",
           paste0("\\multicolumn{", cols - 2, "}{c}{$d$}"))
    write_row(fn, r)
    insert_cmidrule(fn, 3, cols, "lr")
    r <- c("$n$", "$p$", square_string(names(df)[-(1:2)]))
    write_row(fn, r)
    insert_midrule(fn)

    n_unique <- unique(df$samplesize)
    r <- rep("", cols)
    for (i in seq_along(n_unique)) {
        df2 <- df %>% filter(samplesize == n_unique[i])
        r[1] <- paste0("\\multirow{", nrow(df2), "}{*}{", n_unique[i], "}")

        for (j in seq_len(nrow(df2))) {
            if (j > 1) {
                r[1] <- "" # don't write multirow twice
            }
            r[2] <- as.character(df2$p[j])
            r[-(1:2)] <- df2[j, -(1:2)]
            write_row(fn, r)
        }
        if (i < length(n_unique)) {
            insert_cmidrule(fn, 1, cols, "lr")
        }
    }
    insert_bottomrule(fn)
    end_table(fn)
    make_tableviewer(fn, fontsize = "small",
                     caption = paste0("Rejection rates for nominal ", level,
                                      " test at ", loc, "."))
}
