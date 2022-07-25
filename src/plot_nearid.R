# Figure 4: Null rejection for nearly point-identified designs
plot_nearid <- function(fn = "nearid.tex") {
    r <- clean_data(levels = seq(from = 0, to = .2, by = .001)) %>%
        filter(p >= (d+2), pointid | (!pointid & (loc == "lb"))) %>%
        arrange(d, p, samplesize)

    do.call(tikz, c(file = fn, TIKZOPTIONS, width = 7, height = 4.5))

    p <- ggplot(r, aes(x = level, y = rejectpr, color = as.factor(proc))) +
        labs(x = "Nominal level",
             y = "Rejection probability",
             color = "\\textbf{Test}") +
        geom_line(lwd = LWD) +
        geom_abline(intercept=0, slope=1, lwd = LWD,
                    color="black", linetype="dotted") +
        scale_y_continuous(breaks =
                           seq(0, max(max(r$rejectpr), max(r$level)),
                               by = .05),
                           labels = numform::ff_num(zero = 0, digits = 2)) +
        scale_x_continuous(labels = numform::ff_num(zero = 0, digits = 2)) +
        scale_color_grey(start = 0, end = .7) +
        facet_grid(samplesize ~ factor(dgp, levels = unique(dgp)),
                   labeller = sslabeller) +
        theme_alex() +
        guides(color = guide_legend(ncol = 4)) +
        theme(legend.position = "bottom",
              legend.direction = "vertical")
    print(p)

    dev.off()
    tools::texi2pdf(fn, clean = TRUE)
}

