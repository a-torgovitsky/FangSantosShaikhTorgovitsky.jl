# Figure 3
plot_idset <- function() {
    r <- load_data()
    fn <- "identified-sets.tex"
    do.call(tikz, c(file = fn, TIKZOPTIONS, width = 7, height = 4.5))

    label_d <- function(s) paste0("$d = ", s, "$")
    eevalstar <- r$eevalstar[1]

    p <- ggplot(r, aes(x = eeval, y = truth)) +
        geom_ribbon(aes(ymin=lb, ymax=ub, fill = as.factor(dimbeta)),
                    alpha = .5) +
        geom_line(lwd = 1.5) +
        labs(x = "Elasticity ($t$)",
             y = "Distribution function (truth and pointwise bounds)") +
        geom_vline(xintercept = eevalstar,
                   linetype = "dotted", color = "black") +
        facet_wrap(~ dimx, labeller = labeller(dimx = label_d)) +
        scale_fill_grey(start = .5, end = 1) +
        theme_alex() +
        theme(legend.position = "none")
    print(p)

    dev.off()
    tools::texi2pdf(fn, clean = TRUE)
}
