# Figure 5
plot_power <- function(level = .10) {
    r <- clean_data(levels = level)
    autoname = "$\\lambda_{n}^{\\rm b}$"
    r <- r %>% mutate(lambda = recode(lambda,
                                      "mult" = "$\\lambda_{n}^{\\rm r}$",
                                      "auto" = autoname))

    fn <- paste0("power_level", round(100*level), ".tex")
    do.call(tikz, c(file = fn, TIKZOPTIONS, width = 7, height = 4.5))

    p <- ggplot(r, aes(x = t, y = rejectpr, color = as.factor(lambda))) +
        labs(x = "Null hypothesis ($\\gamma$)",
             y = "Rejection probability",
             color = "$\\lambda$") +
        geom_line(lwd = LWD) +
        geom_vline(aes(xintercept = lb),
                   linetype = "dotted", color = "black") +
        geom_vline(aes(xintercept = ub),
                   linetype = "dotted", color = "black") +
        geom_hline(yintercept = level,
                   linetype = "dotted", color = "black") +
        scale_y_continuous(breaks = seq(0, 1, by = .10),
                           labels = numform::ff_num(zero = 0, digits = 2)) +
        scale_x_continuous(labels = numform::ff_num(zero = 0, digits = 2)) +
        scale_color_grey(start = 0, end = .7) +
        guides(color = guide_legend(ncol = length(unique(r$lambda)))) +
        facet_grid(samplesize ~ factor(dgp, levels = unique(dgp)),
                   labeller = sslabeller) +
        theme_alex() +
        theme(legend.position = "bottom",
              legend.direction = "vertical")

    print(p)
    dev.off()
    tools::texi2pdf(fn, clean = TRUE)
}
