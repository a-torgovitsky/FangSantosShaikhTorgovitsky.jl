# 2021-08-03
#
# This is for diagnostic purposes
# Not the prettiest graph because I expect to delete it
plot_errors <- function() {
    load_data() %>%
        select(design, rep, proc, t, lambda, starts_with("r_f_")) %>%
        group_by(design, rep, proc, t, lambda) %>%
        summarize_all(mean) %>%
        pivot_longer(starts_with("r_f_")) %>%
        merge_designinfo %>%
        replace_na(list(lambda = "na", value = 0)) ->
        r

    error_plot <- function(df) {
        p <- ggplot(data = df, aes(x = name, y = value)) +
            geom_col() +
            ylim(0, .25) +
            theme_alex() +
            labs(x = "Error type",
                 y = "Average number of fails in each MC rep",
                 title = paste0(df$proc[1], ", t = ", df$t[1],
                                ", lambda = ", df$lambda[1])) +
            facet_grid(samplesize ~ factor(dgp, levels = unique(dgp)))
        return(p)
    }
    r %>%
        group_by(dgp, t, proc, lambda) %>%
        arrange(.by_group = TRUE) %>%
        group_split %>%
        lapply(error_plot) -> plist
    mp <- ggarrange(plotlist = plist, nrow = 1, ncol = 1)
    ggexport(mp, filename = "errors.pdf")
}
