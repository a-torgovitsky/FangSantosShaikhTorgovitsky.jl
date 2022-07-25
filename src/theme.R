LWD <- 2
TIKZOPTIONS <- list(pointsize = 12, standAlone = TRUE, verbose = FALSE)
options(tikzLatexPackages = c(getOption("tikzLatexPackages"),
                              "\\usepackage{amsmath}"))

theme_alex <- function () {
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.title.align = 0.5
    )
}

sslabeller = labeller(samplesize = function(s) paste0("$n = ", s, "$"))
