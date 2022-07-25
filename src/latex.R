make_tableviewer <- function(fn, compile = TRUE,
                             fontsize = "normalsize", caption = "") {
    fntb <- file.path(dirname(fn), paste0("view-", basename(fn)))

    write("\\documentclass[letterpaper, 12pt]{article}", fntb)
    awrite("\\usepackage{amsmath, amsthm, amsfonts, amssymb, bm}", fntb)
    awrite("\\usepackage{fullpage,booktabs,multirow}", fntb)
    awrite("\\thispagestyle{empty}\n", fntb)
    awrite("\\begin{document}\n", fntb)
    awrite("\\begin{table}[t]", fntb)
    awrite("\t\\centering", fntb)
    awrite(paste0("\t\\", fontsize), fntb)
    awrite(paste0("\t\\input{", basename(fn), "}"), fntb)
    awrite(paste0("\t\\caption{", caption, "}"), fntb)
    awrite("\\end{table}", fntb)
    awrite("\\end{document}", fntb)

    if (compile) {
        tools::texi2dvi(fntb, pdf = TRUE, clean = TRUE)
    }
}

awrite <- function (stuff, fn) {
    write(stuff, fn, append = TRUE)
}

start_table <- function(fn, spec) {
    write(paste0("\\begin{tabular}{@{}", spec, "@{}}"),
          fn, append = FALSE)
}

end_table <- function(fn) {
    awrite("\\end{tabular}", fn)
}

insert_toprule <- function(fn) {
    awrite("\t\\toprule", fn)
}

insert_bottomrule <- function(fn) {
    awrite("\t\\bottomrule", fn)
}

insert_midrule <- function(fn) {
    awrite("\t\\midrule", fn)
}

insert_cmidrule <- function(fn, start, stop, spec="") {
    s <- "\t\\cmidrule"
    if (spec != '') {
        s <- paste0(s, "(", spec, ")")
    }
    s <- paste0(s, "{", as.integer(start), "-", as.integer(stop), "}")
    awrite(s, fn)
}

write_row <- function(fn, r) {
    s <- ""
    for (i in seq_along(r)) {
        if (i == 1) {
            s <- "\t"
        } else {
            s <- paste0(s, " & ")
        }
        s <- paste0(s, r[i])
    }
    s <- paste0(s, " \\\\")
    awrite(s, fn)
}
