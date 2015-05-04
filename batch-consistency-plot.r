# 1-20-10 Qunhua Li This program first plots correspondence curve and IDR threshold plot (i.e. number of selected peaks vs
# IDR) for each pair of sample usage: Rscript batch-consistency-plot-merged.r [npairs] [output.dir] [input.file.prefix 1, 2, 3
# ...] [npairs]: integer, number of consistency analyses (e.g. if 2 replicates, npairs=1, if 3 replicates, npairs=3
# [output.dir]: output directory for plot [input.file.prefix 1, 2, 3]: prefix for the output from batch-consistency-analysis2.
# They are the input files for merged analysis see below for examples (i.e. saved.file.prefix). It can be multiple files

# args <- commandArgs(trailingOnly=T)


lines = c("MCF10A", "MCF7", "MDA231")
mods = c("H3K4AC", "H3K4ME3")
prefixes = list()
for (l in lines) {
    for (m in mods) {
        prefixes[[length(prefixes) + 1]] = paste(wd, "/IDR_results/", l, "_", m, sep = "")
    }
}
prefixes = unlist(prefixes)

print('plotting IDR curves to IDR_results-plot.pdf')
plot_IDR_results(prefixes, output.file.prefix = "IDR_results")
