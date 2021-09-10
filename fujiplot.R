#### debug options ####

# setwd("/mnt/obi0/mhomilius/projects/fujiplot/")
# fullargs <- c(
#   "/usr/lib/R/bin/exec/R",
#   "--no-echo",
#   "--no-restore",
#   "--file=fujiplot.R",
#   "--args",
#   "/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/circos/input_cell_type_gated_quantile_all.txt",
#   "/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/circos/traitlist_cell_type_gated_quantile_all.txt",
#   "/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/circos/output/cell_type_gated_quantile_all"
# )
# 
# args <- c(
#   "/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/circos/input_cell_type_gated_quantile_all.txt",
#   "/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/circos/traitlist_cell_type_gated_quantile_all.txt",
#   "/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/circos/output/cell_type_gated_quantile_all"
# )

####

fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

#print(fullargs)
#stop()

script_name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))
script_dir <- dirname(script_name)

suppressMessages(library(dplyr))
library(stringr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(forcats)


# constants
.VERSION = "1.0.3"
CIRCOS_CONF = file.path(script_dir, "config",      "circos.conf")
CIRCOS_PATH = "circos"
# CIRCOS_DEBUG_GROUP = "text,textplace"
CIRCOS_DEBUG_GROUP = "summary"
OUTPUT_BARPLOT = TRUE
SCATTER_BACKGROUND_COLOR_ALPHA = 0.3
LARGE_POINT_SIZE = 18
SMALL_POINT_SIZE = 12

# intermediate files
COLOR_CONF =              file.path(script_dir, "config",      "color.conf")
SCATTER_BACKGROUND_CONF = file.path(script_dir, "config",      "scatter_background.conf")
HIGHLIGHT_DATA =          file.path(script_dir, "data_tracks", "highlights.txt")
SCATTER_DATA =            file.path(script_dir, "data_tracks", "scatter.txt")
STACKED_DATA =            file.path(script_dir, "data_tracks", "stacked.txt")
LABEL_DATA =              file.path(script_dir, "data_tracks", "label.txt")

################################################################################
writeLines(c("*********************************************************************",
             "* Fuji plot -- a circos representation of multiple GWAS results",
     sprintf("* Version %s", .VERSION),
             "* Masahiro Kanai (mkanai@g.harvard.edu)",
             "* Harvard Medical School / RIKEN IMS / Osaka Univerisity",
             "* GNU General Public License v3",
             "*********************************************************************"
           ))

if (identical(args, character(0))) {
  args = file.path(script_dir, "input_example", c("input.txt", "traitlist.txt"))
}
input_fname = normalizePath(args[1])
traitlist_fname = normalizePath(args[2])
if (length(args) > 2){
  output_dir = normalizePath(args[3])
} else {
  output_dir = file.path(script_dir, 'output_example')
}

################################################################################
# helper func
most_common = function(x) {tail(names(sort(table(x))), 1)}

################################################################################
# load data
message("Loading input files...")

df = read.table(input_fname, T, sep = '\t', as.is = T, quote = '', comment.char = '')
traitlist = read.table(traitlist_fname, T, sep = '\t', as.is = T, quote = '', comment.char = '', fileEncoding='utf-8')

n_loci = length(unique(df$LOCUS_ID))
writeLines(c(
     sprintf("* Input data: %s", input_fname),
     sprintf("* Number of significant SNPs: %d", nrow(df)),
     sprintf("* Number of unique loci: %d", n_loci),
             "",
     sprintf("* Trait list: %s", traitlist_fname),
     sprintf("* Number of traits: %d (%s)", nrow(traitlist), str_c(traitlist$TRAIT, collapse = ',')),
     sprintf("* Number of categories: %d (%s)", length(unique(traitlist$CATEGORY)), str_c(unique(traitlist$CATEGORY), collapse = ',')),
             "",
     sprintf("* Output dir: %s", output_dir)
           ))

if ( ! dir.exists(output_dir) ){
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

input_traits = traitlist$TRAIT
if (!all(df$TRAIT %in% input_traits)) {
    missing_traits = setdiff(df$TRAIT, input_traits)
    n_missing = length(missing_traits)
    stop(sprintf("TRAIT columns mismatch.\n%d trait%s in %s %s missing from %s (%s).",
                 n_missing, ifelse(n_missing > 1, "s", ""), input_fname,
                 ifelse(n_missing > 1, "are", "is"), traitlist_fname,
                 str_c(missing_traits, collapse=",")))
}

traitlist = traitlist %>% filter(TRAIT %in% df$TRAIT) %>%
                           mutate(idx = 1:n(),
                                  category_lower = str_replace_all(str_to_lower(CATEGORY), '[^a-z0-9_]', '_')) %>%
                           mutate(parameters = str_c('fill_color=', category_lower))
excluded_traits = setdiff(input_traits, traitlist$TRAIT)
writeLines(c("",
     sprintf("Excluded %d traits because of no significant SNPs (%s).", length(excluded_traits), str_c(excluded_traits, collapse = ',')),
             ""
           ))

################################################################################
message("Generating configuration and data files for circos...")

# output color config
str_c_comma = function(x){str_c(x, collapse = ",")}
cols = traitlist %>% select(category_lower, COLOR) %>%
                     unique() %>%
                     mutate(rgb = apply(t(col2rgb(COLOR)), 1, str_c_comma),
                            rgba = apply(floor(t((1 - SCATTER_BACKGROUND_COLOR_ALPHA) * 255 + SCATTER_BACKGROUND_COLOR_ALPHA * col2rgb(COLOR))), 1, str_c_comma))
writeLines(c("<colors>",
             sprintf("\t%s = %s\n\talpha_%s = %s", cols$category_lower, cols$rgb, cols$category_lower, cols$rgba),
             "</colors>"), COLOR_CONF)
message(sprintf("* Color configuration: %s", COLOR_CONF))


################################################################################
# output scatterplot background config
colsep = table(factor(traitlist$category_lower, levels=unique(traitlist$category_lower)))
bg = data.frame(category_lower = unique(traitlist$category_lower),
                y0 = nrow(traitlist) - cumsum(colsep) - 0.5,
                y1 = nrow(traitlist) - cumsum(c(0,colsep))[1:length(colsep)] - 0.5)
writeLines(c("<backgrounds>",
             sprintf("<background>\n\tcolor = alpha_%s\n\ty0 = %.1f\n\ty1 = %.1f\n</background>", bg$category_lower, bg$y0, bg$y1),
             "</backgrounds>"), SCATTER_BACKGROUND_CONF)
message(sprintf("* Scatter background configuration: %s", SCATTER_BACKGROUND_CONF))


################################################################################
# output pleiotropy highlight data
nsnps_per_locus = df %>% group_by(LOCUS_ID) %>% summarize(n = n())
#nsnps_per_locus = df %>% group_by(LOCUS_ID, MARKER) %>% select(MARKER) %>% summarize(n = n_distinct())
df = df %>% mutate(CHR = str_c("hs", CHR),
                   nsnps = nsnps_per_locus$n[match(LOCUS_ID, nsnps_per_locus$LOCUS_ID)])


# CHANGED: showing all loci
# inter_categorical = df %>% group_by(LOCUS_ID) %>% summarize(CHR = most_common(CHR),
#                                                             BP = most_common(BP),
#                                                             n = length(unique(CATEGORY))) %>%
#                                                   filter(n > 1)
inter_categorical = df %>% group_by(LOCUS_ID) %>% summarize(CHR = most_common(CHR),
                                                            BP = most_common(BP),
                                                            n = length(unique(CATEGORY)))

write.table(inter_categorical[c("CHR", "BP", "BP")], HIGHLIGHT_DATA, sep = "\t", row.names = F, col.names = F, quote = F)
message(sprintf("* Highlights data (inter-categorical pleiotropic loci): %s", HIGHLIGHT_DATA))


################################################################################
# output outer scatter plot data
scatter = merge(df, traitlist, by = "TRAIT", all.x = T)
scatter$value = nrow(traitlist) - scatter$idx

# CHANGED use large or small points based on whether there is an annotated gene or intergenic region
#scatter$parameters = str_c(scatter$parameters, str_c('z=', scatter$nsnps), str_c('glyph_size=', ifelse(scatter$nsnps > 1, LARGE_POINT_SIZE, SMALL_POINT_SIZE)), sep = ",")
scatter$parameters = str_c(scatter$parameters, str_c('z=', scatter$nsnps), str_c('glyph_size=', ifelse(!is.na(scatter$GENE.x) & !is.na(scatter$GENE.y), LARGE_POINT_SIZE, SMALL_POINT_SIZE)), sep = ",")

scatter = scatter[order(scatter$nsnps, decreasing=T),]
write.table(scatter[c("CHR", "BP", "BP", "value", "parameters")], SCATTER_DATA, sep = "\t", row.names = F, col.names = F, quote = F)
message(sprintf("* Scatter plot data (significant loci): %s", SCATTER_DATA))


################################################################################

# TODO change or drop the stacked data viz

# output inner stacked scatter plot data
stacked = list()
stacked_y = rep(0, n_loci)
names(stacked_y) = 1:n_loci

for (i in 1:nrow(traitlist)) {
  x = subset(scatter, idx == i)
  x$value = stacked_y[x$LOCUS_ID]
  stacked_y[x$LOCUS_ID] = stacked_y[x$LOCUS_ID] + 1
  stacked[[i]] = x
}

stacked = do.call(rbind, stacked)
stacked$parameters = str_c(stacked$parameters, str_c('z=', stacked$nsnps), sep = ",")
stacked = stacked[order(stacked$nsnps, decreasing=T),]
write.table(stacked[c("CHR", "BP", "BP", "value", "parameters")], STACKED_DATA, sep = "\t", row.names = F, col.names = F, quote = F)
message(sprintf("* Stacked bar plot data (# significant SNPs per locus): %s", STACKED_DATA))

################################################################################

# output label data
# CHANGED: dropping NA labels
# label = df %>% filter(LOCUS_ID %in% inter_categorical$LOCUS_ID) %>%
#                  group_by(LOCUS_ID) %>%
#                    summarize(CHR = most_common(CHR),
#                              BP = most_common(BP),
#                              GENE = most_common(GENE))
label = df %>% filter(LOCUS_ID %in% inter_categorical$LOCUS_ID) %>%
  group_by(LOCUS_ID) %>%
  filter(!is.na(GENE)) %>%
  summarize(CHR = most_common(CHR),
            BP = most_common(BP),
            GENE = most_common(GENE))
write.table(label[c("CHR", "BP", "BP", "GENE")], LABEL_DATA, sep = "\t", row.names = F, col.names = F, quote = F)
message(sprintf("* Label data (name of inter-categorical pleiotropic loci): %s", LABEL_DATA))


################################################################################
# call circos
cmd = sprintf("%s -conf %s %s", CIRCOS_PATH, CIRCOS_CONF, ifelse(CIRCOS_DEBUG_GROUP == "", "", sprintf("-debug_group %s", CIRCOS_DEBUG_GROUP)))
writeLines(c("",
             "Calling circos to plot...",
     sprintf("* Call: %s", cmd),
             ""
           ))
setwd(script_dir) # circos config files are specified with the relative paths
system(cmd)

# move the output file from circos to the specified location
if (output_dir != file.path(script_dir, 'output')){
  for(ext in c('png', 'svg')){
    system(sprintf("mv %s %s", file.path(script_dir, 'output', sprintf('circos.%s', ext)), file.path(output_dir, sprintf('circos.%s', ext))))
  }
}

# clean-up the intermediate files
for(f in c(
  COLOR_CONF,
  SCATTER_BACKGROUND_CONF,
  HIGHLIGHT_DATA,
  SCATTER_DATA,
  STACKED_DATA,
  LABEL_DATA
)){
  if (file.exists(f)) file.remove(f)
}

################################################################################
# output bar plot
if (OUTPUT_BARPLOT) {
  bar = df %>% group_by(TRAIT) %>%
               summarize(total = length(MARKER),
                         pleiotropic = length(MARKER[nsnps > 1]),
                         inter_categorical = length(MARKER[LOCUS_ID %in% inter_categorical$LOCUS_ID])) %>%
               mutate(single = total - pleiotropic,
                      intra_categorical = pleiotropic - inter_categorical)
  bar = merge(bar, traitlist, by = "TRAIT")
  bar = bar[order(bar$idx, decreasing=T),]
  #rownames(bar) = bar$TRAIT

  cairo_pdf(file.path(output_dir, "barplot.pdf"), width = 8, height = 4, family = "Helvetica")
    barplot(t(bar[,c("inter_categorical", "intra_categorical", "single")]), ylim = c(0, 10), space = 0, col = c("black", "grey50", "white"))
  . = dev.off()
}


################################################################################
# additional plots

bar_df <- data.table(bar)
bar_df_m <-  melt(bar_df, measure.vars=c('total', 'inter_categorical', 'intra_categorical', 'single'))

bar_df_m$fill <- '#ffffff'
bar_df_m[variable == 'inter_categorical']$fill <- bar_df_m[variable == 'inter_categorical']$COLOR

bar_df_m$fill <- relevel(as.factor(bar_df_m$fill), '#ffffff')
bar_df_m$TRAIT <- factor(bar_df_m$TRAIT, levels=levels(fct_inorder(bar$TRAIT)))
bar_df_m$CATEGORY2 <- fct_rev(fct_inorder(bar_df_m$CATEGORY2))
bar_df_m$CATEGORY <- fct_rev(fct_inorder(bar_df_m$CATEGORY))

ggplot(bar_df_m[variable == 'total'], aes(TRAIT, value)) +
  geom_bar(aes(fill=I(fill), color=I(COLOR)), stat='identity', position='stack') +
  theme_classic() +
  theme(panel.spacing.y = unit(0, "lines"),
        strip.background = element_blank(),
        #strip.text.y.left = element_text(angle = 0),
        strip.text.y.left = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.line.y=element_blank()
  ) +
  scale_y_continuous(limits=c(0, 10), breaks = c(0, 5, 10)) +
  facet_grid(CATEGORY~., scales = "free_y", space = "free_y", switch = "y") +
  coord_flip()


ggsave(file.path(output_dir, "barplot.pdf"), width = 2, height = 3)

# get the unique number of SNPs / MARKERS per TRAIT
# should not use the number of rows
nsnps_per_trait = df %>%
  group_by(TRAIT) %>%
  summarize(
    marker_ct = n_distinct(MARKER),
    locus_ct = n_distinct(LOCUS_ID),
    gene_ct = n_distinct(GENE)
  )

nsnps_per_trait <- merge(nsnps_per_trait, unique(bar_df[, c('TRAIT', 'CATEGORY', 'COLOR')]))
nsnps_per_trait$TRAIT <- factor(nsnps_per_trait$TRAIT, levels=levels(fct_inorder(bar$TRAIT)))
#nsnps_per_trait$CATEGORY2 <- fct_rev(fct_inorder(nsnps_per_trait$CATEGORY2))
#nsnps_per_trait$CATEGORY <- fct_rev(fct_inorder(nsnps_per_trait$CATEGORY))
nsnps_per_trait$CATEGORY <- fct_rev(factor(nsnps_per_trait$CATEGORY, levels=levels(fct_inorder(bar$CATEGORY))))

ggplot(nsnps_per_trait, aes(TRAIT, marker_ct)) +
  geom_bar(aes(fill=I(COLOR), color=I(COLOR)), stat='identity', position='stack') +
  #coord_cartesian() +
  theme_classic() +
  theme(panel.spacing.y = unit(0, "lines"),
        strip.background = element_blank(),
        #strip.text.y.left = element_text(angle = 0),
        strip.text.y.left = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.line.y=element_blank()
  ) +
  scale_y_continuous(
    #limits=c(0, 10),
    breaks = c(1, 5, 10),
    #trans='log10'
    ) +
  facet_grid(CATEGORY~., scales = "free_y", space = "free_y", switch = "y") +
  coord_flip(ylim = c(0, 10))
ggsave(file.path(output_dir, "barplot_marker_count.pdf"), width = 1, height = 2)

ggplot(nsnps_per_trait, aes(TRAIT, gene_ct)) +
  geom_bar(aes(fill=I(COLOR), color=I(COLOR)), stat='identity', position='stack') +
  #coord_cartesian() +
  theme_classic() +
  theme(panel.spacing.y = unit(0, "lines"),
        strip.background = element_blank(),
        #strip.text.y.left = element_text(angle = 0),
        strip.text.y.left = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.line.y=element_blank()
  ) +
  scale_y_continuous(
    #limits=c(0, 10),
    breaks = c(1, 5, 10),
    #trans='log10'
  ) +
  facet_grid(CATEGORY~., scales = "free_y", space = "free_y", switch = "y") +
  coord_flip(ylim = c(0, 10))
ggsave(file.path(output_dir, "barplot_gene_count.pdf"), width = 1, height = 2)

ggplot(nsnps_per_trait, aes(TRAIT, marker_ct)) +
  geom_bar(aes(fill=I(COLOR), color=I(COLOR)), stat='identity', position='stack') +
  #coord_cartesian() +
  theme_classic() +
  theme(panel.spacing.y = unit(0, "lines"),
        strip.background = element_blank(),
        #strip.text.y.left = element_text(angle = 0),
        strip.text.y.left = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.line.y=element_blank()
  ) +
  scale_y_continuous(
    #limits=c(0, 10),
    breaks = c(10, 5, 1),
    #trans='log10'
    #limits=rev
  ) +
  facet_grid(CATEGORY~., scales = "free_y", space = "free_y", switch = "y") +
  coord_flip(ylim = c(10, 1))
ggsave(file.path(output_dir, "barplot_marker_count_rev.pdf"), width = 1, height = 2)

cols_df <- data.table(cols)
cols_df$cols_idx <- as.numeric(rownames(cols))
cols_df$colsep <- colsep
cols_df$category_upper <- toupper(cols_df$category_lower)

labels_plot <- ggplot(cols_df, aes(x=1, y=-cols_idx, label=category_upper, color=I(COLOR))) +
  geom_text_repel(
    #box.padding = 0.5,
    #nudge_y = 1,
    #segment.curvature = -0.1,
    #segment.ncp = 3,
    #segment.angle = 20,
    max.overlaps = Inf,
    min.segment.length = 0,
    #force        = 0.1,
    nudge_x      = -0.001,
    direction    = "x",
    #hjust        = 1,
    #segment.size = 0.2
  ) +
  theme_void(base_size = 20)
ggsave(file.path(output_dir, "circos_labels.pdf"), width = 3, height = 2)

writeLines(c("", "",
             sprintf("* Final circos outputs: %s.{png,svg}.", file.path(output_dir, 'circos')),
             "",
     sprintf("Finished at %s.", Sys.time())
           ))
