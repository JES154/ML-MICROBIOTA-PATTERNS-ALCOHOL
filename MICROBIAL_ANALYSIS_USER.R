# ============================================
# PACKAGES
# ============================================
pkgs <- c("phyloseq","dplyr","tidyr","ggplot2","data.table","forcats","vegan",
          "reshape2","ggridges","GGally","tidytext","scales","patchwork","ggpubr")
new <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Experiment paths ----
exp1 <- list(otu="C:/BIOINFORMATICA_TESIS_RESPALDO/FULL_CONTROLS/qiime2_output/exported/feature-table.tsv",
             tax="C:/BIOINFORMATICA_TESIS_RESPALDO/FULL_CONTROLS/qiime2_output/exported/taxonomy.tsv",
             meta="C:/BIOINFORMATICA_TESIS_RESPALDO/FULL_CONTROLS/META_CONTROLS.tsv", name="HC")

exp2 <- list(otu="C:/BIOINFORMATICA_TESIS_RESPALDO/ENA_PRJNA867698/PRJNA867698_ALCOHOL/qiime2_output2/exported/feature-table.tsv",
             tax="C:/BIOINFORMATICA_TESIS_RESPALDO/ENA_PRJNA867698/PRJNA867698_ALCOHOL/qiime2_output2/exported/taxonomy.tsv",
             meta="C:/BIOINFORMATICA_TESIS_RESPALDO/ENA_PRJNA867698/MD_60ALCOHOL.tsv", name="HD")

exp3 <- list(otu="C:/BIOINFORMATICA_TESIS_RESPALDO/ENA_PRJNA517050/ONLY_ALCOHOL/qiime2_output/exported/feature-table.tsv",
             tax="C:/BIOINFORMATICA_TESIS_RESPALDO/ENA_PRJNA517050/ONLY_ALCOHOL/qiime2_output/exported/taxonomy.tsv",
             meta="C:/BIOINFORMATICA_TESIS_RESPALDO/ENA_PRJNA517050/MD_118ALCOHOL.tsv", name="VHD")

# ============================================
# FUNCTION: BUILD PHYLOSEQ OBJECT
# ============================================
get_phy <- function(paths, grp_label) {
  require(data.table); require(phyloseq)
  cat("\n", rep("=",50), "\nProcessing experiment: ", grp_label, "\n", rep("=",50), "\n", sep="")
  
  tryCatch({
    # 1. OTU table
    cat("1. Reading OTU table...\n")
    otu <- fread(paths$otu, skip=1, header=TRUE)
    colnames(otu)[1] <- "FeatureID"
    sn <- colnames(otu)[-1]
    otu_mat <- t(as.matrix(otu[,-1])); rownames(otu_mat) <- sn; colnames(otu_mat) <- otu$FeatureID
    cat("   Samples:", length(sn), "| Features:", ncol(otu_mat), "\n")
    OTU <- otu_table(otu_mat, taxa_are_rows=FALSE)
    
    # 2. Taxonomy
    cat("2. Reading taxonomy...\n")
    tax <- fread(paths$tax, header=TRUE)
    if("Feature ID" %in% colnames(tax)) setnames(tax, "Feature ID", "Feature.ID")
    tax_mat <- do.call(cbind, tstrsplit(tax$Taxon, "; ", fill=NA))
    colnames(tax_mat) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    rownames(tax_mat) <- tax$Feature.ID
    for(i in 1:ncol(tax_mat)) {
      tax_mat[,i] <- gsub("^[kpcofgs]__","", tax_mat[,i])
      tax_mat[is.na(tax_mat[,i]) | tax_mat[,i]=="", i] <- "Unclassified"
    }
    TAX <- tax_table(tax_mat)
    
    # 3. Metadata
    cat("3. Reading metadata...\n")
    meta <- fread(paths$meta, header=TRUE)
    sid_col <- NULL
    for(col in colnames(meta)) {
      if(grepl("sample|Sample|ID|id|run|Run", col, ignore.case=TRUE)) {
        if(length(intersect(as.character(meta[[col]]), sn)) > 0) { sid_col <- col; break }
      }
    }
    if(is.null(sid_col)) sid_col <- colnames(meta)[1]
    meta_df <- as.data.frame(meta)
    meta_df[[sid_col]] <- as.character(meta_df[[sid_col]])
    rownames(meta_df) <- meta_df[[sid_col]]
    meta_df$Group <- grp_label
    meta_df <- meta_df[rownames(meta_df) %in% sn, , drop=FALSE]
    if(nrow(meta_df)==0) { cat("ERROR: No common samples between OTU and metadata\n"); return(NULL) }
    
    # 4. Align data
    cat("4. Aligning data...\n")
    common <- intersect(rownames(otu_mat), rownames(meta_df))
    otu_filt <- otu_mat[common, , drop=FALSE]
    meta_df  <- meta_df[common, , drop=FALSE]
    feats    <- intersect(colnames(otu_filt), rownames(TAX))
    OTU <- otu_table(otu_filt[, feats, drop=FALSE], taxa_are_rows=FALSE)
    TAX <- tax_table(tax_mat[feats, , drop=FALSE])
    META <- sample_data(meta_df)
    
    # 5. Build phyloseq
    ps <- phyloseq(OTU, TAX, META)
    cat("Phyloseq created! Samples:", nsamples(ps), "| Features:", ntaxa(ps),
        "| Total reads:", sum(sample_sums(ps)), "\n")
    return(ps)
    
  }, error=function(e) { cat("ERROR:", e$message, "\n"); return(NULL) })
}

# ============================================
# BUILD PHYLOSEQ OBJECTS
# ============================================
phy1 <- get_phy(exp1, exp1$name)
phy2 <- get_phy(exp2, exp2$name)
phy3 <- get_phy(exp3, exp3$name)

# ---- Summary ----
summarize_exp <- function(ps, name) {
  if(is.null(ps)) { cat("FAILED:", name, "\n"); return(FALSE) }
  cat(name, "- Samples:", nsamples(ps), "| Features:", ntaxa(ps),
      "| Reads:", sum(sample_sums(ps)), "| Group:", unique(sample_data(ps)$Group), "\n")
  return(TRUE)
}
ok <- c(summarize_exp(phy1,"HC"), summarize_exp(phy2,"HD"), summarize_exp(phy3,"VHD"))

# ---- Merge ----
valid <- Filter(Negate(is.null), list(phy1, phy2, phy3))
phy_all <- do.call(merge_phyloseq, valid)
phy_rel <- transform_sample_counts(phy_all, function(x) x/sum(x))
cat("Combined object - Samples:", nsamples(phy_all), "\n")
print(table(sample_data(phy_all)$Group))
save(phy_all, phy_rel, file="phyloseq_objects.RData")

# ============================================
# PART 2: RIDGE PLOTS & PHYLUM ANALYSIS
# ============================================
paleta_grupos <- c("HC"="#E57373", "HD"="#2ca02c", "VHD"="#1f77b4")
orden_personalizado <- c("HC", "HD", "VHD")
top_n <- 5

# ---- Prepare phylum data ----
phy_phylum_abs <- tax_glom(phy_all, taxrank="Phylum")
phy_phylum_rel <- tax_glom(phy_rel, taxrank="Phylum")

phylum_abs_df <- psmelt(phy_phylum_abs) %>%
  select(Sample, Group, Phylum, Abundance) %>%
  mutate(Group=factor(Group, levels=orden_personalizado), Log_Abundance=log10(Abundance+1))

phylum_rel_df <- psmelt(phy_phylum_rel) %>%
  select(Sample, Group, Phylum, Abundance) %>%
  rename(Relative_Abundance=Abundance) %>%
  mutate(Group=factor(Group, levels=orden_personalizado), Log_Relative=log10(Relative_Abundance*100+0.001))

top_phyla <- phylum_abs_df %>% group_by(Phylum) %>%
  summarise(Total=sum(Abundance)) %>% arrange(desc(Total)) %>% slice_head(n=top_n) %>% pull(Phylum)

phylum_abs_top <- filter(phylum_abs_df, Phylum %in% top_phyla)
phylum_rel_top <- filter(phylum_rel_df, Phylum %in% top_phyla)

phylum_order_abs <- phylum_abs_top %>% group_by(Phylum) %>%
  summarise(Mean=mean(Log_Abundance)) %>% arrange(desc(Mean)) %>% pull(Phylum)
phylum_abs_top$Phylum <- factor(phylum_abs_top$Phylum, levels=phylum_order_abs)

# ---- Ridge plot: relative ----
ridge_plot_rel <- ggplot(phylum_rel_top, aes(x=Relative_Abundance*100, y=Phylum, fill=Group)) +
  geom_density_ridges(aes(height=after_stat(density)), scale=0.95, alpha=0.75,
                      rel_min_height=0.01, color="black", linewidth=0.25, panel_scaling=FALSE) +
  scale_fill_manual(values=paleta_grupos, name="Experimental group") +
  scale_x_continuous(name="Relative abundance (%)", expand=expansion(mult=c(0.01,0.02)),
                     breaks=scales::pretty_breaks(n=5), labels=scales::percent_format(scale=1)) +
  scale_y_discrete(name=NULL, expand=expansion(add=c(0.2,0.4))) +
  labs(title=NULL, subtitle=NULL) +
  theme_ridges(font_size=11, grid=TRUE) +
  theme(axis.title.x=element_text(size=11,face="bold"), axis.text=element_text(size=9),
        legend.position="right", legend.title=element_text(face="bold",size=10),
        legend.text=element_text(size=9), plot.margin=ggplot2::margin(8,10,8,10))
ridge_plot_rel

# ---- Phylum statistics ----
phylum_stats <- phylum_abs_top %>%
  group_by(Phylum, Group) %>%
  summarise(N=n(), Mean_Abs=mean(Abundance), SD_Abs=sd(Abundance),
            Median_Abs=median(Abundance), Mean_Log=mean(Log_Abundance), SD_Log=sd(Log_Abundance), .groups="drop")
phylum_stats_rel <- phylum_rel_top %>%
  group_by(Phylum, Group) %>%
  summarise(Mean_Rel=mean(Relative_Abundance*100), SD_Rel=sd(Relative_Abundance*100),
            Median_Rel=median(Relative_Abundance*100), .groups="drop")
phylum_stats_complete <- left_join(phylum_stats, phylum_stats_rel, by=c("Phylum","Group"))
print(as.data.frame(head(phylum_stats_complete, 15)))
write.csv(phylum_stats_complete, "phylum_stats_per_group.csv", row.names=FALSE)

# ============================================
# PART 3: PHYLUM COUNTS BY GROUP
# ============================================
phy_phylum_abs2 <- tax_glom(phy_all, taxrank="Phylum")
phylum_counts   <- psmelt(phy_phylum_abs2)

phylum_summary_abs <- phylum_counts %>%
  group_by(Group, Phylum) %>% summarise(Total_Count=sum(Abundance), .groups="drop") %>%
  mutate(Group=factor(Group, levels=orden_personalizado)) %>% arrange(Group, desc(Total_Count))

# Print per-group top 10
for(grp in levels(phylum_summary_abs$Group)) {
  gd <- filter(phylum_summary_abs, Group==grp) %>% arrange(desc(Total_Count))
  cat("\nGroup:", grp, "\n")
  for(i in 1:min(10,nrow(gd)))
    cat(sprintf("  %-25s: %10d\n", substr(gd$Phylum[i],1,25), gd$Total_Count[i]))
  cat(sprintf("  %-25s: %10d\n", "TOTAL", sum(gd$Total_Count)))
}

phylum_percent <- phylum_summary_abs %>%
  group_by(Group) %>%
  mutate(Percent=Total_Count/sum(Total_Count)*100, Cumulative_Percent=cumsum(Percent)) %>%
  arrange(Group, desc(Percent))

top_phyla_overall <- phylum_summary_abs %>% group_by(Phylum) %>%
  summarise(Total=sum(Total_Count)) %>% arrange(desc(Total)) %>% slice_head(n=15) %>% pull(Phylum)

phylum_top <- filter(phylum_summary_abs, Phylum %in% top_phyla_overall)
phylum_top_percent <- filter(phylum_percent, Phylum %in% top_phyla_overall)
phylum_order <- phylum_top %>% group_by(Phylum) %>%
  summarise(Total=sum(Total_Count)) %>% arrange(desc(Total)) %>% pull(Phylum)
phylum_top$Phylum <- factor(phylum_top$Phylum, levels=phylum_order)
phylum_top_percent$Phylum <- factor(phylum_top_percent$Phylum, levels=phylum_order)

group_totals <- phylum_summary_abs %>% group_by(Group) %>%
  summarise(Total_Group=sum(Total_Count)) %>% arrange(factor(Group,levels=orden_personalizado))
subtitle_abs <- paste("Total reads: HC =", format(group_totals$Total_Group[1],big.mark=","),
                      "| HD =", format(group_totals$Total_Group[2],big.mark=","),
                      "| VHD =", format(group_totals$Total_Group[3],big.mark=","))

# INFORMATIVE AND SPECIFIC INFORMATION ABOUT SHARED AND UNIQUE PHYLUMS ACROSS GROUPS
# ---- Phylum statistics summary ----
phyla_per_group <- phylum_summary_abs %>% filter(Total_Count>0) %>% group_by(Group) %>%
  summarise(Unique_Phyla=n(), Total_Reads=sum(Total_Count),
            Mean_Reads=mean(Total_Count), Median_Reads=median(Total_Count)) %>%
  mutate(Group=factor(Group, levels=orden_personalizado)) %>% arrange(Group)
print(phyla_per_group)

phyla_by_group <- lapply(setNames(levels(phylum_summary_abs$Group), levels(phylum_summary_abs$Group)),
                         function(g) filter(phylum_summary_abs, Group==g, Total_Count>0) %>% pull(Phylum))
cat("Phyla shared across all groups:\n")
cat(paste(sort(Reduce(intersect, phyla_by_group)), collapse=", "), "\n")
for(grp in names(phyla_by_group)) {
  unicos <- setdiff(phyla_by_group[[grp]], unique(unlist(phyla_by_group[names(phyla_by_group)!=grp])))
  if(length(unicos)>0) cat("Unique in", grp, ":", paste(sort(unicos), collapse=", "), "\n")
  else cat("No unique phyla in group:", grp, "\n")
}

# ---- Export phylum results ----
if(!dir.exists("results")) dir.create("results")
write.csv(phylum_summary_abs, "results/phylum_counts_by_group_absolute.csv", row.names=FALSE)
write.csv(phylum_percent,     "results/phylum_counts_by_group_percent.csv",  row.names=FALSE)
write.csv(phyla_per_group,    "results/phylum_statistics_summary.csv",       row.names=FALSE)

# ============================================
# PART 4: DIVERSITY ANALYSIS
# ============================================

# ---- 4.1 Alpha diversity ----
sample_depths <- sample_sums(phy_all)
cat("Sequencing depth - Min:", min(sample_depths), "| Max:", max(sample_depths), "| Mean:", round(mean(sample_depths)), "\n")
if(max(sample_depths)/min(sample_depths) > 10)
  cat("WARNING: Large depth variation (>10x). Consider rarefaction.\n")

alpha_diversity <- estimate_richness(phy_all, measures=c("Observed","Chao1","Shannon","Simpson")) %>%
  mutate(Group=factor(sample_data(phy_all)$Group, levels=orden_personalizado),
         Sample=rownames(.), Read_Depth=sample_sums(phy_all)[rownames(.)]) %>%
  arrange(Group)

alpha_summary <- alpha_diversity %>% group_by(Group) %>%
  summarise(n_samples=n(),
            across(c(Observed,Chao1,Shannon,Simpson),
                   list(mean=~mean(.,na.rm=TRUE), sd=~sd(.,na.rm=TRUE)), .names="{.col}_{.fn}"),
            Mean_Read_Depth=mean(Read_Depth), .groups="drop")
print(as.data.frame(alpha_summary))

# ---- Statistical tests ----
perform_statistical_tests <- function(data, metric, metric_name) {
  cat("\n", metric_name, ":\n", sep="")
  kw <- kruskal.test(data[[metric]] ~ data$Group)
  cat("  Kruskal-Wallis: χ² =", round(kw$statistic,3), ", p =", format.pval(kw$p.value,digits=3), "\n")
  if(kw$p.value < 0.05) {
    cat("  Significant differences (p < 0.05)\n")
    if(requireNamespace("FSA",quietly=TRUE)) {
      library(FSA); print(dunnTest(data[[metric]] ~ data$Group, method="bonferroni")$res)
    } else {
      groups <- levels(data$Group)
      for(i in 1:(length(groups)-1)) for(j in (i+1):length(groups)) {
        wt <- wilcox.test(data[[metric]][data$Group==groups[i]],
                          data[[metric]][data$Group==groups[j]], exact=FALSE)
        sig <- cut(wt$p.value, breaks=c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*","ns"))
        cat(sprintf("  %s vs %s: p = %.4f %s\n", groups[i], groups[j], wt$p.value, sig))
      }
    }
  } else cat("  No significant differences\n")
  return(kw$p.value)
}

p_values <- list(
  Observed = perform_statistical_tests(alpha_diversity, "Observed", "Observed Species"),
  Chao1    = perform_statistical_tests(alpha_diversity, "Chao1",    "Chao1 Index"),
  Shannon  = perform_statistical_tests(alpha_diversity, "Shannon",  "Shannon Index"),
  Simpson  = perform_statistical_tests(alpha_diversity, "Simpson",  "Simpson Index")
)

# ---- Alpha diversity plots ----
metric_names <- c(Observed="Observed Species (count)", Chao1="Chao1 Index (estimated richness)",
                  Shannon="Shannon Index (units)", Simpson="Simpson Index (1-D)")
alpha_long <- alpha_diversity %>%
  select(Sample, Group, Observed, Chao1, Shannon, Simpson) %>%
  pivot_longer(cols=c(Observed,Chao1,Shannon,Simpson), names_to="Metric", values_to="Value") %>%
  mutate(Metric_Name=factor(metric_names[Metric], levels=metric_names))

p5 <- ggplot(alpha_long, aes(x=Group, y=Value, fill=Group)) +
  geom_boxplot(outlier.shape=NA, alpha=0.70, width=0.65) +
  geom_jitter(width=0.12, alpha=0.6, size=1.2, color="gray30") +
  scale_x_discrete(limits=orden_personalizado) +
  scale_fill_discrete(limits=orden_personalizado, guide="none") +
  facet_wrap(~Metric_Name, scales="free_y", ncol=2) +
  labs(title=NULL, subtitle=NULL, x="Experimental group", y="Diversity index value") +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle=30,hjust=1,size=9), axis.text.y=element_text(size=9),
        axis.title=element_text(size=10,face="bold"), strip.text=element_text(size=9,face="bold"),
        strip.background=element_rect(fill="gray92"), plot.title=element_text(hjust=0.5,face="bold",size=13),
        plot.subtitle=element_text(hjust=0.5,size=9), plot.margin=ggplot2::margin(10,12,10,12))
print(p5)


# ---- 4.2 Beta diversity ----
bray_dist   <- distance(phy_all, method="bray")
pcoa_result <- ordinate(phy_all, method="PCoA", distance="bray")
pcoa_scores <- as.data.frame(pcoa_result$vectors[,1:3])
colnames(pcoa_scores) <- paste0("PCoA",1:3)
pcoa_scores$Sample <- rownames(pcoa_scores)
pcoa_scores$Group  <- factor(sample_data(phy_all)$Group, levels=orden_personalizado)
eigvals <- pcoa_result$values$Eigenvalues
variance_explained <- round(eigvals[1:3]/sum(eigvals[eigvals>0])*100, 2)
cat("PCoA variance explained:", variance_explained[1], "%", variance_explained[2], "%", variance_explained[3], "%\n")

metadata_df     <- as.data.frame(as.matrix(sample_data(phy_all)))
permanova_result <- adonis2(bray_dist ~ Group, data=metadata_df)
cat("PERMANOVA: R² =", round(permanova_result$R2[1],4), "| F =", round(permanova_result$F[1],3),
    "| p =", format.pval(permanova_result$`Pr(>F)`[1],digits=4), "\n")

disp_result    <- betadisper(bray_dist, group=metadata_df$Group)
permdisp_result <- permutest(disp_result, permutations=999)
permdisp_f_value <- NA; permdisp_p_value <- NA
if("tab" %in% names(permdisp_result)) {
  permdisp_f_value <- permdisp_result$tab$F[1]
  permdisp_p_value <- permdisp_result$tab$`Pr(>F)`[1]
} else if("statistic" %in% names(permdisp_result)) {
  permdisp_f_value <- permdisp_result$statistic
  permdisp_p_value <- tryCatch(permdisp_result$permutations$p[1], error=function(e) NA)
}
cat("PERMDISP: F =", round(permdisp_f_value,3), "| p =", format.pval(permdisp_p_value,digits=4), "\n")

# ---- Beta diversity plots ----
p9 <- ggplot(pcoa_scores, aes(x=PCoA1, y=PCoA2)) +
  stat_ellipse(aes(fill=Group), geom="polygon", type="t", level=0.95, alpha=0.35, color=NA, show.legend=FALSE) +
  stat_ellipse(aes(color=Group), type="t", level=0.95, linewidth=0.7, show.legend=FALSE) +
  geom_point(aes(fill=Group), shape=21, size=3, color="black", stroke=0.35) +
  scale_fill_discrete(limits=orden_personalizado) + scale_color_discrete(limits=orden_personalizado) +
  labs(x=paste0("PCoA1 (",variance_explained[1],"%)"), y=paste0("PCoA2 (",variance_explained[2],"%)")) +
  theme_bw(base_size=11) +
  theme(axis.title=element_text(size=11,face="bold"), axis.text=element_text(size=9),
        panel.grid.minor=element_blank(), legend.position="none", plot.margin=ggplot2::margin(6,6,6,6))
if(nrow(pcoa_scores)<=30)
  p9 <- p9 + geom_text(aes(label=Sample), size=2.5, vjust=-0.7, check_overlap=TRUE)
print(p9)

p10 <- ggplot(pcoa_scores, aes(x=PCoA1, y=PCoA3, fill=Group)) +
  geom_point(size=3, alpha=0.85, shape=21, color="black", stroke=0.3) +
  stat_ellipse(aes(color=Group), type="t", level=0.95, linewidth=0.8, alpha=0.15) +
  scale_fill_discrete(limits=orden_personalizado) + scale_color_discrete(limits=orden_personalizado) +
  labs(title="Beta diversity: PCoA (PCoA1 vs PCoA3)",
       subtitle=paste0("PCoA1: ",variance_explained[1],"% | PCoA3: ",variance_explained[3],"%"),
       x=paste0("PCoA1 (",variance_explained[1],"%)"), y=paste0("PCoA3 (",variance_explained[3],"%)"),
       fill="Experimental group") +
  theme_bw(base_size=10) +
  theme(axis.title=element_text(size=10,face="bold"), axis.text=element_text(size=9),
        legend.position="right", legend.title=element_text(size=9,face="bold"),
        legend.text=element_text(size=8), plot.title=element_text(hjust=0.5,face="bold",size=13),
        plot.subtitle=element_text(hjust=0.5,size=9),
        panel.grid.major=element_line(color="gray92"), panel.grid.minor=element_blank(),
        plot.margin=ggplot2::margin(10,12,10,12))
print(p10)

# ---- Export diversity results ----
write.csv(alpha_diversity,              "results/alpha_diversity_metrics.csv",           row.names=FALSE)
write.csv(alpha_summary,               "results/alpha_diversity_summary.csv",            row.names=FALSE)
write.csv(pcoa_scores,                 "results/beta_diversity_pcoa_scores.csv",         row.names=FALSE)
write.csv(as.data.frame(as.matrix(bray_dist)), "results/beta_diversity_distance_matrix.csv")
write.csv(data.frame(Metric=names(p_values), Kruskal_Wallis_p=unlist(p_values),
                     Significant=unlist(p_values)<0.05),
          "results/statistical_test_results.csv", row.names=FALSE)
write.csv(as.data.frame(permanova_result), "results/permanova_results.csv")

# ---- Executive summary ----
cat("\nSTATISTICAL SIGNIFICANCE:\n")
for(m in names(p_values))
  cat(sprintf("  %s %-12s: p = %.4f %s\n", ifelse(p_values[[m]]<0.05,"✓","✗"), m,
              p_values[[m]], ifelse(p_values[[m]]<0.05,"(SIGNIFICANT)","")))
cat("  PERMANOVA: p =", format.pval(permanova_result$`Pr(>F)`[1],digits=4),
    ifelse(permanova_result$`Pr(>F)`[1]<0.05,"(SIGNIFICANT)",""), "\n")
if(!is.na(permdisp_p_value))
  cat("  PERMDISP:  p =", format.pval(permdisp_p_value,digits=4),
      ifelse(permdisp_p_value<0.05,"(SIGNIFICANT - heterogeneous variances)","(homogeneous variances)"), "\n")

cat("\nALPHA DIVERSITY PATTERNS:\n")
for(metric in c("Observed","Chao1","Shannon","Simpson")) {
  mm <- alpha_summary[[paste0(metric,"_mean")]]
  gs <- alpha_summary$Group
  trend <- if(all(diff(mm)>0)) paste("↑ Increases from",gs[1],"to",gs[length(gs)])
  else if(all(diff(mm)<0)) paste("↓ Decreases from",gs[1],"to",gs[length(gs)])
  else paste("↗ Maximum in", gs[which.max(mm)])
  cat(sprintf("  %-10s: %s (HC=%.1f, HD=%.1f, VHD=%.1f)\n", metric, trend, mm[1], mm[2], mm[3]))
}
cat("\nBETA DIVERSITY SUMMARY:\n")
cat("  PCoA1:", variance_explained[1], "% | PCoA2:", variance_explained[2], "% | Total:", sum(variance_explained[1:2]), "%\n")
cat("  PERMANOVA R²:", round(permanova_result$R2[1],3), "(",round(permanova_result$R2[1]*100,1),"% explained)\n")

# ============================================
# INTEGRATED FIGURE (PUBLICATION)
# ============================================
fig_diversity <- (ridge_plot_rel | p9) / p5 +
  plot_layout(heights=c(1,1.25), guides="collect") +
  plot_annotation(tag_levels="A") &
  theme(legend.position="right",
        legend.title=element_text(size=10,face="bold"), legend.text=element_text(size=9),
        plot.title=element_text(face="bold",size=15,hjust=0.5),
        plot.subtitle=element_text(size=11,hjust=0.5,color="gray35"),
        plot.tag=element_text(face="bold",size=14), plot.margin=ggplot2::margin(10,10,10,10))
fig_diversity

# SAVE FIGURES

ggsave("ridge_plot_abundancia_relativa.png", ridge_plot_rel, width=12, height=10, dpi=300)
ggsave("results/phylum_grouped_relative.png",  p5, width=11, height=9,  dpi=300)
ggsave("results/alpha_diversity_boxplots.png",          p9,  width=12, height=10, dpi=300)
ggsave("results/beta_diversity_pcoa_1_2.png", p13, width=10, height=8, dpi=300)
ggsave("results/beta_diversity_pcoa_1_3.png", p14, width=10, height=8, dpi=300)
ggsave("Figure_1_Integrated_Microbiota_Diversity.png",  fig_diversity, width=18, height=14, dpi=600,  units="in", bg="white")