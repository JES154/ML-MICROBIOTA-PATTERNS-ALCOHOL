# =========================================
# ---- PACKAGES ----
# =========================================
pkgs <- c("phyloseq","dplyr","tidyr","ggplot2","data.table","forcats","vegan","reshape2",
          "pheatmap","scales","stringr","ggnewscale","ggridges","tibble","patchwork",
          "rstatix","ggpubr")
new <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new)) install.packages(new)
lapply(pkgs, library, character.only = TRUE)

# =========================================
# ---- FUNCTIONS ----
# =========================================

# Load PICRUSt2 outputs
load_picrust_data <- function(picrust_path, experiment_name) {
  metagenome <- fread(file.path(picrust_path, "pathways_out", "path_abun_unstrat.tsv.gz"), header=TRUE, sep="\t") %>% as.data.frame()
  enzyme     <- fread(file.path(picrust_path, "EC_metagenome_out", "pred_metagenome_unstrat.tsv.gz"), header=TRUE, sep="\t") %>% as.data.frame()
  rownames(metagenome) <- metagenome[,1]; metagenome <- metagenome[,-1]
  rownames(enzyme)     <- enzyme[,1];     enzyme     <- enzyme[,-1]
  list(name=experiment_name, metagenome=metagenome, enzyme=enzyme)
}

# Get top N features by mean abundance
get_top_features <- function(feature_table, experiment_name, top_n=15) {
  avg_abun    <- rowMeans(feature_table)
  top_features <- sort(avg_abun, decreasing=TRUE)[1:top_n]
  data.frame(Feature=names(top_features), Abundance=top_features, Experiment=experiment_name)
}

# =========================================
# ---- EXPERIMENT PATHS ----
# =========================================
exp1 <- list(
  otu="HC/feature-table_controls.tsv",
  tax="HC/taxonomy_controls.tsv",
  meta="HC/META_CONTROLS.tsv",
  name="HC",
  picrust="HC/picrust2_out_HC"
)
exp2 <- list(
  otu="HD/feature-table_HD.tsv",
  tax="HD/taxonomy_HD.tsv",
  meta="HD/MD_60ALCOHOL.tsv",
  name="HD",
  picrust="HD/picrust2_out_HD"
)
exp3 <- list(
  otu="VHD/feature-table_VHD.tsv",
  tax="VHD/taxonomy_VHD.tsv",
  meta="VHD/MD_118ALCOHOL.tsv",
  name="VHD",
  picrust="/VHD/picrust2_out_VHD"
)

# =========================================
# ---- LOAD DATA ----
# =========================================
picrust_exp1 <- load_picrust_data(exp1$picrust, exp1$name)
picrust_exp2 <- load_picrust_data(exp2$picrust, exp2$name)
picrust_exp3 <- load_picrust_data(exp3$picrust, exp3$name)

# =========================================
# ---- MACHINE LEARNING: PATHWAYS vs ENZYMES ----
# =========================================
pkgs_ml <- c("caret","randomForest","e1071","xgboost","openxlsx","tibble","ggplot2","dplyr","MLmetrics")
new_ml  <- pkgs_ml[!(pkgs_ml %in% installed.packages()[,"Package"])]
if(length(new_ml)) install.packages(new_ml)
lapply(pkgs_ml, library, character.only = TRUE)

# =========================================
# ---- CROSS-VALIDATION SETTINGS ----
# =========================================
# multiClassSummary computes: Accuracy, Kappa, F1, Sensitivity, Specificity (per class)
ctrl <- trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats         = 3,
  classProbs      = TRUE,
  savePredictions = "final",
  summaryFunction = multiClassSummary
)

# =========================================
# ---- MODEL TRAINING ----
# =========================================
train_models <- function(features, labels, ctrl) {
  set.seed(123)
  
  # ----- RANDOM FOREST -----
  rf_model   <- train(x=features, y=labels, method="rf", trControl=ctrl, importance=TRUE)
  rf_metrics <- rf_model$results[which.max(rf_model$results$Accuracy),
                                 c("Accuracy","Kappa","Mean_F1","Mean_Sensitivity","Mean_Specificity")]
  
  # ----- SVM RADIAL -----
  set.seed(123)
  svm_model   <- train(x=features, y=labels, method="svmRadial", trControl=ctrl,
                       preProcess=c("center","scale"), tuneLength=5)
  svm_metrics <- svm_model$results[which.max(svm_model$results$Accuracy),
                                   c("Accuracy","Kappa","Mean_F1","Mean_Sensitivity","Mean_Specificity")]
  
  # ----- XGBOOST -----
  set.seed(123)
  xgb_model   <- train(x=features, y=labels, method="xgbTree", trControl=ctrl, tuneLength=3)
  xgb_metrics <- xgb_model$results[which.max(xgb_model$results$Accuracy),
                                   c("Accuracy","Kappa","Mean_F1","Mean_Sensitivity","Mean_Specificity")]
  
  # ----- COMPARISON -----
  perf_df <- rbind(
    cbind(Algorithm="Random Forest", rf_metrics),
    cbind(Algorithm="SVM Radial",    svm_metrics),
    cbind(Algorithm="XGBoost",       xgb_metrics)
  )
  perf_df        <- as.data.frame(perf_df)
  perf_df[,2:6]  <- sapply(perf_df[,2:6], as.numeric)
  perf_df$Global_Score <- rowMeans(perf_df[, c("Accuracy","Kappa","Mean_F1","Mean_Sensitivity","Mean_Specificity")])
  
  best_algo <- perf_df$Algorithm[which.max(perf_df$Global_Score)]
  message("[Best algorithm (Global Score): ", best_algo, "]")
  print(perf_df)
  
  # ---- Feature importances ----
  list(
    models = list(rf=rf_model, svm=svm_model, xgb=xgb_model),
    imps   = list(rf=varImp(rf_model)$importance, svm=varImp(svm_model)$importance, xgb=varImp(xgb_model)$importance),
    metrics = perf_df,
    best    = best_algo
  )
}

# =========================================
# ---- PRINT ML SUMMARY ----
# =========================================
print_model_summary <- function(models_list, best_algo_name, dataset_type) {
  cat("\n======================================\n")
  cat("[", dataset_type, "] Machine Learning Summary\n", sep="")
  cat("Best algorithm (Global Score):", best_algo_name, "\n\n")
  print(models_list$metrics)
  
  algo_key    <- switch(best_algo_name, "Random Forest"="rf", "SVM Radial"="svm", "XGBoost"="xgb")
  top_features <- head(models_list$imps[[algo_key]][order(-rowMeans(models_list$imps[[algo_key]])), , drop=FALSE], 10)
  cat("\nTop 10 features (", dataset_type, ", ", best_algo_name, "):\n", sep="")
  print(top_features)
  cat("======================================\n")
}

# =========================================
# ---- COMBINE PICRUSt2 DATA FOR ML ----
# =========================================
combine_picrust_data <- function(picrust_list, type=c("pathways","enzyme")) {
  type <- match.arg(type)
  combined_tables <- list()
  group_labels    <- c()
  
  for(pic in picrust_list) {
    tbl          <- if(type=="pathways") pic$metagenome else pic$enzyme
    group_labels <- c(group_labels, rep(pic$name, ncol(tbl)))
    combined_tables[[pic$name]] <- tbl
  }
  
  all_features <- unique(unlist(lapply(combined_tables, rownames)))
  
  # Align tables and fill missing features with 0
  tables_aligned <- lapply(combined_tables, function(tbl) {
    missing <- setdiff(all_features, rownames(tbl))
    if(length(missing) > 0) tbl[missing, ] <- 0
    tbl[all_features, , drop=FALSE]
  })
  
  combined    <- do.call(cbind, tables_aligned)
  combined_df <- as.data.frame(t(combined))   # samples as rows
  combined_df$Group <- factor(group_labels)
  combined_df
}

# =========================================
# ---- CREATE PICRUSt2 FEATURE NAME MAPPING ----
# =========================================
create_picrust_mapping <- function(picrust_path) {
  path_raw <- fread(file.path(picrust_path, "pathways_out",    "path_abun_unstrat.tsv.gz"))
  ec_raw   <- fread(file.path(picrust_path, "EC_metagenome_out","pred_metagenome_unstrat.tsv.gz"))
  list(
    pathways = data.frame(FeatureID=path_raw[[1]], FeatureName=path_raw[[1]], stringsAsFactors=FALSE),
    enzymes  = data.frame(FeatureID=ec_raw[[1]],  FeatureName=ec_raw[[1]],  stringsAsFactors=FALSE)
  )
}

# =========================================
# ---- RESOLVE EC NAMES FROM KEGG ----
# =========================================
library(KEGGREST)

get_ec_names <- function(ec_ids) {
  sapply(ec_ids, function(ec) {
    ec_clean <- gsub("EC:", "", ec)
    res <- tryCatch(keggGet(paste0("ec:", ec_clean))[[1]]$NAME, error=function(e) NA)
    ifelse(is.null(res), ec, paste0(ec, " - ", res[1]))
  })
}

# =========================================
# ---- PATHWAYS ANALYSIS ----
# =========================================
path_data      <- combine_picrust_data(list(picrust_exp1, picrust_exp2, picrust_exp3), "pathways")
path_data$Group <- factor(path_data$Group)
labels_path    <- factor(make.names(path_data$Group))
features_path  <- path_data[, !(colnames(path_data) %in% "Group")]

models_path    <- train_models(features_path, labels_path, ctrl)
best_algo_path <- models_path$best
message("[PATHWAYS] Best algorithm: ", best_algo_path)

algo_key      <- switch(best_algo_path, "Random Forest"="rf", "SVM Radial"="svm", "XGBoost"="xgb")
best_imp_path <- models_path$imps[[algo_key]]
path_top20    <- head(best_imp_path[order(-rowMeans(best_imp_path)), , drop=FALSE], 30)
print_model_summary(models_path, best_algo_path, "Pathways")

# =========================================
# ---- ENZYMES ANALYSIS ----
# =========================================
enz_data      <- combine_picrust_data(list(picrust_exp1, picrust_exp2, picrust_exp3), "enzyme")
enz_data$Group <- factor(enz_data$Group)
labels_enz    <- factor(make.names(enz_data$Group))
features_enz  <- enz_data[, !(colnames(enz_data) %in% "Group")]

models_enz    <- train_models(features_enz, labels_enz, ctrl)
best_algo_enz <- models_enz$best
message("[ENZYMES] Best algorithm: ", best_algo_enz)

algo_key     <- switch(best_algo_enz, "Random Forest"="rf", "SVM Radial"="svm", "XGBoost"="xgb")
best_imp_enz <- models_enz$imps[[algo_key]]
enz_top20    <- head(best_imp_enz[order(-rowMeans(best_imp_enz)), , drop=FALSE], 30)
print_model_summary(models_enz, best_algo_enz, "Enzymes")

# =========================================
# ---- PLOT: MODEL PERFORMANCE ----
# =========================================
metrics_all <- bind_rows(
  models_enz$metrics  %>% mutate(Level="Enzymes"),
  models_path$metrics %>% mutate(Level="Pathways")
)

metrics_long <- metrics_all %>%
  pivot_longer(cols=c(Accuracy,Kappa,Mean_F1,Mean_Sensitivity,Mean_Specificity,Global_Score),
               names_to="Metric", values_to="Value") %>%
  mutate(
    Metric = case_when(
      Metric == "Mean_F1"          ~ "F1 Score",
      Metric == "Mean_Sensitivity" ~ "Sensitivity",
      Metric == "Mean_Specificity" ~ "Specificity",
      Metric == "Global_Score"     ~ "Global Score",
      TRUE ~ Metric
    ),
    Metric = factor(Metric, levels=c("Accuracy","Kappa","F1 Score","Sensitivity","Specificity","Global Score"))
  )

p <- ggplot(metrics_long, aes(x=Algorithm, y=Value, fill=Algorithm)) +
  geom_col(position="dodge", width=0.65, color="black", linewidth=0.25) +
  facet_grid(Level ~ Metric, scales="free_y") +
  scale_fill_brewer(palette="Set2") +
  scale_y_continuous(expand=expansion(mult=c(0,0.03))) +
  theme_bw(base_size=9) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(size=7),
    axis.title       = element_text(size=9,face="bold"),
    plot.title       = element_text(size=11,face="bold",hjust=0.5),
    strip.text       = element_text(face="bold",size=8),
    strip.background = element_rect(fill="grey92",color="black",linewidth=0.4),
    legend.position  = "bottom",
    legend.title     = element_text(face="bold",size=8),
    legend.text      = element_text(size=7),
    legend.key.size  = unit(0.4,"cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color="grey88",linewidth=0.25),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color="black",linewidth=0.5)
  ) +
  labs(title=NULL, y="Score", x=NULL, fill="Algorithm")
print(p)

ggsave("model_performance.png",  p, width=12, height=6, dpi=600, bg="white")


# =========================================
# ---- PERMANOVA: FUNCTIONAL DIVERSITY ----
# =========================================
calc_permanova <- function(feature_df, group_col="Group", method="bray") {
  group   <- feature_df[[group_col]]
  num_df  <- feature_df[, sapply(feature_df, is.numeric)]
  valid   <- rowSums(num_df, na.rm=TRUE) > 0
  num_df  <- num_df[valid, , drop=FALSE]
  group   <- group[valid]
  num_df[is.na(num_df)] <- 0
  
  if(length(group) != nrow(num_df))
    stop("Inconsistency: number of groups != number of samples after cleaning.")
  if(nrow(num_df) < 3 || length(unique(group)) < 2) {
    warning("Too few samples or groups to run PERMANOVA.")
    return(NA)
  }
  
  vegan::adonis2(vegan::vegdist(num_df, method=method) ~ group, permutations=999)
}

permanova_enz <- calc_permanova(enz_data, group_col="Group", method="bray")
message("\n[PERMANOVA - ENZYMES]")
print(permanova_enz)

# =========================================
# ---- BETA DIVERSITY + PERMANOVA PLOT ----
# =========================================
plot_beta_div <- function(df, group_col="Group", method="bray", title="Beta Diversity") {
  group  <- df[[group_col]]
  num_df <- df[, sapply(df, is.numeric), drop=FALSE]
  valid  <- rowSums(num_df, na.rm=TRUE) > 0
  num_df <- num_df[valid, , drop=FALSE]
  group  <- droplevels(group[valid])
  num_df[is.na(num_df)] <- 0
  num_df <- sweep(num_df, 1, rowSums(num_df), "/")   # normalize to proportions
  
  dist_mat <- vegan::vegdist(num_df, method=method)
  pcoa_res <- cmdscale(dist_mat, eig=TRUE, k=2)
  pcoa_df  <- as.data.frame(pcoa_res$points)
  colnames(pcoa_df) <- c("PC1","PC2")
  pcoa_df$Group <- group
  
  var_exp <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)
  perm    <- vegan::adonis2(dist_mat ~ group, permutations=999)
  r2      <- round(perm$R2[1], 3)
  pval    <- signif(perm$`Pr(>F)`[1], 3)
  
  ggplot(pcoa_df, aes(x=PC1, y=PC2, color=Group)) +
    geom_point(size=4, alpha=0.9) +
    stat_ellipse(aes(fill=Group), geom="polygon", alpha=0.15, color=NA) +
    theme_minimal(base_size=14) +
    theme(panel.grid.major=element_line(color="gray80",linewidth=0.3),
          panel.grid.minor=element_line(color="gray90",linewidth=0.2),
          legend.position="right",
          plot.title=element_text(face="bold",hjust=0.5)) +
    labs(title=title, x=paste0("PC1 (",var_exp[1],"%)"), y=paste0("PC2 (",var_exp[2],"%)")) +
    annotate("text", x=min(pcoa_df$PC1), y=max(pcoa_df$PC2),
             label=paste0("R² = ",r2,"\nP = ",pval), hjust=0, vjust=1, size=5, color="black")
}

p_path <- plot_beta_div(path_data, "Group", "bray", "PCoA (Bray-Curtis) - Pathways")
print(p_path)
p_enz  <- plot_beta_div(enz_data,  "Group", "bray", "PCoA (Bray-Curtis) - Enzymes")
print(p_enz)

ggsave("p_path.png",  p_path, width=12, height=6, dpi=600, bg="white")
ggsave("p_enz.png",  p_enz, width=12, height=6, dpi=600, bg="white")

# =========================================
# ---- FEATURE NAME MAP ----
# =========================================
feature_name_map <- c(
  "METH-ACETATE-PWY"                        = "Methanogenesis from acetate",
  "3-HYDROXYPHENYLACETATE-DEGRADATION-PWY"  = "4-hydroxyphenylacetate degradation",
  "P124-PWY"                                = "Bifidobacterium shunt",
  "PWY-7374"                                = "5,8-dihydroxy-2-naphthoate biosynthesis I",
  "P122-PWY"                                = "Heterolactic fermentation",
  "EC:1.5.1.28"                             = "Opine dehydrogenase",
  "EC:1.12.1.2"                             = "Hydrogen dehydrogenase",
  "EC:2.8.3.17"                             = "3-(aryl)acryloyl-CoA:(R)-3-(aryl)lactate CoA-transferase",
  "EC:1.11.2.4"                             = "Fatty-acid peroxygenase",
  "EC:2.4.1.144"                            = "beta-1,4-mannosyl-glycoprotein 4-beta-N-acetylglucosaminyltransferase"
)

# =========================================
# ---- VIOLIN PLOTS WITH IMPORTANCE HIGHLIGHT ----
# =========================================
group_colors <- c("HC"="#E57373", "HD"="#029E73", "VHD"="#0173B2")

plot_violins_with_importance_highlight_log <- function(
    feature_df,
    top_features_df,
    feature_name_map  = NULL,
    group_col         = "Group",
    title_text        = NULL,
    ridge_colors      = c(group_colors),
    n_features        = 20,
    show_importance   = TRUE,
    show_legend       = TRUE,
    show_stats        = TRUE,
    importance_colors = c(low="lightgrey", high="black")
) {
  # ---- Top features ----
  importance_data <- top_features_df %>%
    rownames_to_column("Feature") %>%
    slice_head(n=n_features) %>%
    mutate(Feature_short=if(!is.null(feature_name_map)) recode(Feature, !!!feature_name_map) else Feature) %>%
    arrange(desc(Overall))
  ordered_features <- importance_data$Feature
  
  # ---- Long format ----
  df_long <- feature_df[, c(group_col, ordered_features), drop=FALSE] %>%
    reshape2::melt(id.vars=group_col, variable.name="Feature", value.name="Abundance") %>%
    mutate(
      Feature_short        = if(!is.null(feature_name_map)) recode(Feature, !!!feature_name_map) else Feature,
      Abundance_transformed = log10(Abundance + 1)
    )
  
  panel_top   <- max(df_long$Abundance_transformed, na.rm=TRUE)
  panel_range <- diff(range(df_long$Abundance_transformed, na.rm=TRUE))
  
  # ---- Dunn test ----
  if(show_stats) {
    dunn_results <- df_long %>%
      group_by(Feature) %>%
      dunn_test(as.formula(paste("Abundance_transformed ~", group_col)), p.adjust.method="bonferroni") %>%
      ungroup() %>%
      mutate(
        Feature_short = if(!is.null(feature_name_map)) recode(Feature, !!!feature_name_map) else Feature,
        p_raw  = p,
        p_adj  = p.adj,
        p.signif = case_when(
          p.adj <= 0.0001 ~ "****",
          p.adj <= 0.001  ~ "***",
          p.adj <= 0.01   ~ "**",
          p.adj <= 0.05   ~ "*",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(p.signif)) %>%
      add_xy_position(x="Feature_short", dodge=0.85) %>%
      group_by(Feature_short) %>%
      mutate(y.position = panel_top + panel_range * 0.22 + (row_number()-1) * panel_range * 0.07) %>%
      ungroup()
  }
  
  # ---- Medians ----
  df_stats <- df_long %>%
    group_by(Feature_short, !!sym(group_col)) %>%
    summarise(median_val=median(Abundance_transformed, na.rm=TRUE), .groups="drop")
  
  # ---- Base plot ----
  p <- ggplot(df_long,
              aes(x=reorder(Feature_short, match(Feature, ordered_features)),
                  y=Abundance_transformed, fill=!!sym(group_col))) +
    geom_violin(alpha=0.65, trim=FALSE, color="gray30", linewidth=0.3,
                scale="width", position=position_dodge(width=0.85)) +
    geom_point(data=df_stats,
               aes(x=Feature_short, y=median_val, fill=!!sym(group_col)),
               position=position_dodge(width=0.85),
               shape=21, size=2, color="white", stroke=0.4, show.legend=FALSE) +
    scale_fill_manual(values=ridge_colors, name="Group")
  
  # ---- Significance brackets ----
  if(show_stats && nrow(dunn_results) > 0)
    p <- p + stat_pvalue_manual(dunn_results, label="p.signif", tip.length=0.01,
                                bracket.size=0.3, size=3)
  
  # ---- Importance bar ----
  if(show_importance) {
    bar_y <- panel_top + panel_range * 0.55
    p <- p +
      ggnewscale::new_scale_fill() +
      geom_tile(data=importance_data,
                aes(x=Feature_short, y=bar_y, fill=Overall),
                inherit.aes=FALSE, width=0.7, height=panel_range*0.05,
                color="gray30", linewidth=0.2) +
      scale_fill_gradient(low=importance_colors["low"], high=importance_colors["high"],
                          name="Feature\nImportance")
  }
  
  # ---- Final formatting ----
  p + scale_x_discrete(expand=expansion(add=0.6),
                       labels=function(x) stringr::str_wrap(x, width=15)) +
    scale_y_continuous(expand=expansion(mult=c(0.02,0.05))) +
    labs(title=title_text,
         x="Features (ordered by importance)",
         y=expression(bold("Relative Abundance (" * log[10](x+1) * ")"))) +
    theme_minimal(base_size=10) +
    theme(
      plot.title   = element_text(face="bold", hjust=0.5),
      axis.title   = element_text(face="bold"),
      axis.text.x  = element_text(size=8, angle=30, hjust=1),
      panel.border = element_rect(color="gray30", fill=NA),
      legend.position = ifelse(show_legend,"right","none")
    )
}

# =========================================
# ---- GENERATE VIOLIN PLOTS ----
# =========================================
p_pathways_hybrid <- plot_violins_with_importance_highlight_log(
  feature_df        = path_data,
  top_features_df   = path_top20,
  feature_name_map  = feature_name_map,
  group_col         = "Group",
  title_text        = NULL,
  ridge_colors      = group_colors,
  n_features        = 5,
  show_importance   = TRUE,
  show_legend       = TRUE,
  importance_colors = c(low="#E8F4F8", high="#1E40AF")
)

p_enzymes_hybrid <- plot_violins_with_importance_highlight_log(
  feature_df        = enz_data,
  top_features_df   = enz_top20,
  feature_name_map  = feature_name_map,
  group_col         = "Group",
  title_text        = NULL,
  ridge_colors      = group_colors,
  n_features        = 5,
  show_importance   = TRUE,
  show_legend       = TRUE,
  importance_colors = c(low="#FEE2E2", high="#DC2626")
)

# ---- Combined multi-panel figure ----
combined_plot <- p_pathways_hybrid + p_enzymes_hybrid +
  plot_layout(ncol=1) +
  plot_annotation(tag_levels="A",
                  theme=theme(plot.tag=element_text(size=14, face="bold")))
print(combined_plot)

ggsave("Figure_Combined_Pathways_Enzymes.png",  combined_plot, width=12, height=10, dpi=600, bg="white")

# ggsave("Figure_Combined_Pathways_Enzymes.pdf",  combined_plot, device=cairo_pdf, width=174, height=200, units="mm")
# ggsave("Figure_Combined_Pathways_Enzymes.tiff", combined_plot, width=174, height=200, units="mm", dpi=600, compression="lzw")

# =========================================
# ---- LONG FORMAT FOR STATISTICAL TESTS ----
# =========================================
create_long_df <- function(feature_data, top_features_df, feature_name_map=NULL, n_top=5) {
  top_feats <- top_features_df %>%
    mutate(mean_importance=rowMeans(.)) %>%
    arrange(desc(mean_importance)) %>%
    head(n_top) %>%
    rownames()
  
  df_long <- feature_data[, top_feats, drop=FALSE] %>%
    mutate(Group=feature_data$Group) %>%
    reshape2::melt(id.vars="Group", variable.name="Feature", value.name="Abundance") %>%
    mutate(
      Feature_short         = if(!is.null(feature_name_map)) recode(Feature, !!!feature_name_map) else Feature,
      Abundance_transformed = log10(Abundance + 1)
    )
  df_long
}

df_path_long <- create_long_df(path_data, path_top20, n_top=5)
df_enz_long  <- create_long_df(enz_data,  enz_top20,  n_top=5)

# ---- Kruskal-Wallis ----
kw_path <- df_path_long %>% group_by(Feature_short) %>%
  do(kruskal_test(data=., Abundance_transformed ~ Group)) %>% ungroup()
kw_enz  <- df_enz_long  %>% group_by(Feature_short) %>%
  do(kruskal_test(data=., Abundance_transformed ~ Group)) %>% ungroup()

# ---- Dunn post-hoc ----
dunn_path <- df_path_long %>% group_by(Feature_short) %>%
  do(dunn_test(data=., Abundance_transformed ~ Group, p.adjust.method="bonferroni")) %>% ungroup()
dunn_enz  <- df_enz_long  %>% group_by(Feature_short) %>%
  do(dunn_test(data=., Abundance_transformed ~ Group, p.adjust.method="bonferroni")) %>% ungroup()

# ---- Effect size (eta²) ----
eta_path <- df_path_long %>% group_by(Feature_short) %>%
  do(kruskal_effsize(data=., Abundance_transformed ~ Group)) %>% ungroup()
eta_enz  <- df_enz_long  %>% group_by(Feature_short) %>%
  do(kruskal_effsize(data=., Abundance_transformed ~ Group)) %>% ungroup()

# ---- Summary tables ----
summary_path <- kw_path %>% select(Feature_short, p) %>% rename(KW_pvalue=p) %>%
  left_join(eta_path %>% select(Feature_short, effsize), by="Feature_short") %>% rename(Eta2=effsize)
summary_enz  <- kw_enz  %>% select(Feature_short, p) %>% rename(KW_pvalue=p) %>%
  left_join(eta_enz  %>% select(Feature_short, effsize), by="Feature_short") %>% rename(Eta2=effsize)

message("[Statistics generated successfully]")
print(summary_path)
print(summary_enz)

# ---- Export to Excel ----
if(!require("writexl")) install.packages("writexl")
library(writexl)
write_xlsx(summary_path, path="Summary_Pathways_Statistics.xlsx")
write_xlsx(summary_enz,  path="Summary_Enzymes_Statistics.xlsx")
write_xlsx(dunn_path,    path="Posthoc_Dunn_Pathways.xlsx")
write_xlsx(dunn_enz,     path="Posthoc_Dunn_Enzymes.xlsx")
