# ============================================================================
# PROYECTO BIOINFORMÁTICA: Senescencia en hMSC
# FASE 1: Descubrimiento de firmas transcriptómicas
# Autores: Eduardo Tello, Anjali Castro
# ============================================================================

# --- LIBRERÍAS ---
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(kableExtra)
  library(gridExtra)
  library(RColorBrewer)
  library(matrixStats)
  library(patchwork)
  library(ggrepel)
})

# --- PARÁMETROS GLOBALES ---
FDR_THRESHOLD <- 0.05
FC_THRESHOLD  <- 1.0

# ============================================================================
# 1. CARGA Y PREPROCESAMIENTO
# ============================================================================

# --- GSE39035: Senescencia cronológica (Edad × Pasaje) ---
gse39035    <- getGEO("GSE39035", GSEMatrix = TRUE)[[1]]
expr390_raw <- exprs(gse39035)
meta390_raw <- pData(gse39035)

# Metadata estructurada
meta_39035 <- meta390_raw %>%
  transmute(
    Sample    = geo_accession,
    Age_years = as.numeric(gsub("\\D", "", `age:ch1`)),
    Passage   = as.numeric(`passage:ch1`),
    Subject   = factor(`subject id:ch1`),
    CyDye     = factor(label_ch1)
  ) %>%
  dplyr::mutate(
    GroupEdad = factor(ifelse(Age_years >= 60, "Old", "Young"), 
                       levels = c("Young", "Old")),
    GroupPassage = factor(ifelse(Passage <= 4, "Low", "High"), 
                          levels = c("Low", "High"))
  )

# --- GSE7888: Senescencia replicativa (Early-Mid-Late) ---
gse7888     <- getGEO("GSE7888", GSEMatrix = TRUE)[[1]]
expr788_raw <- exprs(gse7888)
meta788_raw <- pData(gse7888)

meta_7888 <- meta788_raw %>%
  transmute(
    Sample  = geo_accession,
    Source  = source_name_ch1,
    Lot     = factor(str_extract(source_name_ch1, "lot#[A-Za-z0-9]+")),
    Passage = as.numeric(str_extract(source_name_ch1, "\\d+(?=(st|nd|rd|th) passage)"))
  ) %>%
  dplyr::mutate(
    GroupPassage = factor(dplyr::case_when(
      Passage <= 5  ~ "Early",
      Passage >= 22 ~ "Late",
      TRUE          ~ "Mid"
    ), levels = c("Early","Mid","Late"))
  )

# ============================================================================
# 2. CONTROL DE CALIDAD PRE-NORMALIZACIÓN
# ============================================================================

cat("\n=== CONTROL DE CALIDAD PRE-NORMALIZACIÓN ===\n")

# --- QC1: Boxplots de distribuciones crudas ---
par(mfrow = c(1, 2), mar = c(8, 4, 4, 2))
boxplot(expr390_raw, las = 2, main = "GSE39035 - Pre-normalización",
        ylab = "log2(intensidad)", cex.axis = 0.6, col = "lightblue")
boxplot(expr788_raw, las = 2, main = "GSE7888 - Pre-normalización",
        ylab = "log2(intensidad)", cex.axis = 0.6, col = "lightcoral")
par(mfrow = c(1, 1))

# --- QC2: Density plots ---
par(mfrow = c(1, 2))
plotDensities(expr390_raw, legend = FALSE, main = "GSE39035 - Densidades crudas")
plotDensities(expr788_raw, legend = FALSE, main = "GSE7888 - Densidades crudas")
par(mfrow = c(1, 1))

# --- QC3: Test de confounding batch-biología ---
cat("\n--- Test Chi-cuadrado: Batch vs Variables Biológicas ---\n")
confound_table_age <- table(meta_39035$Subject, meta_39035$GroupEdad)
confound_table_pass <- table(meta_39035$Subject, meta_39035$GroupPassage)

chi_age  <- chisq.test(confound_table_age)
chi_pass <- chisq.test(confound_table_pass)

cat(sprintf("Subject vs Edad: Chi² = %.2f, p = %.4f\n", chi_age$statistic, chi_age$p.value))
cat(sprintf("Subject vs Pasaje: Chi² = %.2f, p = %.4f\n", chi_pass$statistic, chi_pass$p.value))

if(chi_age$p.value < 0.05 | chi_pass$p.value < 0.05) {
  warning("ALERTA: Posible confounding entre batch y variables biológicas")
} else {
  cat("✓ No hay confounding significativo batch-biología\n")
}

# ============================================================================
# 3. ANOTACIÓN GÉNICA
# ============================================================================

# --- GPL13607 (Agilent) para GSE39035 ---
gpl_agilent   <- getGEO("GPL13607", AnnotGPL = TRUE)
annot_agilent <- Table(gpl_agilent)[, c("ID","GeneName")]
colnames(annot_agilent) <- c("ID","SYMBOL")

annot_agilent <- annot_agilent %>%
  dplyr::filter(SYMBOL != "", SYMBOL != "---") %>%
  dplyr::mutate(ID = as.character(ID))

# Colapso por mediana (manejo de probes múltiples)
expr390_annot <- expr390_raw %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  inner_join(annot_agilent, by = "ID") %>%
  dplyr::select(-ID) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(where(is.numeric), ~median(.x, na.rm = TRUE))) %>%
  as.data.frame() %>%
  column_to_rownames("SYMBOL")

# --- GPL570 (Affymetrix) para GSE7888 ---
gpl_affy   <- getGEO("GPL570", AnnotGPL = TRUE)
annot_affy <- Table(gpl_affy)[, c("ID","Gene symbol")]
colnames(annot_affy) <- c("ID","SYMBOL")

annot_affy <- annot_affy %>%
  dplyr::filter(SYMBOL != "", SYMBOL != "---") %>%
  dplyr::mutate(ID = as.character(ID))

expr788_annot <- expr788_raw %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  inner_join(annot_affy, by = "ID") %>%
  dplyr::select(-ID) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(where(is.numeric), ~median(.x, na.rm = TRUE))) %>%
  as.data.frame() %>%
  column_to_rownames("SYMBOL")

# ============================================================================
# 4. NORMALIZACIÓN Y FILTRADO
# ============================================================================

# --- Normalización cuantílica (between-array) ---
expr390_norm <- normalizeBetweenArrays(expr390_annot, method = "quantile")
expr390_norm <- expr390_norm[, meta_39035$Sample]

expr788_log  <- log2(expr788_annot + 1)  # Affymetrix requiere log2 manual
expr788_norm <- normalizeBetweenArrays(expr788_log, method = "quantile")
expr788_norm <- expr788_norm[, meta_7888$Sample]

# --- QC POST-NORMALIZACIÓN ---
cat("\n=== CONTROL DE CALIDAD POST-NORMALIZACIÓN ===\n")

par(mfrow = c(1, 2), mar = c(8, 4, 4, 2))
boxplot(expr390_norm, las = 2, main = "GSE39035 - Post-normalización",
        ylab = "log2(intensidad)", cex.axis = 0.6, col = "lightgreen")
boxplot(expr788_norm, las = 2, main = "GSE7888 - Post-normalización",
        ylab = "log2(intensidad)", cex.axis = 0.6, col = "lightgreen")
par(mfrow = c(1, 1))

# --- Filtrado por expresión (mediana > log2(100)) ---
keep390 <- rowMedians(expr390_norm) > log2(100)
keep788 <- rowMedians(expr788_norm) > log2(100)

expr390_filt <- expr390_norm[keep390, ]
expr788_filt <- expr788_norm[keep788, ]

cat(sprintf("\nGSE39035: %d/%d genes retenidos (%.1f%%)\n", 
            sum(keep390), length(keep390), 100*sum(keep390)/length(keep390)))
cat(sprintf("GSE7888: %d/%d genes retenidos (%.1f%%)\n", 
            sum(keep788), length(keep788), 100*sum(keep788)/length(keep788)))

# ============================================================================
# 5. ANÁLISIS DE SENSIBILIDAD DE UMBRALES
# ============================================================================

cat("\n=== ANÁLISIS DE SENSIBILIDAD DE UMBRALES ===\n")

# Función para contar DEGs bajo diferentes umbrales
sensitivity_analysis <- function(tt_full, name) {
  fc_thresholds  <- c(0.5, 1.0, 1.5, 2.0)
  fdr_thresholds <- c(0.01, 0.05, 0.1)
  
  cat(sprintf("\n--- %s ---\n", name))
  cat("FDR\\FC\t0.5\t1.0\t1.5\t2.0\n")
  
  for(fdr in fdr_thresholds) {
    cat(sprintf("%.2f\t", fdr))
    for(fc in fc_thresholds) {
      n_degs <- sum(tt_full$adj.P.Val < fdr & abs(tt_full$logFC) >= fc, na.rm = TRUE)
      cat(sprintf("%d\t", n_degs))
    }
    cat("\n")
  }
}

# ============================================================================
# 6. ANÁLISIS EXPLORATORIO (PCA)
# ============================================================================

pca_df <- function(mat, meta) {
  pca <- prcomp(t(mat), scale. = TRUE)
  df  <- as.data.frame(pca$x) %>%
    rownames_to_column("Sample") %>%
    left_join(meta, by = "Sample")
  list(pca = pca, df = df)
}

theme_pca <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face="bold", hjust=0),
    axis.title = element_text(face="bold"),
    panel.border = element_rect(color="grey70", fill=NA, linewidth=0.4),
    panel.grid.minor = element_blank()
  )

# --- PCA GSE39035 (con corrección batch) ---
expr390_corr <- removeBatchEffect(
  expr390_filt,
  batch  = meta_39035$Subject,
  design = model.matrix(~ GroupEdad + GroupPassage, meta_39035)
)

pca390 <- pca_df(expr390_corr, meta_39035)
var390 <- summary(pca390$pca)$importance[2, 1:2] * 100

p_390 <- ggplot(pca390$df, aes(PC1, PC2, color = GroupEdad, shape = GroupPassage)) +
  geom_point(size = 3) +
  labs(
    title = "A) GSE39035 — Edad × Pasaje",
    x = sprintf("PC1 (%.1f%%)", var390[1]),
    y = sprintf("PC2 (%.1f%%)", var390[2])
  ) +
  scale_color_manual(values = c("Young"="#3b82f6", "Old"="#ef4444")) +
  theme_pca

# --- PCA GSE7888 ---
expr788_corr <- removeBatchEffect(
  expr788_filt,
  batch  = meta_7888$Lot,
  design = model.matrix(~ GroupPassage, meta_7888)
)

pca788 <- pca_df(expr788_corr, meta_7888)
var788 <- summary(pca788$pca)$importance[2, 1:2] * 100

p_788 <- ggplot(pca788$df, aes(PC1, PC2, color = GroupPassage)) +
  geom_point(size = 3) +
  labs(
    title = "B) GSE7888 — Early–Mid–Late",
    x = sprintf("PC1 (%.1f%%)", var788[1]),
    y = sprintf("PC2 (%.1f%%)", var788[2])
  ) +
  scale_color_manual(values = c("Early"="#22c55e","Mid"="#f97316","Late"="#ef4444")) +
  theme_pca

# Figura combinada
print((p_390 | p_788) + plot_layout(guides="collect"))

# ============================================================================
# 7. ANÁLISIS DIFERENCIAL (limma con duplicateCorrelation)
# ============================================================================

fit_limma_block <- function(expr_mat, design, contrast_vec, block, contrast_name,
                            fdr_cut = FDR_THRESHOLD, fc_cut = FC_THRESHOLD) {
  
  # Correlación intra-bloque (para medidas repetidas/pareadas)
  cor_fit <- duplicateCorrelation(expr_mat, design, block = block)
  
  cat(sprintf("\n--- %s ---\n", contrast_name))
  cat(sprintf("Correlación intra-bloque (ρ) = %.3f\n", cor_fit$consensus.correlation))
  
  # Modelo lineal con estructura de correlación
  fit <- lmFit(expr_mat, design,
               block       = block,
               correlation = cor_fit$consensus.correlation)
  
  # Contraste y moderación Bayesiana
  fit2 <- contrasts.fit(fit, contrast_vec)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  
  # Tabla completa
  tt <- topTable(fit2, number = Inf, sort.by = "P") %>%
    rownames_to_column("Gene")
  
  # DEGs significativos
  sig <- tt %>%
    dplyr::filter(adj.P.Val < fdr_cut, abs(logFC) >= fc_cut)
  
  cat(sprintf("DEGs significativos: %d (Up: %d, Down: %d)\n", 
              nrow(sig), 
              sum(sig$logFC > 0), 
              sum(sig$logFC < 0)))
  
  list(
    fit2        = fit2,
    full_table  = tt,
    sig_table   = sig,
    rho         = cor_fit$consensus.correlation,
    contrast    = contrast_name
  )
}

# --- CONTRASTE 1: Edad (Old vs Young) ---
mt390 <- meta_39035
X390  <- expr390_filt

design_age <- model.matrix(~ GroupEdad + GroupPassage, data = mt390)
cont_age   <- makeContrasts(Old_vs_Young = GroupEdadOld, levels = design_age)

res_age      <- fit_limma_block(X390, design_age, cont_age, mt390$Subject, "Edad (Old vs Young)")
DEG_age_full <- res_age$full_table
DEG_age      <- res_age$sig_table

# --- CONTRASTE 2: Pasaje (High vs Low) ---
design_pass <- model.matrix(~ GroupPassage + GroupEdad, data = mt390)
cont_pass   <- makeContrasts(High_vs_Low = GroupPassageHigh, levels = design_pass)

res_pass      <- fit_limma_block(X390, design_pass, cont_pass, mt390$Subject, "Pasaje (High vs Low)")
DEG_pass_full <- res_pass$full_table
DEG_pass      <- res_pass$sig_table

# --- CONTRASTE 3: Aditivo (Edad + Pasaje) ---
design_add <- model.matrix(~ GroupEdad + GroupPassage, data = mt390)
cont_add   <- makeContrasts(Additive = GroupEdadOld + GroupPassageHigh, levels = design_add)

res_add      <- fit_limma_block(X390, design_add, cont_add, mt390$Subject, "Aditivo (Edad + Pasaje)")
DEG_add_full <- res_add$full_table
DEG_add      <- res_add$sig_table

# --- CONTRASTE 4: Interacción (Edad × Pasaje) ---
design_int <- model.matrix(~ GroupEdad * GroupPassage, data = mt390)
colnames(design_int) <- make.names(colnames(design_int))
cont_int <- makeContrasts(Interaction = GroupEdadOld.GroupPassageHigh, levels = design_int)

res_int      <- fit_limma_block(X390, design_int, cont_int, mt390$Subject, "Interacción Edad×Pasaje")
DEG_int_full <- res_int$full_table
DEG_int      <- res_int$sig_table

# --- GSE7888: Late vs Early ---
X788  <- expr788_filt
mt788 <- meta_7888

mt_le <- mt788 %>% dplyr::filter(GroupPassage %in% c("Early","Late")) %>% droplevels()
X_le  <- X788[, mt_le$Sample]
mt_le$GroupPassage <- factor(mt_le$GroupPassage, levels = c("Early","Late"))

design_le <- model.matrix(~ GroupPassage, data = mt_le)
cont_le   <- makeContrasts(Late_vs_Early = GroupPassageLate, levels = design_le)

res_le        <- fit_limma_block(X_le, design_le, cont_le, mt_le$Lot, "Late vs Early")
DEG_7888_full <- res_le$full_table
DEG_7888      <- res_le$sig_table

# --- Análisis de sensibilidad ---
sensitivity_analysis(DEG_age_full, "Edad (Old vs Young)")
sensitivity_analysis(DEG_pass_full, "Pasaje (High vs Low)")
sensitivity_analysis(DEG_7888_full, "Late vs Early")

# ============================================================================
# 8. VALIDACIÓN ESTADÍSTICA DEL MODELO
# ============================================================================
# 8. VALIDACIÓN ESTADÍSTICA DEL MODELO
# ============================================================================

cat("\n=== VALIDACIÓN DE SUPUESTOS DEL MODELO LINEAL ===\n")

# --- Diagnóstico 1: Mean-variance trend ---
par(mfrow = c(2, 2))
plotSA(res_age$fit2, main = "Edad: Mean-variance trend")
plotSA(res_pass$fit2, main = "Pasaje: Mean-variance trend")
plotSA(res_le$fit2, main = "Replicativa: Mean-variance trend")

# --- Diagnóstico 2: Distribución de p-valores ---
# CORRECCIÓN: Usar histogramas de p-valores en lugar de Q-Q plot de residuos
hist(DEG_age_full$P.Value, breaks = 50, main = "Edad: Distribución p-valores",
     xlab = "p-valor", col = "lightblue", border = "white")
abline(h = nrow(DEG_age_full)/50, col = "red", lty = 2, lwd = 2)
text(0.8, nrow(DEG_age_full)/50 + 20, "Uniforme esperado", col = "red", pos = 3)

par(mfrow = c(1, 1))

cat("✓ Si plotSA muestra línea horizontal → varianza constante\n")
cat("✓ Si histograma p-valores tiene pico en 0 + uniforme → modelo adecuado\n")

# ============================================================================
# 9. TABLA RESUMEN DE DEGs
# ============================================================================

count_degs <- function(deg_sig, contrast_name) {
  total <- nrow(deg_sig)
  up <- sum(deg_sig$logFC > 0)
  down <- sum(deg_sig$logFC < 0)
  
  tibble(
    Contraste     = contrast_name,
    Total_DEGs    = total,
    Upregulated   = up,
    Downregulated = down,
    Porc_Up       = ifelse(total > 0, round(100 * up / total, 1), 0),
    Porc_Down     = ifelse(total > 0, round(100 * down / total, 1), 0)
  )
}

tabla_deg_resumen <- dplyr::bind_rows(
  count_degs(DEG_age,   "1. Edad (Old vs Young)"),
  count_degs(DEG_pass,  "2. Pasaje (High vs Low)"),
  count_degs(DEG_add,   "3. Aditivo (Edad + Pasaje)"),
  count_degs(DEG_int,   "4. Interacción Edad×Pasaje"),
  count_degs(DEG_7888,  "5. Late vs Early")
)

print(tabla_deg_resumen %>%
        kable(caption = "Número de DEGs por contraste (FDR<0.05, |logFC|≥1)",
              col.names = c("Contraste", "Total DEGs", "↑ Up", "↓ Down", "% Up", "% Down")) %>%
        kable_styling(bootstrap_options = c("striped","hover")))

# ============================================================================
# 10. VOLCANO PLOTS
# ============================================================================

plot_volcano <- function(df, titulo, logFC_cut = FC_THRESHOLD, fdr_cut = FDR_THRESHOLD) {
  df2 <- df %>%
    dplyr::mutate(
      status = dplyr::case_when(
        adj.P.Val < fdr_cut & logFC >=  logFC_cut ~ "Up",
        adj.P.Val < fdr_cut & logFC <= -logFC_cut ~ "Down",
        TRUE ~ "NoSig"
      ),
      status = factor(status, levels = c("Up","Down","NoSig"))
    )
  
  ggplot(df2, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
    geom_point(alpha = 0.6, size = 1.6) +
    scale_color_manual(
      values = c(Up = "#ef4444", Down = "#3b82f6", NoSig = "#d1d5db"),
      labels = c("Upregulated", "Downregulated", "No significativo")
    ) +
    geom_vline(xintercept = c(-logFC_cut, logFC_cut), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "grey50") +
    labs(title = titulo, x = "log2(Fold Change)", y = expression(-log[10](FDR))) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold"),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.4),
      legend.position = "bottom"
    )
}

v_age  <- plot_volcano(DEG_age_full,  "A) Edad (Old vs Young)")
v_pass <- plot_volcano(DEG_pass_full, "B) Pasaje (High vs Low)")
v_late <- plot_volcano(DEG_7888_full, "C) Senescencia replicativa (Late vs Early)")

print((v_age | v_pass | v_late) + plot_layout(guides = "collect"))

# ============================================================================
# 11. HEATMAPS (Top 30 DEGs)
# ============================================================================

paleta <- colorRampPalette(rev(brewer.pal(7, "RdBu")))(101)

# --- Top 30 Edad ---
top_age   <- DEG_age_full %>% dplyr::slice_min(adj.P.Val, n=30) %>% dplyr::pull(Gene)
mat_age   <- expr390_filt[top_age, meta_39035$Sample]
mat_age_z <- t(scale(t(mat_age)))

# CORRECCIÓN: Eliminar rownames antes de column_to_rownames
ann_age <- meta_39035 %>%
  dplyr::select(Sample, GroupEdad, GroupPassage) %>%
  as.data.frame() %>%                              # Asegurar que es data.frame
  `rownames<-`(NULL) %>%                           # Eliminar rownames existentes
  column_to_rownames("Sample")

p_h1 <- pheatmap(
  mat_age_z,
  annotation_col = ann_age,
  color          = paleta,
  show_colnames  = FALSE,
  fontsize_row   = 7,
  main           = "A) Edad (Old vs Young)",
  silent         = TRUE
)

# --- Top 30 Replicativa ---
top_rep   <- DEG_7888_full %>% dplyr::slice_min(adj.P.Val, n=30) %>% dplyr::pull(Gene)
mat_rep   <- expr788_filt[top_rep, meta_7888$Sample]
mat_rep_z <- t(scale(t(mat_rep)))

# CORRECCIÓN: Eliminar rownames antes de column_to_rownames
ann_rep <- meta_7888 %>%
  dplyr::select(Sample, GroupPassage, Lot) %>%
  as.data.frame() %>%                              # Asegurar que es data.frame
  `rownames<-`(NULL) %>%                           # Eliminar rownames existentes
  column_to_rownames("Sample")

p_h2 <- pheatmap(
  mat_rep_z,
  annotation_col = ann_rep,
  color          = paleta,
  show_colnames  = FALSE,
  fontsize_row   = 7,
  main           = "B) Replicativa (Late vs Early)",
  silent         = TRUE
)

grid.arrange(p_h1$gtable, p_h2$gtable, ncol = 2)

# ============================================================================
# 12. FIRMA CORE DE SENESCENCIA
# ============================================================================

# Intersección de genes con dirección concordante
age_up   <- DEG_age %>% dplyr::filter(logFC > 0) %>% dplyr::pull(Gene)
age_down <- DEG_age %>% dplyr::filter(logFC < 0) %>% dplyr::pull(Gene)

rep_up   <- DEG_7888 %>% dplyr::filter(logFC > 0) %>% dplyr::pull(Gene)
rep_down <- DEG_7888 %>% dplyr::filter(logFC < 0) %>% dplyr::pull(Gene)

core_up   <- intersect(age_up, rep_up)
core_down <- intersect(age_down, rep_down)
core_all  <- c(core_up, core_down)

cat(sprintf("\n=== FIRMA CORE DE SENESCENCIA: %d genes ===\n", length(core_all)))
cat(sprintf("Up-regulated en ambos: %d\n", length(core_up)))
cat(sprintf("Down-regulated en ambos: %d\n", length(core_down)))

if (length(core_all) > 0) {
  core_signature <- tibble(
    Gene       = core_all,
    Direccion  = ifelse(core_all %in% core_up, "Up en ambos", "Down en ambos"),
    logFC_Edad = DEG_age_full$logFC[match(core_all, DEG_age_full$Gene)],
    logFC_Rep  = DEG_7888_full$logFC[match(core_all, DEG_7888_full$Gene)],
    FDR_Edad   = DEG_age_full$adj.P.Val[match(core_all, DEG_age_full$Gene)],
    FDR_Rep    = DEG_7888_full$adj.P.Val[match(core_all, DEG_7888_full$Gene)]
  ) %>%
    dplyr::arrange(desc(abs(logFC_Edad)))
  
  print(core_signature %>%
          dplyr::mutate(
            logFC_Edad = sprintf("%.2f", logFC_Edad),
            logFC_Rep  = sprintf("%.2f", logFC_Rep),
            FDR_Edad   = formatC(FDR_Edad, format = "e", digits = 2),
            FDR_Rep    = formatC(FDR_Rep, format = "e", digits = 2)
          ) %>%
          kable(caption = "Firma CORE de senescencia (dirección concordante)") %>%
          kable_styling(bootstrap_options = c("striped","hover")))
  
  # Guardar firma CORE
  write.csv(core_signature, "CORE_senescence_signature.csv", row.names = FALSE)
  cat("\n✓ Firma CORE guardada en 'CORE_senescence_signature.csv'\n")
}

# ============================================================================
# 13. EXPORTAR RESULTADOS PARA FASE 2
# ============================================================================

cat("\n=== EXPORTANDO DATOS PARA FASE 2 ===\n")

# Listas de DEGs por contraste
degs_export <- list(
  DEG_age_full  = DEG_age_full,
  DEG_age_sig   = DEG_age,
  DEG_pass_full = DEG_pass_full,
  DEG_pass_sig  = DEG_pass,
  DEG_rep_full  = DEG_7888_full,
  DEG_rep_sig   = DEG_7888,
  core_signature = if(exists("core_signature")) core_signature else NULL
)

saveRDS(degs_export, "FASE1_DEGs_export.rds")
cat("✓ DEGs exportados a 'FASE1_DEGs_export.rds'\n")

# Matrices de expresión normalizadas
expr_export <- list(
  expr390_norm = expr390_norm,
  expr390_filt = expr390_filt,
  expr788_norm = expr788_norm,
  expr788_filt = expr788_filt,
  meta_39035   = meta_39035,
  meta_7888    = meta_7888
)

saveRDS(expr_export, "FASE1_expr_matrices.rds")
cat("✓ Matrices de expresión exportadas a 'FASE1_expr_matrices.rds'\n")

# ============================================================================
# 14. RESUMEN EJECUTIVO
# ============================================================================

cat("\n" ,"=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("                     RESUMEN EJECUTIVO FASE 1\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

cat(sprintf("Dataset GSE39035: %d muestras, %d genes post-filtrado\n", 
            ncol(expr390_filt), nrow(expr390_filt)))
cat(sprintf("Dataset GSE7888:  %d muestras, %d genes post-filtrado\n", 
            ncol(expr788_filt), nrow(expr788_filt)))

cat("\nDEGs por contraste (FDR<0.05, |logFC|≥1):\n")
print(tabla_deg_resumen %>% dplyr::select(Contraste, Total_DEGs, Upregulated, Downregulated))

cat(sprintf("\nFirma CORE: %d genes con dirección concordante\n", length(core_all)))
cat(sprintf("Correlación intra-bloque (ρ) promedio: %.3f\n", 
            mean(c(res_age$rho, res_pass$rho, res_le$rho))))

cat("\n✓ Todos los gráficos generados exitosamente\n")
cat("✓ Datos exportados para Fase 2\n")
cat("\n" ,"=" %>% rep(70) %>% paste(collapse=""), "\n\n")

# ============================================================================
# FIN DE FASE 1
# ============================================================================
