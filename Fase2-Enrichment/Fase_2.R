# ============================================================================
# FASE 2: ANÁLISIS DE ENRIQUECIMIENTO E INTEGRACIÓN INMUNOMODULADORA
# Proyecto: Predicción y Rescate de Función Inmunomoduladora en hMSC
# Autores: Eduardo Tello & Anjali Castro
# Curso: Bioinformática Avanzada UTEC 2025-2
# ============================================================================


# ============================================================================
# CONFIGURACIÓN INICIAL
# ============================================================================

# Librerías principales
library(tidyverse)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(igraph)
library(ggraph)
library(pROC)
library(ComplexHeatmap)
library(circlize)
library(patchwork)

# Configuración de gráficos
theme_set(theme_minimal(base_size = 13))

# Crear directorios
dir.create("Fase2-Enrichment/results", recursive = TRUE, showWarnings = FALSE)
dir.create("Fase2-Enrichment/figures", recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. CARGA DE DATOS FASE 1
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("CARGANDO DATOS DE FASE 1\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Cargar DEGs
degs_data <- readRDS("Fase1-Discovery/results/FASE1_DEGs_export.rds")
DEG_age_full  <- degs_data$DEG_age_full
DEG_7888_full <- degs_data$DEG_rep_full
core_signature <- degs_data$core_signature

# Cargar expresión
expr_data <- readRDS("Fase1-Discovery/results/FASE1_expr_matrices.rds")
expr390_filt <- expr_data$expr390_filt
expr788_filt <- expr_data$expr788_filt
meta_39035   <- expr_data$meta_39035
meta_7888    <- expr_data$meta_7888

cat("✅ Datos cargados:\n")
cat(sprintf("   - DEG_age: %d genes\n", nrow(DEG_age_full)))
cat(sprintf("   - DEG_7888: %d genes\n", nrow(DEG_7888_full)))
cat(sprintf("   - Firma CORE: %d genes\n", length(core_signature)))
cat(sprintf("   - Expresión GSE39035: %d genes × %d muestras\n", 
            nrow(expr390_filt), ncol(expr390_filt)))
cat(sprintf("   - Expresión GSE7888: %d genes × %d muestras\n",
            nrow(expr788_filt), ncol(expr788_filt)))

# ============================================================================
# 2. CONSTRUCCIÓN DEL CATÁLOGO INMUNOMODULADOR
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("CONSTRUYENDO CATÁLOGO DE GENES INMUNOMODULADORES\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# 2.1 Genes de literatura (curados manualmente)
genes_lit <- c(
  # Checkpoints inmunosupresores
  "CD274", "PDCD1LG2", "IDO1", "IDO2", "HAVCR2", "LAG3", "CTLA4",
  "CD80", "CD86", "BTLA", "VSIR", "TNFRSF14",
  
  # HLA
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
  "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
  
  # Citoquinas inmunosupresoras
  "IL10", "TGFB1", "TGFB2", "TGFB3", "IL35", "IL27",
  
  # Citoquinas pro-inflamatorias (SASP)
  "IL6", "IL8", "CXCL8", "IL1A", "IL1B", "TNF",
  "CCL2", "CCL5", "CCL20", "CXCL1", "CXCL2", "CXCL10",
  
  # Enzimas inmunomoduladoras
  "ARG1", "ARG2", "NOS2", "PTGS2", "ENTPD1", "NT5E",
  
  # Factores de crecimiento
  "VEGFA", "HGF", "FGF2", "PDGFA", "PDGFB",
  
  # Moléculas adhesión
  "ICAM1", "VCAM1", "SELE", "SELP",
  
  # Factores transcripción relacionados
  "FOXP3", "STAT1", "STAT3", "STAT6", "NFKB1", "NFKB2", "RELA", "RELB"
)

cat(sprintf("📚 Genes curados de literatura: %d\n", length(genes_lit)))

# 2.2 Gene Ontology: Immune system process
cat("\n🔍 Descargando genes de Gene Ontology (GO:0002376)...\n")

go_immune <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = "GO:0002376",
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "GOALL"
) %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::pull(SYMBOL) %>%
  unique()

cat(sprintf("   GO immune genes: %d\n", length(go_immune)))

# 2.3 MSigDB C7: Immunologic signatures
cat("\n🔍 Descargando MSigDB C7 (Immunologic signatures)...\n")

msig_c7 <- msigdbr(species = "Homo sapiens", category = "C7")
genes_c7 <- unique(msig_c7$gene_symbol)

cat(sprintf("   MSigDB C7 genes: %d\n", length(genes_c7)))

# 2.4 Catálogo combinado
genes_inmune <- unique(c(genes_lit, go_immune, genes_c7))

cat("\n═══════════════════════════════════════════════════════════\n")
cat(sprintf("✅ CATÁLOGO FINAL: %d genes inmunomoduladores únicos\n", 
            length(genes_inmune)))
cat("═══════════════════════════════════════════════════════════\n")

# ============================================================================
# 3. INTERSECCIÓN: DEGs ∩ CATÁLOGO INMUNE
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("INTERSECCIÓN: DEGs × CATÁLOGO INMUNE\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# 3.1 Dataset Edad (GSE39035)
immu_age_sig <- DEG_age_full %>%
  dplyr::filter(Gene %in% genes_inmune, 
                adj.P.Val < 0.05, 
                abs(logFC) >= 1) %>%
  dplyr::arrange(adj.P.Val)

cat(sprintf("📊 GSE39035 (Edad):\n"))
cat(sprintf("   - DEGs inmunes totales: %d\n", nrow(immu_age_sig)))
cat(sprintf("   - Upregulated: %d\n", sum(immu_age_sig$logFC > 0)))
cat(sprintf("   - Downregulated: %d\n", sum(immu_age_sig$logFC < 0)))

# 3.2 Dataset Replicativa (GSE7888)
immu_7888_sig <- DEG_7888_full %>%
  dplyr::filter(Gene %in% genes_inmune,
                adj.P.Val < 0.05,
                abs(logFC) >= 1) %>%
  dplyr::arrange(adj.P.Val)

cat(sprintf("\n📊 GSE7888 (Replicativa):\n"))
cat(sprintf("   - DEGs inmunes totales: %d\n", nrow(immu_7888_sig)))
cat(sprintf("   - Upregulated: %d\n", sum(immu_7888_sig$logFC > 0)))
cat(sprintf("   - Downregulated: %d\n", sum(immu_7888_sig$logFC < 0)))

# 3.3 Verificación de genes clave
key_genes <- c("CD274", "IDO1", "IL6", "CXCL8", "CCL2", "TGFB1", 
               "HLA-A", "HLA-G", "ICAM1", "VCAM1")

cat("\n🔬 VERIFICACIÓN DE GENES CLAVE:\n")
cat("─────────────────────────────────────────────────────────\n")

for(gene in key_genes) {
  in_age <- gene %in% immu_age_sig$Gene
  in_rep <- gene %in% immu_7888_sig$Gene
  
  lfc_age <- ifelse(in_age, 
                    sprintf("%.2f", immu_age_sig$logFC[immu_age_sig$Gene == gene]),
                    "NS")
  lfc_rep <- ifelse(in_rep,
                    sprintf("%.2f", immu_7888_sig$logFC[immu_7888_sig$Gene == gene]),
                    "NS")
  
  cat(sprintf("%-10s | Edad: %6s | Replicativa: %6s\n", gene, lfc_age, lfc_rep))
}

# ============================================================================
# 4. ANÁLISIS CUANTITATIVO DE OVERLAP (MEJORA #1)
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("ANÁLISIS CUANTITATIVO DE CONVERGENCIA (Test Hipergeométrico)\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Universo: genes inmunes medidos en ambos datasets
genes_universe_390 <- intersect(genes_inmune, rownames(expr390_filt))
genes_universe_788 <- intersect(genes_inmune, rownames(expr788_filt))
genes_universe <- intersect(genes_universe_390, genes_universe_788)

N <- length(genes_universe)
M <- nrow(immu_age_sig)
n <- nrow(immu_7888_sig)
k <- length(intersect(immu_age_sig$Gene, immu_7888_sig$Gene))

# Test hipergeométrico
p_overlap <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
expected_overlap <- (M * n) / N
fold_enrichment <- k / expected_overlap

cat(sprintf("Universo (N):               %d genes inmunes\n", N))
cat(sprintf("DEGs Edad (M):              %d genes\n", M))
cat(sprintf("DEGs Replicativa (n):       %d genes\n", n))
cat(sprintf("Overlap observado (k):      %d genes\n", k))
cat(sprintf("Overlap esperado (azar):    %.1f genes\n", expected_overlap))
cat(sprintf("Fold-enrichment:            %.2fx\n", fold_enrichment))
cat(sprintf("P-valor (hipergeométrico):  %.2e\n", p_overlap))

if(p_overlap < 0.001) {
  cat("\n✅ CONCLUSIÓN: Overlap ALTAMENTE SIGNIFICATIVO (p < 0.001)\n")
  cat("   → Edad y Replicativa convergen en genes inmunes.\n")
  cat("   → Validación robusta de la hipótesis de convergencia funcional.\n")
} else if(p_overlap < 0.05) {
  cat("\n✅ CONCLUSIÓN: Overlap significativo (p < 0.05)\n")
  cat("   → Evidencia de convergencia parcial.\n")
} else {
  cat("\n⚠️ CONCLUSIÓN: Overlap NO significativo (p > 0.05)\n")
  cat("   → Edad y Replicativa afectan genes inmunes distintos.\n")
  cat("   → Heterogeneidad en mecanismos de senescencia.\n")
}

# Genes en overlap
genes_overlap <- intersect(immu_age_sig$Gene, immu_7888_sig$Gene)

cat(sprintf("\n📋 Genes en overlap (n=%d):\n", length(genes_overlap)))
if(length(genes_overlap) > 0 & length(genes_overlap) <= 20) {
  cat(paste(genes_overlap, collapse = ", "), "\n")
} else if(length(genes_overlap) > 20) {
  cat(paste(genes_overlap[1:20], collapse = ", "), "...\n")
}

# ============================================================================
# 5. ENRIQUECIMIENTO FUNCIONAL: ORA (GO + KEGG)
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("ENRIQUECIMIENTO FUNCIONAL: ORA (Over-Representation Analysis)\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# 5.1 Convertir a ENTREZ IDs
cat("🔄 Convirtiendo símbolos a ENTREZ IDs...\n")

symbol2entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(c(immu_age_sig$Gene, immu_7888_sig$Gene)),
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "SYMBOL"
) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

immu_age_entrez <- symbol2entrez %>%
  dplyr::filter(SYMBOL %in% immu_age_sig$Gene) %>%
  dplyr::pull(ENTREZID)

immu_7888_entrez <- symbol2entrez %>%
  dplyr::filter(SYMBOL %in% immu_7888_sig$Gene) %>%
  dplyr::pull(ENTREZID)

# Universo
universe_390_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = genes_universe_390,
  columns = "ENTREZID",
  keytype = "SYMBOL"
) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

universe_788_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = genes_universe_788,
  columns = "ENTREZID",
  keytype = "SYMBOL"
) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

cat(sprintf("   - Edad: %d → %d ENTREZ IDs\n", 
            nrow(immu_age_sig), length(immu_age_entrez)))
cat(sprintf("   - Replicativa: %d → %d ENTREZ IDs\n",
            nrow(immu_7888_sig), length(immu_7888_entrez)))

# 5.2 GO Biological Process - Edad
cat("\n🔍 GO Biological Process (Edad)...\n")

ego_age <- enrichGO(
  gene = immu_age_entrez,
  universe = universe_390_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 300,
  readable = TRUE
)

if(!is.null(ego_age) && nrow(ego_age@result) > 0) {
  cat(sprintf("   ✅ Términos enriquecidos: %d\n", nrow(ego_age@result)))
} else {
  cat("   ⚠️ No se encontraron términos enriquecidos\n")
}

# 5.3 GO Biological Process - Replicativa
cat("\n🔍 GO Biological Process (Replicativa)...\n")

ego_7888 <- enrichGO(
  gene = immu_7888_entrez,
  universe = universe_788_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 300,
  readable = TRUE
)

if(!is.null(ego_7888) && nrow(ego_7888@result) > 0) {
  cat(sprintf("   ✅ Términos enriquecidos: %d\n", nrow(ego_7888@result)))
} else {
  cat("   ⚠️ No se encontraron términos enriquecidos\n")
}

# 5.4 KEGG Pathways - Edad
cat("\n🔍 KEGG Pathways (Edad)...\n")

ek_age <- enrichKEGG(
  gene = immu_age_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 300
)

if(!is.null(ek_age) && nrow(ek_age@result) > 0) {
  ek_age <- setReadable(ek_age, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  cat(sprintf("   ✅ Pathways enriquecidos: %d\n", nrow(ek_age@result)))
} else {
  cat("   ⚠️ No se encontraron pathways enriquecidos\n")
}

# 5.5 KEGG Pathways - Replicativa
cat("\n🔍 KEGG Pathways (Replicativa)...\n")

ek_7888 <- enrichKEGG(
  gene = immu_7888_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 300
)

if(!is.null(ek_7888) && nrow(ek_7888@result) > 0) {
  ek_7888 <- setReadable(ek_7888, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  cat(sprintf("   ✅ Pathways enriquecidos: %d\n", nrow(ek_7888@result)))
} else {
  cat("   ⚠️ No se encontraron pathways enriquecidos\n")
}

# ============================================================================
# 6. GSEA: HALLMARK GENE SETS
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("GSEA: HALLMARK GENE SETS\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# 6.1 Descargar Hallmark
cat("🔄 Descargando Hallmark gene sets...\n")

msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

cat(sprintf("   ✅ Hallmark: %d gene sets\n", 
            length(unique(msig_h$gs_name))))

# 6.2 Ranking - Edad
cat("\n🔄 Construyendo ranking (Edad)...\n")

rank_age <- DEG_age_full %>%
  dplyr::mutate(score = sign(logFC) * (-log10(P.Value))) %>%
  dplyr::arrange(desc(score))

# Convertir a ENTREZ
rank_age_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rank_age$Gene,
  columns = "ENTREZID",
  keytype = "SYMBOL"
) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::left_join(rank_age, by = c("SYMBOL" = "Gene")) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_max(abs(score), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(score))

rank_age_vec <- rank_age_entrez$score
names(rank_age_vec) <- rank_age_entrez$ENTREZID

cat(sprintf("   ✅ Ranking: %d genes\n", length(rank_age_vec)))

# 6.3 GSEA - Edad
cat("\n🔍 Ejecutando GSEA (Edad)...\n")

gsea_age_h <- GSEA(
  geneList = rank_age_vec,
  TERM2GENE = msig_h,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 15,
  maxGSSize = 500,
  seed = 123,
  verbose = FALSE
)

if(!is.null(gsea_age_h) && nrow(gsea_age_h@result) > 0) {
  n_sig <- sum(gsea_age_h@result$p.adjust < 0.25)
  cat(sprintf("   ✅ Gene sets totales: %d\n", nrow(gsea_age_h@result)))
  cat(sprintf("   ✅ Significativos (FDR<0.25): %d\n", n_sig))
} else {
  cat("   ⚠️ No se encontraron gene sets enriquecidos\n")
}

# 6.4 Ranking - Replicativa
cat("\n🔄 Construyendo ranking (Replicativa)...\n")

rank_7888 <- DEG_7888_full %>%
  dplyr::mutate(score = sign(logFC) * (-log10(P.Value))) %>%
  dplyr::arrange(desc(score))

rank_7888_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rank_7888$Gene,
  columns = "ENTREZID",
  keytype = "SYMBOL"
) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::left_join(rank_7888, by = c("SYMBOL" = "Gene")) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_max(abs(score), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(score))

rank_7888_vec <- rank_7888_entrez$score
names(rank_7888_vec) <- rank_7888_entrez$ENTREZID

cat(sprintf("   ✅ Ranking: %d genes\n", length(rank_7888_vec)))

# 6.5 GSEA - Replicativa
cat("\n🔍 Ejecutando GSEA (Replicativa)...\n")

gsea_7888_h <- GSEA(
  geneList = rank_7888_vec,
  TERM2GENE = msig_h,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 15,
  maxGSSize = 500,
  seed = 123,
  verbose = FALSE
)

if(!is.null(gsea_7888_h) && nrow(gsea_7888_h@result) > 0) {
  n_sig <- sum(gsea_7888_h@result$p.adjust < 0.25)
  cat(sprintf("   ✅ Gene sets totales: %d\n", nrow(gsea_7888_h@result)))
  cat(sprintf("   ✅ Significativos (FDR<0.25): %d\n", n_sig))
} else {
  cat("   ⚠️ No se encontraron gene sets enriquecidos\n")
}

# 6.6 Vías clave esperadas
hallmark_expected <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING"
)

cat("\n🔬 VERIFICACIÓN DE VÍAS CLAVE (GSEA Edad):\n")
cat("─────────────────────────────────────────────────────────\n")

if(!is.null(gsea_age_h) && nrow(gsea_age_h@result) > 0) {
  for(pathway in hallmark_expected) {
    if(pathway %in% gsea_age_h@result$ID) {
      res <- gsea_age_h@result[gsea_age_h@result$ID == pathway, ]
      cat(sprintf("%-40s | NES: %6.3f | FDR: %.3f\n",
                  gsub("HALLMARK_", "", pathway),
                  res$NES,
                  res$p.adjust))
    } else {
      cat(sprintf("%-40s | No enriquecido\n", 
                  gsub("HALLMARK_", "", pathway)))
    }
  }
}

# ============================================================================
# 7. ANÁLISIS DE REGULADORES: TF Y miRNA
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("ANÁLISIS DE REGULADORES MAESTROS\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# 7.1 Factores de Transcripción (MSigDB C3:TFT)
cat("🔄 Descargando interacciones TF-Gene (MSigDB C3:TFT:GTRD)...\n")

msig_c3_tft <- msigdbr(species = "Homo sapiens", 
                       category = "C3", 
                       subcategory = "TFT:GTRD")

tf_network <- msig_c3_tft %>%
  dplyr::select(TF = gs_name, Target = gene_symbol) %>%
  dplyr::mutate(TF = gsub("_TARGET_GENES", "", TF))

cat(sprintf("   ✅ Red TF-Gene: %d TFs × %d genes\n",
            length(unique(tf_network$TF)),
            length(unique(tf_network$Target))))

# TFs regulando genes inmunes (Edad)
cat("\n🔍 TFs regulando genes inmunes (Edad)...\n")

tf_age <- tf_network %>%
  dplyr::filter(Target %in% immu_age_sig$Gene) %>%
  dplyr::count(TF, name = "n_targets") %>%
  dplyr::filter(n_targets >= 3) %>%
  dplyr::arrange(desc(n_targets))

cat(sprintf("   ✅ TFs identificados: %d (≥3 targets)\n", nrow(tf_age)))

if(nrow(tf_age) > 0) {
  cat("\n   Top 10 TFs:\n")
  print(head(tf_age, 10))
}

# TFs regulando genes inmunes (Replicativa)
cat("\n🔍 TFs regulando genes inmunes (Replicativa)...\n")

tf_7888 <- tf_network %>%
  dplyr::filter(Target %in% immu_7888_sig$Gene) %>%
  dplyr::count(TF, name = "n_targets") %>%
  dplyr::filter(n_targets >= 3) %>%
  dplyr::arrange(desc(n_targets))

cat(sprintf("   ✅ TFs identificados: %d (≥3 targets)\n", nrow(tf_7888)))

# Verificación de TFs esperados
expected_tfs <- c("RELA", "RELB", "NFKB1", "NFKB2", "STAT1", "STAT3", 
                  "IRF1", "IRF9", "TP53", "JUN", "FOS")

cat("\n🔬 VERIFICACIÓN DE TFs ESPERADOS (Edad):\n")
cat("─────────────────────────────────────────────────────────\n")

for(tf in expected_tfs) {
  if(tf %in% tf_age$TF) {
    n <- tf_age$n_targets[tf_age$TF == tf]
    cat(sprintf("%-10s | ✅ %d targets\n", tf, n))
  } else {
    cat(sprintf("%-10s | ⚠️ No detectado (< 3 targets)\n", tf))
  }
}

# 7.2 miRNAs (MSigDB C3:MIR)
cat("\n🔄 Descargando interacciones miRNA-Gene (MSigDB C3:MIR)...\n")

msig_mir <- msigdbr(species = "Homo sapiens",
                    category = "C3",
                    subcategory = "MIR:MIRDB")

mirna_network <- msig_mir %>%
  dplyr::select(miRNA = gs_name, Target = gene_symbol) %>%
  dplyr::mutate(miRNA = gsub("_.*", "", miRNA))

cat(sprintf("   ✅ Red miRNA-Gene: %d miRNAs × %d genes\n",
            length(unique(mirna_network$miRNA)),
            length(unique(mirna_network$Target))))

# miRNAs regulando genes inmunes (Edad)
cat("\n🔍 miRNAs regulando genes inmunes (Edad)...\n")

mir_age <- mirna_network %>%
  dplyr::filter(Target %in% immu_age_sig$Gene) %>%
  dplyr::count(miRNA, name = "n_targets") %>%
  dplyr::filter(n_targets >= 3) %>%
  dplyr::arrange(desc(n_targets))

cat(sprintf("   ✅ miRNAs identificados: %d (≥3 targets)\n", nrow(mir_age)))

if(nrow(mir_age) > 0) {
  cat("\n   Top 10 miRNAs:\n")
  print(head(mir_age, 10))
}

# miRNAs regulando genes inmunes (Replicativa)
cat("\n🔍 miRNAs regulando genes inmunes (Replicativa)...\n")

mir_7888 <- mirna_network %>%
  dplyr::filter(Target %in% immu_7888_sig$Gene) %>%
  dplyr::count(miRNA, name = "n_targets") %>%
  dplyr::filter(n_targets >= 3) %>%
  dplyr::arrange(desc(n_targets))

cat(sprintf("   ✅ miRNAs identificados: %d (≥3 targets)\n", nrow(mir_7888)))

# Verificación de miRNAs esperados
expected_mirs <- c("MIR146A", "MIR155", "MIR34A", "MIR21", "MIR29A")

cat("\n🔬 VERIFICACIÓN DE miRNAs ESPERADOS (Edad):\n")
cat("─────────────────────────────────────────────────────────\n")

for(mir in expected_mirs) {
  if(mir %in% mir_age$miRNA) {
    n <- mir_age$n_targets[mir_age$miRNA == mir]
    cat(sprintf("%-10s | ✅ %d targets\n", mir, n))
  } else {
    cat(sprintf("%-10s | ⚠️ No detectado (< 3 targets)\n", mir))
  }
}

# ============================================================================
# 8. RED REGULATORIA TF-GENE-miRNA + ANÁLISIS DE CENTRALIDAD (MEJORA #2)
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("RED REGULATORIA TF-GENE-miRNA CON ANÁLISIS DE CENTRALIDAD\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# 8.1 Seleccionar top reguladores
top_tfs <- unique(c(
  tf_age$TF[1:min(10, nrow(tf_age))],
  tf_7888$TF[1:min(10, nrow(tf_7888))]
))

top_mirs <- unique(c(
  mir_age$miRNA[1:min(10, nrow(mir_age))],
  mir_7888$miRNA[1:min(10, nrow(mir_7888))]
))

# Genes core inmunes (top por prioridad)
core_immune_genes <- unique(c(
  immu_age_sig$Gene[1:min(20, nrow(immu_age_sig))],
  immu_7888_sig$Gene[1:min(20, nrow(immu_7888_sig))],
  genes_overlap
))

cat(sprintf("🔧 Componentes de la red:\n"))
cat(sprintf("   - TFs: %d\n", length(top_tfs)))
cat(sprintf("   - miRNAs: %d\n", length(top_mirs)))
cat(sprintf("   - Genes: %d\n", length(core_immune_genes)))

# 8.2 Construir edges
edges_tf <- tf_network %>%
  dplyr::filter(TF %in% top_tfs, Target %in% core_immune_genes) %>%
  dplyr::select(from = TF, to = Target) %>%
  dplyr::mutate(Type = "TF", Weight = 2)

edges_mir <- mirna_network %>%
  dplyr::filter(miRNA %in% top_mirs, Target %in% core_immune_genes) %>%
  dplyr::select(from = miRNA, to = Target) %>%
  dplyr::mutate(Type = "miRNA", Weight = 1)

edges_all <- dplyr::bind_rows(edges_tf, edges_mir) %>%
  dplyr::distinct()

cat(sprintf("\n🔗 Edges construidos:\n"))
cat(sprintf("   - TF→Gene: %d\n", nrow(edges_tf)))
cat(sprintf("   - miRNA→Gene: %d\n", nrow(edges_mir)))
cat(sprintf("   - Total: %d\n", nrow(edges_all)))

# 8.3 Construir nodos
all_nodes <- unique(c(edges_all$from, edges_all$to))

nodes <- data.frame(
  name = all_nodes,
  node_type = ifelse(all_nodes %in% top_tfs, "TF",
                     ifelse(all_nodes %in% top_mirs, "miRNA", "Gene")),
  stringsAsFactors = FALSE
)

# Añadir logFC promedio para genes
gene_logfc <- bind_rows(
  immu_age_sig %>% dplyr::select(Gene, logFC),
  immu_7888_sig %>% dplyr::select(Gene, logFC)
) %>%
  group_by(Gene) %>%
  summarise(logFC = mean(logFC), .groups = "drop")

nodes <- nodes %>%
  left_join(gene_logfc, by = c("name" = "Gene"))

# 8.4 Crear grafo
g <- graph_from_data_frame(edges_all, vertices = nodes, directed = TRUE)

cat(sprintf("\n📊 Red construida:\n"))
cat(sprintf("   - Nodos: %d (%d TFs, %d miRNAs, %d Genes)\n",
            vcount(g),
            sum(V(g)$node_type == "TF"),
            sum(V(g)$node_type == "miRNA"),
            sum(V(g)$node_type == "Gene")))
cat(sprintf("   - Edges: %d\n", ecount(g)))

# 8.5 CALCULAR MÉTRICAS DE CENTRALIDAD (MEJORA #2)
cat("\n🔍 Calculando métricas de centralidad...\n")

V(g)$betweenness <- betweenness(g, directed = TRUE)
V(g)$closeness <- closeness(g, mode = "out")
V(g)$pagerank <- page_rank(g)$vector

# Extraer métricas
nodes_centrality <- data.frame(
  name = V(g)$name,
  type = V(g)$node_type,
  betweenness = V(g)$betweenness,
  closeness = V(g)$closeness,
  pagerank = V(g)$pagerank,
  stringsAsFactors = FALSE
)

# 8.6 RANKING DE REGULADORES POR PAGERANK
top_regulators <- nodes_centrality %>%
  dplyr::filter(type %in% c("TF", "miRNA")) %>%
  dplyr::arrange(desc(pagerank)) %>%
  dplyr::slice(1:15)

# ============================================================================
# 9. CONSTRUCCIÓN DEL PANEL MSC-IMMUNOSCORE
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("CONSTRUCCIÓN DEL PANEL MSC-IMMUNOSCORE\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# 9.1 Genes regulados (por top TFs/miRNAs)
regulated_genes <- unique(edges_all$to)

cat(sprintf("🔧 Genes regulados por top TFs/miRNAs: %d\n", 
            length(regulated_genes)))

# 9.2 Scoring multi-criterio - Edad
cat("\n🔍 Calculando Priority Score (Edad)...\n")

panel_age <- immu_age_sig %>%
  dplyr::mutate(
    LitSupport = ifelse(Gene %in% genes_lit, 10, 0),
    RegScore = ifelse(Gene %in% regulated_genes, 5, 0),
    PriorityScore = 2 * abs(logFC) + 
      1 * (-log10(adj.P.Val)) + 
      LitSupport + 
      RegScore
  ) %>%
  dplyr::arrange(desc(PriorityScore)) %>%
  dplyr::slice(1:30)

cat(sprintf("   ✅ Top 30 genes (Edad)\n"))

# 9.3 Scoring multi-criterio - Replicativa
cat("\n🔍 Calculando Priority Score (Replicativa)...\n")

panel_7888 <- immu_7888_sig %>%
  dplyr::mutate(
    LitSupport = ifelse(Gene %in% genes_lit, 10, 0),
    RegScore = ifelse(Gene %in% regulated_genes, 5, 0),
    PriorityScore = 2 * abs(logFC) + 
      1 * (-log10(adj.P.Val)) + 
      LitSupport + 
      RegScore
  ) %>%
  dplyr::arrange(desc(PriorityScore)) %>%
  dplyr::slice(1:30)

cat(sprintf("   ✅ Top 30 genes (Replicativa)\n"))

# 9.4 Panel CORE (intersección + top de cada uno)
core_panel <- intersect(panel_age$Gene, panel_7888$Gene)

# Añadir top únicos de cada dataset
genes_panel <- unique(c(
  core_panel,
  panel_age$Gene[1:min(15, nrow(panel_age))],
  panel_7888$Gene[1:min(15, nrow(panel_7888))]
))

cat("\n═══════════════════════════════════════════════════════════\n")
cat(sprintf("✅ PANEL MSC-IMMUNOSCORE: %d genes\n", length(genes_panel)))
cat("═══════════════════════════════════════════════════════════\n")
cat(sprintf("   - Genes core (overlap): %d\n", length(core_panel)))
cat(sprintf("   - Genes adicionales: %d\n", 
            length(genes_panel) - length(core_panel)))

# 9.5 Clasificación funcional
cat("\n🔍 Clasificando genes por función...\n")

# Definir categorías funcionales
immunosuppressive <- c("CD274", "PDCD1LG2", "IDO1", "IDO2", "HAVCR2",
                       "HLA-G", "HLA-E", "IL10", "TGFB1", "TGFB2", "TGFB3",
                       "ARG1", "ARG2", "ENTPD1", "NT5E", "FOXP3")

proinflammatory <- c("IL6", "IL8", "CXCL8", "IL1A", "IL1B", "TNF",
                     "CCL2", "CCL5", "CCL20", "CXCL1", "CXCL2", "CXCL10",
                     "ICAM1", "VCAM1", "NOS2", "PTGS2")

# Asignar pesos
gene_weights <- rep(0, length(genes_panel))
names(gene_weights) <- genes_panel

for(gene in genes_panel) {
  if(gene %in% immunosuppressive) {
    gene_weights[gene] <- -1  # Inmunosupresor (pérdida = malo)
  } else if(gene %in% proinflammatory) {
    gene_weights[gene] <- +1  # Pro-inflamatorio (ganancia = malo)
  }
}

n_immsupp <- sum(gene_weights == -1)
n_proinf <- sum(gene_weights == +1)
n_neutral <- sum(gene_weights == 0)

cat(sprintf("   ✅ Clasificación:\n"))
cat(sprintf("      - Inmunosupresores (w=-1): %d genes\n", n_immsupp))
cat(sprintf("      - Pro-inflamatorios (w=+1): %d genes\n", n_proinf))
cat(sprintf("      - Neutros/Otros (w=0): %d genes\n", n_neutral))

# ============================================================================
# 10. CÁLCULO DEL MSC-IMMUNOSCORE
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("CÁLCULO DEL MSC-IMMUNOSCORE\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Función de cálculo
calculate_immunoscore <- function(expr_mat, genes, weights) {
  genes_present <- intersect(genes, rownames(expr_mat))
  
  if(length(genes_present) == 0) {
    warning("No genes del panel presentes en matriz de expresión")
    return(rep(NA, ncol(expr_mat)))
  }
  
  mat_sub <- expr_mat[genes_present, , drop = FALSE]
  mat_z <- t(scale(t(mat_sub)))
  weights_sub <- weights[genes_present]
  mat_weighted <- sweep(mat_z, 1, weights_sub, "*")
  scores <- colMeans(mat_weighted, na.rm = TRUE)
  
  return(scores)
}

# Calcular scores
cat("🔍 Calculando ImmunoScore para GSE39035...\n")
score_age <- calculate_immunoscore(expr390_filt, genes_panel, gene_weights)

cat("🔍 Calculando ImmunoScore para GSE7888...\n")
score_7888 <- calculate_immunoscore(expr788_filt, genes_panel, gene_weights)

# Añadir a metadata
meta_39035$ImmunoScore <- score_age
meta_7888$ImmunoScore <- score_7888

cat("\n✅ ImmunoScore calculado:\n")
cat(sprintf("   - GSE39035: Rango [%.3f, %.3f]\n",
            min(score_age, na.rm = TRUE),
            max(score_age, na.rm = TRUE)))
cat(sprintf("   - GSE7888: Rango [%.3f, %.3f]\n",
            min(score_7888, na.rm = TRUE),
            max(score_7888, na.rm = TRUE)))

# Comparación por grupos
cat("\n📊 ImmunoScore por grupo (GSE39035):\n")
cat("─────────────────────────────────────────────────────────\n")

score_summary_age <- meta_39035 %>%
  group_by(GroupEdad) %>%
  summarise(
    n = n(),
    Mean = mean(ImmunoScore, na.rm = TRUE),
    SD = sd(ImmunoScore, na.rm = TRUE),
    Min = min(ImmunoScore, na.rm = TRUE),
    Max = max(ImmunoScore, na.rm = TRUE),
    .groups = "drop"
  )

print(score_summary_age)

# Test t
if("Young" %in% meta_39035$GroupEdad && "Old" %in% meta_39035$GroupEdad) {
  test_age <- t.test(
    ImmunoScore ~ GroupEdad,
    data = meta_39035
  )
  
  cat(sprintf("\n📊 Test t (Young vs Old):\n"))
  cat(sprintf("   t = %.3f, p = %.4f\n", 
              test_age$statistic, test_age$p.value))
  
  if(test_age$p.value < 0.05) {
    cat("   ✅ Diferencia significativa entre grupos\n")
  } else {
    cat("   ⚠️ Diferencia NO significativa\n")
  }
}

cat("\n📊 ImmunoScore por pasaje (GSE7888):\n")
cat("─────────────────────────────────────────────────────────\n")

score_summary_7888 <- meta_7888 %>%
  group_by(GroupPassage) %>%
  summarise(
    n = n(),
    Mean = mean(ImmunoScore, na.rm = TRUE),
    SD = sd(ImmunoScore, na.rm = TRUE),
    Min = min(ImmunoScore, na.rm = TRUE),
    Max = max(ImmunoScore, na.rm = TRUE),
    .groups = "drop"
  )

print(score_summary_7888)

# ============================================================================
# 11. VALIDACIÓN ROBUSTA: ROC CURVES (MEJORA #3)
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("VALIDACIÓN ROBUSTA: ROC CURVES\n")
cat("═══════════════════════════════════════════════════════════\n\n")

cat("🔍 Comparando poder predictivo: ImmunoScore vs variables clínicas\n\n")

# 11.1 Preparar datos (GSE39035)
df_roc_age <- meta_39035 %>%
  dplyr::mutate(
    Truth = ifelse(GroupEdad == "Old", 1, 0),
    Age_numeric = as.numeric(Age_years),
    Passage_numeric = as.numeric(Passage)
  ) %>%
  dplyr::filter(!is.na(ImmunoScore), !is.na(Age_numeric))

# 11.2 ROC curves individuales
cat("📈 Calculando ROC curves...\n")

roc_score <- roc(Truth ~ ImmunoScore, data = df_roc_age, quiet = TRUE)
roc_age_var <- roc(Truth ~ Age_numeric, data = df_roc_age, quiet = TRUE)
roc_pass <- roc(Truth ~ Passage_numeric, data = df_roc_age, quiet = TRUE)

# 11.3 Modelo combinado (Edad + Pasaje)
model_combined <- glm(Truth ~ Age_numeric + Passage_numeric,
                      data = df_roc_age, family = binomial)
df_roc_age$Combined_pred <- predict(model_combined, type = "response")
roc_combined <- roc(Truth ~ Combined_pred, data = df_roc_age, quiet = TRUE)

# 11.4 Comparar AUCs
cat("\n═══════════════════════════════════════════════════════════\n")
cat("COMPARACIÓN DE PODER PREDICTIVO (ROC-AUC)\n")
cat("═══════════════════════════════════════════════════════════\n\n")

cat(sprintf("ImmunoScore (panel genes):  AUC = %.3f [%.3f-%.3f]\n",
            auc(roc_score),
            ci(roc_score)[1], ci(roc_score)[3]))
cat(sprintf("Edad sola:                  AUC = %.3f [%.3f-%.3f]\n",
            auc(roc_age_var),
            ci(roc_age_var)[1], ci(roc_age_var)[3]))
cat(sprintf("Pasaje solo:                AUC = %.3f [%.3f-%.3f]\n",
            auc(roc_pass),
            ci(roc_pass)[1], ci(roc_pass)[3]))
cat(sprintf("Edad + Pasaje (combinado):  AUC = %.3f [%.3f-%.3f]\n",
            auc(roc_combined),
            ci(roc_combined)[1], ci(roc_combined)[3]))

# 11.5 Tests estadísticos (DeLong)
cat("\n📊 TESTS ESTADÍSTICOS (DeLong):\n")
cat("─────────────────────────────────────────────────────────\n")

test_vs_age <- roc.test(roc_score, roc_age_var, method = "delong")
test_vs_pass <- roc.test(roc_score, roc_pass, method = "delong")
test_vs_comb <- roc.test(roc_score, roc_combined, method = "delong")

cat(sprintf("ImmunoScore vs Edad:        p = %.4f %s\n",
            test_vs_age$p.value,
            ifelse(test_vs_age$p.value < 0.05, "(*)", "")))
cat(sprintf("ImmunoScore vs Pasaje:      p = %.4f %s\n",
            test_vs_pass$p.value,
            ifelse(test_vs_pass$p.value < 0.05, "(*)", "")))
cat(sprintf("ImmunoScore vs Edad+Pasaje: p = %.4f %s\n",
            test_vs_comb$p.value,
            ifelse(test_vs_comb$p.value < 0.05, "(*)", "")))

# 11.6 Interpretación
cat("\n💡 INTERPRETACIÓN:\n")

if(auc(roc_score) > auc(roc_combined) && test_vs_comb$p.value < 0.05) {
  cat("✅ ImmunoScore supera significativamente a variables clínicas.\n")
  cat("   → Panel transcriptómico aporta valor predictivo único.\n")
  cat("   → Validar en cohorte independiente (GSE35958).\n")
} else if(auc(roc_score) > auc(roc_age_var) && test_vs_age$p.value < 0.05) {
  cat("✅ ImmunoScore supera a edad sola.\n")
  cat("   → Pero no supera significativamente a modelo combinado.\n")
  cat("   → Sugerencia: Integrar ImmunoScore + clínicas.\n")
} else {
  cat("⚠️ ImmunoScore NO supera significativamente a variables clínicas.\n")
  cat("   → Panel requiere refinamiento:\n")
  cat("      • Aumentar tamaño de panel (feature selection)\n")
  cat("      • Validar con LASSO/Elastic Net\n")
  cat("      • Integrar con machine learning (Fase 3)\n")
}

# ============================================================================
# 12. ANÁLISIS DE RESCATE VIRTUAL (MEJORA #4)
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("ANÁLISIS DE RESCATE VIRTUAL (IN SILICO)\n")
cat("═══════════════════════════════════════════════════════════\n\n")

cat("🔬 Simulando overexpression de top TFs → Efecto en ImmunoScore\n\n")

# 12.1 Top TFs por PageRank
top_tfs_rescue <- top_regulators %>%
  dplyr::filter(type == "TF") %>%
  dplyr::slice(1:5) %>%
  dplyr::pull(name)

cat(sprintf("🎯 TFs a evaluar (top 5): %s\n",
            paste(top_tfs_rescue, collapse = ", ")))

# 12.2 Función de perturbación
simulate_overexpression <- function(tf_name, network, panel_genes, weights) {
  
  # Genes target de este TF
  targets <- network %>%
    dplyr::filter(from == tf_name, Type == "TF") %>%
    dplyr::pull(to)
  
  # Targets en nuestro panel
  targets_panel <- intersect(targets, panel_genes)
  
  if(length(targets_panel) == 0) {
    return(data.frame(
      TF = tf_name,
      n_targets_panel = 0,
      n_immunosupp_targets = 0,
      predicted_delta_score = NA,
      rescue_potential = NA,
      interpretation = "No targets en panel"
    ))
  }
  
  # Contar targets inmunosupresores (downregulados en senescencia)
  immunosupp_targets <- targets_panel[weights[targets_panel] == -1]
  n_immunosupp <- length(immunosupp_targets)
  
  # Simular rescate: upregular genes inmunosupresores
  # Lógica: TF activa targets → si son inmunosupresores, reduce penalización
  weights_perturbed <- weights
  
  if(n_immunosupp > 0) {
    # Reducir peso negativo (simular rescate parcial)
    for(target in immunosupp_targets) {
      weights_perturbed[target] <- weights[target] * 0.5  # Rescate 50%
    }
  }
  
  # Calcular delta score
  score_original <- mean(weights[panel_genes])
  score_perturbed <- mean(weights_perturbed[panel_genes])
  delta <- score_perturbed - score_original
  
  # Rescue potential: negativo de delta (queremos bajar score)
  rescue_score <- -delta
  
  # Interpretación
  if(rescue_score > 0.05) {
    interp <- "Alto potencial de rescate"
  } else if(rescue_score > 0.01) {
    interp <- "Moderado potencial"
  } else if(rescue_score > 0) {
    interp <- "Bajo potencial"
  } else {
    interp <- "Sin efecto de rescate"
  }
  
  return(data.frame(
    TF = tf_name,
    n_targets_panel = length(targets_panel),
    n_immunosupp_targets = n_immunosupp,
    predicted_delta_score = delta,
    rescue_potential = rescue_score,
    interpretation = interp,
    stringsAsFactors = FALSE
  ))
}

# 12.3 Evaluar top TFs
rescue_results <- lapply(top_tfs_rescue, function(tf) {
  simulate_overexpression(tf, edges_all, genes_panel, gene_weights)
}) %>% bind_rows() %>%
  dplyr::arrange(desc(rescue_potential))

# 12.4 Resultados
cat("\n═══════════════════════════════════════════════════════════\n")
cat("RESULTADOS: ANÁLISIS DE RESCATE VIRTUAL\n")
cat("═══════════════════════════════════════════════════════════\n\n")

rescue_display <- rescue_results %>%
  dplyr::mutate(
    Prioridad = row_number(),
    `Rescue Potential` = sprintf("%.4f", rescue_potential),
    `Delta Score` = sprintf("%.4f", predicted_delta_score)
  ) %>%
  dplyr::select(Prioridad, TF, 
                `Targets en panel` = n_targets_panel,
                `Targets inmunosup.` = n_immunosupp_targets,
                `Rescue Potential`,
                Interpretación = interpretation)

print(rescue_display)

cat("\n💡 INTERPRETACIÓN:\n")
cat("   - Rescue potential > 0.05: ALTA prioridad experimental\n")
cat("   - Rescue potential 0.01-0.05: MEDIA prioridad\n")
cat("   - Rescue potential < 0.01: BAJA prioridad\n")
cat(sprintf("\n🎯 TOP CANDIDATO: %s\n", rescue_results$TF[1]))
cat(sprintf("   - Targets totales: %d\n", 
            rescue_results$n_targets_panel[1]))
cat(sprintf("   - Targets inmunosupresores: %d\n",
            rescue_results$n_immunosupp_targets[1]))
cat(sprintf("   - Rescue potential: %.4f\n",
            rescue_results$rescue_potential[1]))

cat("\n📋 VALIDACIÓN EXPERIMENTAL SUGERIDA:\n")
cat("   1. Transfectar hMSC senescentes con vector de overexpression\n")
cat("   2. Medir expresión de genes inmunosupresores (qPCR/Western):\n")
cat("      • PD-L1 (CD274)\n")
cat("      • IDO1\n")
cat("      • HLA-G\n")
cat("   3. Ensayos funcionales:\n")
cat("      • MLR (Mixed Lymphocyte Reaction) → Supresión T cells\n")
cat("      • ELISA → IL-6, IL-10, TGF-β\n")
cat("   4. Calcular ImmunoScore post-intervención\n")

# ============================================================================
# 13. VALIDACIÓN CON LITERATURA (MEJORA #5)
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("VALIDACIÓN CON LITERATURA\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# 13.1 Tabla comparativa
literature_genes <- data.frame(
  Gene = c("CD274", "IDO1", "GATA2", "HLA-A", "HLA-G", 
           "IL6", "CXCL8", "CCL2", "TGFB1", "ICAM1"),
  Literatura = c(
    "Down en senescencia (Gao 2023)",
    "Down en senescencia (Gao 2023)",
    "Down 60%, regula PD-L1 (Gao 2023)",
    "Alterado (Yang 2020)",
    "Down en senescencia",
    "Up (SASP, Coppé 2008)",
    "Up (SASP, Coppé 2008)",
    "Up (SASP)",
    "Down (pérdida inmunosup.)",
    "Up (senescencia)"
  ),
  stringsAsFactors = FALSE
)

# Añadir nuestros datos
literature_genes$`Nuestros datos (Edad)` <- sapply(literature_genes$Gene, function(g) {
  if(g %in% immu_age_sig$Gene) {
    lfc <- immu_age_sig$logFC[immu_age_sig$Gene == g]
    direction <- ifelse(lfc > 0, "↑", "↓")
    return(sprintf("%s %.2f", direction, lfc))
  } else {
    return("NS")
  }
})

literature_genes$`Nuestros datos (Rep.)` <- sapply(literature_genes$Gene, function(g) {
  if(g %in% immu_7888_sig$Gene) {
    lfc <- immu_7888_sig$logFC[immu_7888_sig$Gene == g]
    direction <- ifelse(lfc > 0, "↑", "↓")
    return(sprintf("%s %.2f", direction, lfc))
  } else {
    return("NS")
  }
})

# Concordancia
literature_genes$Concordancia <- sapply(1:nrow(literature_genes), function(i) {
  lit <- literature_genes$Literatura[i]
  edad <- literature_genes$`Nuestros datos (Edad)`[i]
  rep <- literature_genes$`Nuestros datos (Rep.)`[i]
  
  if(edad == "NS" && rep == "NS") return("⚠️")
  
  expected_down <- grepl("Down", lit)
  expected_up <- grepl("Up", lit)
  
  edad_down <- grepl("↓", edad)
  edad_up <- grepl("↑", edad)
  rep_down <- grepl("↓", rep)
  rep_up <- grepl("↑", rep)
  
  if((expected_down && (edad_down || rep_down)) ||
     (expected_up && (edad_up || rep_up))) {
    return("✅")
  } else if(edad == "NS" && rep == "NS") {
    return("⚠️")
  } else {
    return("❌")
  }
})

cat("═══════════════════════════════════════════════════════════\n")
cat("TABLA: VALIDACIÓN CON LITERATURA\n")
cat("═══════════════════════════════════════════════════════════\n\n")

print(literature_genes, row.names = FALSE)

cat("\n📊 RESUMEN DE CONCORDANCIA:\n")
cat(sprintf("   ✅ Concordante: %d genes\n", 
            sum(literature_genes$Concordancia == "✅")))
cat(sprintf("   ❌ Discordante: %d genes\n",
            sum(literature_genes$Concordancia == "❌")))
cat(sprintf("   ⚠️ No significativo: %d genes\n",
            sum(literature_genes$Concordancia == "⚠️")))

concordance_rate <- sum(literature_genes$Concordancia == "✅") / 
  nrow(literature_genes) * 100

cat(sprintf("\n💡 Tasa de concordancia: %.1f%%\n", concordance_rate))

if(concordance_rate >= 70) {
  cat("   ✅ Alta concordancia con literatura publicada\n")
  cat("   → Validación robusta de nuestros hallazgos\n")
} else if(concordance_rate >= 50) {
  cat("   ⚠️ Concordancia moderada\n")
  cat("   → Revisar genes discordantes (diferencias dataset/metodología)\n")
} else {
  cat("   ❌ Baja concordancia\n")
  cat("   → Revisar pipeline o considerar heterogeneidad biológica\n")
}

# ============================================================================
# 14. EXPORTACIÓN DE DATOS PARA FASE 3
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("EXPORTACIÓN DE DATOS PARA FASE 3\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# 14.1 Lista completa de objetos
fase2_export <- list(
  # Genes inmunes
  immu_age_sig = immu_age_sig,
  immu_7888_sig = immu_7888_sig,
  genes_inmune_catalog = genes_inmune,
  genes_overlap = genes_overlap,
  
  # Enriquecimiento
  ego_age = ego_age,
  ego_7888 = ego_7888,
  ek_age = ek_age,
  ek_7888 = ek_7888,
  gsea_age_h = gsea_age_h,
  gsea_7888_h = gsea_7888_h,
  
  # Reguladores
  tf_age = tf_age,
  tf_7888 = tf_7888,
  mir_age = mir_age,
  mir_7888 = mir_7888,
  tf_network = tf_network,
  mirna_network = mirna_network,
  
  # Red regulatoria
  regulatory_network = g,
  nodes_centrality = nodes_centrality,
  top_regulators = top_regulators,
  
  # Panel MSC-ImmunoScore
  panel_age = panel_age,
  panel_7888 = panel_7888,
  genes_panel = genes_panel,
  gene_weights = gene_weights,
  core_panel = core_panel,
  
  # ImmunoScore
  score_age = score_age,
  score_7888 = score_7888,
  
  # Validación
  roc_results = list(
    roc_score = roc_score,
    roc_age = roc_age_var,
    roc_pass = roc_pass,
    roc_combined = roc_combined,
    test_vs_combined = test_vs_comb
  ),
  
  rescue_analysis = rescue_results,
  literature_validation = literature_genes,
  
  # Overlap estadístico
  overlap_stats = list(
    N = N,
    M = M,
    n = n,
    k = k,
    p_value = p_overlap,
    fold_enrichment = fold_enrichment
  ),
  
  # Metadata
  meta_39035 = meta_39035,
  meta_7888 = meta_7888,
  
  # Expresión
  expr390_filt = expr390_filt,
  expr788_filt = expr788_filt
)

# 14.2 Guardar
cat("💾 Guardando objetos...\n")

saveRDS(fase2_export, 
        "Fase2-Enrichment/results/FASE2_complete_export.rds")

cat("   ✅ FASE2_complete_export.rds\n")

# 14.3 Panel en CSV
cat("\n💾 Exportando panel MSC-ImmunoScore...\n")

panel_csv <- data.frame(
  Gene = names(gene_weights),
  Weight = gene_weights,
  Type = ifelse(gene_weights == -1, "Immunosuppressive",
                ifelse(gene_weights == 1, "Proinflammatory", "Neutral")),
  InCorePanelOverlap = names(gene_weights) %in% core_panel,
  stringsAsFactors = FALSE
)

write.csv(panel_csv,
          "Fase2-Enrichment/results/MSC_ImmunoScore_panel.csv",
          row.names = FALSE)

cat("   ✅ MSC_ImmunoScore_panel.csv\n")

# 14.4 Resumen de reguladores
cat("\n💾 Exportando top reguladores...\n")

write.csv(top_regulators,
          "Fase2-Enrichment/results/top_regulators_pagerank.csv",
          row.names = FALSE)

cat("   ✅ top_regulators_pagerank.csv\n")

# 14.5 Resultados de rescate
cat("\n💾 Exportando análisis de rescate...\n")

write.csv(rescue_results,
          "Fase2-Enrichment/results/rescue_analysis_results.csv",
          row.names = FALSE)

cat("   ✅ rescue_analysis_results.csv\n")

cat("\n═══════════════════════════════════════════════════════════\n")
cat("✅ EXPORTACIÓN COMPLETADA\n")
cat("═══════════════════════════════════════════════════════════\n\n")

cat("📁 Archivos generados:\n")
cat("   Fase2-Enrichment/results/\n")
cat("   ├── FASE2_complete_export.rds (objeto R completo)\n")
cat("   ├── MSC_ImmunoScore_panel.csv (panel de genes)\n")
cat("   ├── top_regulators_pagerank.csv (TFs/miRNAs prioritarios)\n")
cat("   └── rescue_analysis_results.csv (análisis de rescate)\n")

# ============================================================================
# 15. RESUMEN FINAL
# ============================================================================

cat("\n\n")
cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║                   FASE 2 COMPLETADA                        ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

cat("📊 RESUMEN DE RESULTADOS:\n")
cat("─────────────────────────────────────────────────────────────\n\n")

cat(sprintf("1. CATÁLOGO INMUNOMODULADOR:\n"))
cat(sprintf("   - Genes totales: %d\n", length(genes_inmune)))
cat(sprintf("   - Fuentes: Literatura + GO + MSigDB C7\n\n"))

cat(sprintf("2. DEGs INMUNES:\n"))
cat(sprintf("   - Edad: %d genes (Up: %d, Down: %d)\n",
            nrow(immu_age_sig),
            sum(immu_age_sig$logFC > 0),
            sum(immu_age_sig$logFC < 0)))
cat(sprintf("   - Replicativa: %d genes (Up: %d, Down: %d)\n",
            nrow(immu_7888_sig),
            sum(immu_7888_sig$logFC > 0),
            sum(immu_7888_sig$logFC < 0)))
cat(sprintf("   - Overlap: %d genes (p=%.2e, FE=%.2fx)\n\n",
            k, p_overlap, fold_enrichment))

if(!is.null(ego_age) && nrow(ego_age@result) > 0) {
  cat(sprintf("3. ENRIQUECIMIENTO FUNCIONAL:\n"))
  cat(sprintf("   - GO terms (Edad): %d\n", nrow(ego_age@result)))
  if(!is.null(ek_age) && nrow(ek_age@result) > 0) {
    cat(sprintf("   - KEGG pathways (Edad): %d\n", nrow(ek_age@result)))
  }
  if(!is.null(gsea_age_h) && nrow(gsea_age_h@result) > 0) {
    n_sig_gsea <- sum(gsea_age_h@result$p.adjust < 0.25)
    cat(sprintf("   - GSEA Hallmark (Edad): %d sig (FDR<0.25)\n\n", n_sig_gsea))
  }
}

cat(sprintf("4. REGULADORES MAESTROS:\n"))
cat(sprintf("   - TFs identificados: %d (≥3 targets)\n", nrow(tf_age)))
cat(sprintf("   - miRNAs identificados: %d (≥3 targets)\n", nrow(mir_age)))
cat(sprintf("   - Top regulador (PageRank): %s\n\n",
            top_regulators$name[1]))

cat(sprintf("5. PANEL MSC-IMMUNOSCORE:\n"))
cat(sprintf("   - Genes totales: %d\n", length(genes_panel)))
cat(sprintf("   - Inmunosupresores: %d\n", n_immsupp))
cat(sprintf("   - Pro-inflamatorios: %d\n", n_proinf))
cat(sprintf("   - ImmunoScore AUC: %.3f\n\n", auc(roc_score)))

cat(sprintf("6. ANÁLISIS DE RESCATE:\n"))
cat(sprintf("   - Top candidato TF: %s\n", rescue_results$TF[1]))
cat(sprintf("   - Rescue potential: %.4f\n", 
            rescue_results$rescue_potential[1]))
cat(sprintf("   - Targets inmunosup.: %d\n\n",
            rescue_results$n_immunosupp_targets[1]))

cat(sprintf("7. VALIDACIÓN CON LITERATURA:\n"))
cat(sprintf("   - Genes evaluados: %d\n", nrow(literature_genes)))
cat(sprintf("   - Tasa concordancia: %.1f%%\n\n", concordance_rate))

cat("═════════════════════════════════════════════════════════════\n")
cat("✅ FASE 2 EJECUTADA CON ÉXITO\n")
cat("═════════════════════════════════════════════════════════════\n\n")

cat("🚀 PRÓXIMOS PASOS (FASE 3 - Machine Learning):\n\n")
cat("1. Cargar datos:\n")
cat("   fase2 <- readRDS('Fase2-Enrichment/results/FASE2_complete_export.rds')\n\n")
cat("2. Usar genes_panel como features para clasificador\n\n")
cat("3. Entrenar modelos:\n")
cat("   - Random Forest\n")
cat("   - SVM\n")
cat("   - Logistic Regression con LASSO\n\n")
cat("4. Validar en cohorte independiente (GSE35958)\n\n")
cat("5. Análisis de feature importance\n\n")
cat("6. Reportar métricas finales (AUC, sensibilidad, especificidad)\n\n")

cat("═════════════════════════════════════════════════════════════\n\n")

# FIN DEL SCRIPT



# ============================================================================
# FIGURAS – FASE 2
# ============================================================================

dir.create("Fase2-Enrichment/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Fase2-Enrichment/tables", recursive = TRUE, showWarnings = FALSE)

# --- FIG 1: Overlap ---
library(ggvenn)
venn_data <- list(
  Edad = immu_age_sig$Gene,
  Replicativa = immu_7888_sig$Gene
)

p_overlap <- ggvenn(venn_data, fill_color = c("#3b82f6", "#ef4444"))
ggsave("Fase2-Enrichment/figures/Overlap_Inmune_Edad_vs_Rep.png",
       p_overlap, width = 5, height = 5, dpi = 300)

# --- FIG 2: Enriquecimiento GO ---
if (!is.null(ego_age) && nrow(ego_age@result) > 0) {
  p_go_age <- dotplot(ego_age, showCategory = 20) + ggtitle("GO BP – Edad")
  ggsave("Fase2-Enrichment/figures/GO_Edad.png",
         p_go_age, width = 8, height = 6, dpi = 300)
}

if (!is.null(ego_7888) && nrow(ego_7888@result) > 0) {
  p_go_rep <- dotplot(ego_7888, showCategory = 20) + ggtitle("GO BP – Replicativa")
  ggsave("Fase2-Enrichment/figures/GO_Replicativa.png",
         p_go_rep, width = 8, height = 6, dpi = 300)
}

# --- FIG 3: GSEA ---
if (!is.null(gsea_age_h) && nrow(gsea_age_h@result) > 0) {
  p_gsea_age <- gseaplot2(gsea_age_h, geneSetID = 1)
  ggsave("Fase2-Enrichment/figures/GSEA_Edad.png",
         p_gsea_age, width = 8, height = 6, dpi = 300)
}

if (!is.null(gsea_7888_h) && nrow(gsea_7888_h@result) > 0) {
  p_gsea_rep <- gseaplot2(gsea_7888_h, geneSetID = 1)
  ggsave("Fase2-Enrichment/figures/GSEA_Replicativa.png",
         p_gsea_rep, width = 8, height = 6, dpi = 300)
}

# --- FIG 4: Red regulatoria TF/miRNA ---
p_graph <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(color = Type), alpha = 0.4) +
  geom_node_point(aes(color = node_type, size = pagerank)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3)

ggsave("Fase2-Enrichment/figures/Red_Regulatoria.png",
       p_graph, width = 10, height = 8, dpi = 300)

# --- FIG 5: ROC curves ---
png("Fase2-Enrichment/figures/ROC_ImmunoScore.png",
    width = 1200, height = 900)
plot.roc(roc_score, col = "#10b981", lwd = 3, main = "ROC: ImmunoScore vs Edad")
lines.roc(roc_age_var, col = "#3b82f6", lwd = 2)
lines.roc(roc_pass, col = "#ef4444", lwd = 2)
lines.roc(roc_combined, col = "purple", lwd = 2)
legend("bottomright",
       legend = c("ImmunoScore", "Edad", "Pasaje", "Combinado"),
       col = c("#10b981", "#3b82f6", "#ef4444", "purple"),
       lwd = 3)
dev.off()

# --- FIG 6: ImmunoScore distribuciones ---
p_immuno <- ggplot(meta_39035, aes(GroupEdad, ImmunoScore, fill = GroupEdad)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() + ggtitle("ImmunoScore por Edad")

ggsave("Fase2-Enrichment/figures/ImmunoScore_Edad.png",
       p_immuno, width = 5, height = 6, dpi = 300)

write.csv(immu_age_sig,
          "Fase2-Enrichment/tables/Inmune_DEGs_Edad.csv",
          row.names = FALSE)

write.csv(immu_7888_sig,
          "Fase2-Enrichment/tables/Inmune_DEGs_Replicativa.csv",
          row.names = FALSE)

write.csv(top_regulators,
          "Fase2-Enrichment/tables/Top_Regulators_PageRank.csv",
          row.names = FALSE)

write.csv(rescue_results,
          "Fase2-Enrichment/tables/Rescue_Virtual_TFs.csv",
          row.names = FALSE)

write.csv(score_summary_age,
          "Fase2-Enrichment/tables/Scores_Edad.csv")

write.csv(score_summary_7888,
          "Fase2-Enrichment/tables/Scores_Replicativa.csv")

cat("✓ Tablas de Fase 2 guardadas en 'Fase2-Enrichment/tables'\n")
