############################################################
# VALIDACIÓN EXTERNA — GSE35958 (Osteoporosis vs Control)
############################################################

library(GEOquery)
library(limma)
library(dplyr)
library(tibble)
library(stringr)
library(matrixStats)

#############################
# 1. Descargar dataset
#############################

gse <- getGEO("GSE35958", GSEMatrix = TRUE)[[1]]

expr_raw  <- exprs(gse)
meta_raw  <- pData(gse)

#############################
# 2. Descargar anotación GPL570
#############################

gpl <- getGEO("GPL570", AnnotGPL = TRUE)
annot <- Table(gpl)

# detectar columna correcta
symbol_col_candidates <- c(
  "Gene symbol", "gene_symbol", "GeneSymbol", "GENE_SYMBOL", "Symbol"
)

symbol_col <- intersect(symbol_col_candidates, colnames(annot))
symbol_col <- symbol_col[1]

cat("Usando columna de símbolo génico:", symbol_col, "\n")

annot_clean <- annot %>%
  dplyr::select(ID, !!sym(symbol_col)) %>%
  dplyr::rename(GENE_SYMBOL = !!sym(symbol_col)) %>%
  dplyr::filter(GENE_SYMBOL != "", GENE_SYMBOL != "---") %>%
  dplyr::mutate(ID = as.character(ID))

#############################
# 3. Unir expresión con anotación
#############################

expr_df <- expr_raw %>%
  as.data.frame() %>%
  rownames_to_column("ID")

expr_annot <- inner_join(expr_df, annot_clean, by = "ID")

# mover símbolo al principio
expr_annot <- expr_annot %>%
  relocate(GENE_SYMBOL)

# eliminar columna ID
expr_annot <- expr_annot %>%
  dplyr::select(-ID)

# Checar si todas las demás son NUMÉRICAS
non_num <- sapply(expr_annot[,-1], function(x) !is.numeric(x))
if (any(non_num)) {
  stop("Hay columnas no-numéricas en la expresión después de merge.")
}

#############################
# 4. Colapsar probes → gen
#############################

expr_gene <- expr_annot %>%
  group_by(GENE_SYMBOL) %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  as.data.frame()

rownames(expr_gene) <- expr_gene$GENE_SYMBOL
expr_gene <- expr_gene[, -1]

#############################
# FIX PARA COLUMNAS "list"
#############################

# detectar columnas tipo list
list_cols <- sapply(expr_gene, is.list)

if (any(list_cols)) {
  cat("Columnas list detectadas y convertidas a numeric:\n")
  print(names(expr_gene)[list_cols])
  
  expr_gene[list_cols] <- lapply(expr_gene[list_cols], function(col) {
    suppressWarnings(as.numeric(unlist(col)))
  })
}

# convertir todo lo demás a numeric
expr_gene <- expr_gene %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))

# verificar que ya todo es numérico
stopifnot(all(sapply(expr_gene, is.numeric)))


#############################
# 5. Detectar escala
#############################

q95 <- quantile(as.numeric(as.matrix(expr_gene)), 0.95)

if (q95 > 20) {
  cat("Dataset está en escala lineal. Aplicando log2(x+1)...\n")
  expr_gene <- log2(expr_gene + 1)
} else {
  cat("Dataset ya está en escala log2.\n")
}

#############################
# 6. Normalizar
#############################

expr_norm <- normalizeBetweenArrays(expr_gene)

#############################
# 7. Filtrar baja expresión
#############################

expr_filt <- expr_norm
cat("Genes finales:", nrow(expr_filt), "\n")

#############################
# 8. Metadata y grupos
#############################

meta <- meta_raw %>%
  dplyr::mutate(
    Sample = geo_accession,
    Group = ifelse(
      grepl("osteopor|OP", title, ignore.case = TRUE),
      "OP", "Control"
    ),
    Label = ifelse(Group == "OP", 1L, 0L)
  )

#############################
# 9. Alinear con matriz
#############################

common <- intersect(colnames(expr_filt), meta$Sample)

expr_filt <- expr_filt[, common]
meta <- meta %>%
  dplyr::filter(Sample %in% common) %>%
  dplyr::arrange(match(Sample, colnames(expr_filt)))

stopifnot(all(colnames(expr_filt) == meta$Sample))

#############################
# 10. Exportar
#############################

dir.create("ML_validation_GSE35958", showWarnings = FALSE)

write.csv(
  expr_filt %>% as.data.frame() %>% rownames_to_column("Gene"),
  "ML_validation_GSE35958/X_test_GSE35958.csv",
  row.names = FALSE
)

write.csv(
  meta %>% dplyr::select(Sample, Group, Label, title),
  "ML_validation_GSE35958/y_test_GSE35958.csv",
  row.names = FALSE
)

write.csv(
  meta,
  "ML_validation_GSE35958/metadata_GSE35958.csv",
  row.names = FALSE
)

cat("\n✔ GSE35958 listo para ML (control positivo)\n")





###############################################################
###############################################################
###############################################################



###############################################################
# Construcción de X_train y y_train a partir de:
# - GSE39035 (Edad: Young vs Old)
# - GSE7888  (Pasaje: Early/Mid vs Late)
#
# SIN INTERSECCIÓN CON PANEL MSC-ImmunoScore
###############################################################

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)

data_dir <- "Fase2_inputs"   # <-- Ajusta la carpeta si es necesario

###############################################################
# 1. CARGAR EXPRESIÓN NORMALIZADA + METADATA
###############################################################

expr_39035 <- read.csv(
  file.path(data_dir, "GSE39035_expression_normalized_filtered.csv"),
  row.names = 1, check.names = FALSE
)

meta_39035 <- read.csv(
  file.path(data_dir, "GSE39035_metadata.csv"),
  check.names = FALSE
)

expr_7888 <- read.csv(
  file.path(data_dir, "GSE7888_expression_normalized_filtered.csv"),
  row.names = 1, check.names = FALSE
)

meta_7888 <- read.csv(
  file.path(data_dir, "GSE7888_metadata.csv"),
  check.names = FALSE
)


###############################################################
# 2. DEFINIR LABELS (tu definición final)
###############################################################

### ---- GSE39035 ----
# GroupEdad debe tener: Young, Old

meta_39035 <- meta_39035 %>%
  mutate(
    Sample = Sample,
    Label = case_when(
      GroupEdad %in% c("Young") ~ 0,
      GroupEdad %in% c("Old")   ~ 1,
      TRUE ~ NA_integer_
    )
  )

table(meta_39035$GroupEdad, meta_39035$Label)


### ---- GSE7888 ----
# GroupPassage debe tener: Early / Mid / Late

meta_7888 <- meta_7888 %>%
  mutate(
    Sample = Sample,
    Label = case_when(
      GroupPassage %in% c("Early")       ~ 0,
      GroupPassage %in% c("Late")        ~ 1,
      TRUE ~ NA_integer_
    )
  )

table(meta_7888$GroupPassage, meta_7888$Label)


###############################################################
# 3. COMBINAR EXPRESIÓN (genes en filas → muestras en columnas)
###############################################################

# Intersección de genes expresados en ambos datasets
common_genes <- intersect(rownames(expr_39035), rownames(expr_7888))
length(common_genes)

expr_39035_f <- expr_39035[common_genes, ]
expr_7888_f  <- expr_7888[common_genes, ]

# Construir matrices (muestras x genes)
X_39035 <- t(expr_39035_f)
X_7888  <- t(expr_7888_f)

# Unir
X_train <- rbind(X_39035, X_7888)
X_train <- as.data.frame(X_train)

###############################################################
# 4. COMBINAR METADATA + LABELS (y_train)
###############################################################

meta_39035$Dataset <- "GSE39035"
meta_7888$Dataset  <- "GSE7888"

meta_all <- bind_rows(meta_39035, meta_7888)

# Reordenar metadata según X_train
meta_all <- meta_all[match(rownames(X_train), meta_all$Sample), ]

stopifnot(all(rownames(X_train) == meta_all$Sample))

y_train <- meta_all %>%
  select(Sample, Dataset, Label, everything())

table(y_train$Label)

###############################################################
# 5. EXPORTAR PARA PYTHON / ML
###############################################################

out_dir <- "ML_training_data"
if (!dir.exists(out_dir)) dir.create(out_dir)

write.csv(X_train,
          file.path(out_dir, "X_train_fullgenes.csv"))

write.csv(y_train,
          file.path(out_dir, "y_train_labels.csv"),
          row.names = FALSE)

write.table(
  rownames(expr_39035),
  file.path(out_dir, "genes_full_list.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

cat("✔ X_train_fullgenes.csv creado\n")
cat("✔ y_train_labels.csv creado\n")
cat("✔ genes_full_list.txt creado\n")
