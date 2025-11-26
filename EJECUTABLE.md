# 🧬 Predicción y Rescate de la Función Inmunomoduladora en hMSC mediante Análisis Transcriptómico y Machine Learning

---

# **ABSTRACT**

La eficacia terapéutica de las células madre mesenquimales humanas (hMSC) depende críticamente de su capacidad inmunomoduladora, la cual decae con la edad del donante y el número de pasajes in vitro durante su expansión. Actualmente no existe un panel transcriptómico estándar que permita cuantificar esta pérdida funcional ni regularizadores moleculares capaces de revertirla.

En este proyecto desarrollamos un pipeline integral de análisis transcriptómico (microarrays) dividido en **tres fases**:

1. **Fase 1:** Descubrimiento de firmas diferenciales de senescencia cronológica y replicativa.
2. **Fase 2:** Análisis de enriquecimiento, identificación de genes inmunomoduladores alterados, construcción de una red TF–miRNA–gene y desarrollo de un panel preliminar **MSC-ImmunoScore**.
3. **Fase 3:** Entrenamiento de un modelo de Machine Learning capaz de clasificar hMSC funcionales vs. senescentes, validado en un dataset externo.

Los resultados muestran >170 genes inmunes alterados en senescencia, pérdida consistente de PD-L1/IDO1, activación de SASP, y un panel transcriptómico con fuerte capacidad predictiva. El modelo final obtuvo **92.7 % de accuracy (AUC = 0.978)** en validación cruzada y **78 % de accuracy** en un dataset externo. Este trabajo aporta un pipeline reproducible, un panel preliminar y reguladores candidatos para orientar estrategias de rescate funcional en hMSC para aplicaciones clínicas.

---

# 📌 1. Descripción General

Este repositorio contiene el pipeline completo para:

* Procesar e integrar datasets transcriptómicos de hMSC
* Identificar firmas de senescencia
* Detectar pérdida de inmunomodulación
* Priorizar genes reguladores
* Construir un modelo predictivo basado en Machine Learning
* Validarlo en datos externos
---

# 🎯 2. Objetivos del Proyecto

### **Objetivo general**

Predecir la pérdida de inmunomodulación en hMSC mediante análisis transcriptómico y Machine Learning, e identificar reguladores con potencial de rescate funcional.

### **Objetivos específicos**

* Definir una **firma core de senescencia** a partir de múltiples datasets.
* Identificar **genes inmunomoduladores alterados** en senescencia.
* Analizar vías y rutas clave mediante GO/KEGG/GSEA.
* Construir una red regulatoria TF–miRNA–gene.
* Derivar un **MSC-ImmunoScore** basado en expresión génica.
* Entrenar y validar un **clasificador funcional** en hMSC.

---

# 🧩 3. Estructura del Repositorio

```
bioinformatica-avanzada-hmsc-prediccion/

├── README.md                     # Documentación principal
│
│
├── Fase0           
│
├
├── Fase1-Discovery/              # Preprocesamiento y firma de senescencia
│   ├── Fase1_script.R
│   ├── results/
│   │   ├── FASE1_DEGs_export.rds
│   │   ├── FASE1_expr_matrices.rds
│   │   └── CORE_senescence_signature.csv
│   ├── figures/
│   └── tables/
│
├── Fase2-Enrichment/             # Inmunomodulación + Red regulatoria
│   ├── Fase2_script.R
│   ├── results/
│   │   ├── FASE2_complete_export.rds
│   │   ├── MSC_ImmunoScore_panel.csv
│   │   ├── top_regulators_pagerank.csv
│   │   └── rescue_analysis_results.csv
│   ├── figures/
│   └── tables/
│
└── Fase3 - Modelo de clasificación Machine Learning/
    ├── train_model.py
    ├── validate_model.py
    ├── X_features_matrix.csv
    ├── y_labels.csv
    ├── metadata_samples.csv
    ├── best_model_msc_senescence.pkl
    ├── scaler_msc_senescence.pkl
    ├── feature_importance.png
    ├── .gitignore
    ├── ControlNegativo/
    ├── ControlPositivo/
    └── validation/
        ├── X_test_GSE35958.csv
        ├── y_test_GSE35958.csv
        └── metadata_GSE35958.csv
```

---

# 🧪 4. Datasets Utilizados

| Dataset      | Descripción              | Uso                     |
| ------------ | ------------------------ | ----------------------- |
| **GSE39035** | hMSC jóvenes vs ancianas | Senescencia por edad    |
| **GSE7888**  | Pasajes early–mid–late   | Senescencia replicativa |
| **GSE35958** | Donantes ancianos        | Validación externa ML   |

Todos son microarrays Affymetrix Human Genome U133 series.

---

# ⚙️ 5. Requisitos e Instalación

### **R ≥ 4.2**

Paquetes principales:

```
limma
GEOquery
tidyverse
clusterProfiler
ComplexHeatmap
msigdbr
enrichplot
igraph
ggraph
pROC
patchwork
```

### **Python ≥ 3.10**

```
pip install -r requirements.txt
```

Incluye:

* scikit-learn
* numpy
* pandas
* matplotlib
* joblib

---

# ▶️ 6. Cómo Ejecutar el Proyecto

## **FASE 1 — Descubrimiento (R)**

Procesamiento, normalización, DEGs y firma core.

```r
source("Fase1-Discovery/Fase1_script.R")
```

Genera:

* FASE1_DEGs_export.rds
* FASE1_expr_matrices.rds
* Figuras: PCA, volcano, heatmaps
* Tablas con DEGs por contraste

---

## **FASE 2 — Inmunomodulación y Red Reguladora (R)**

Filtrado inmune, GO/KEGG/GSEA, TF/miRNA, MSC-ImmunoScore.

```r
source("Fase2-Enrichment/Fase2_script.R")
```

Genera:

* FASE2_complete_export.rds
* Panel MSC-ImmunoScore
* Red regulatoria TF–miRNA–gene
* ROC curves
* Enriquecimiento GO/KEGG/GSEA
* Análisis de rescate virtual

---

## **FASE 3 — Machine Learning (Python)**

Entrena un modelo para predecir funcionalidad inmunomoduladora.

### Entrenar modelo

```bash
cd "Fase 3 - Modelo de clasificación Machine Learning"
python train_model.py
```

### Validar en dataset externo

```bash
python validate_model.py
```

---

# 📊 7. Resultados Principales

### **Senescencia**

* Pasaje tiene mayor impacto transcriptómico que edad.
* Firma core compartida entre datasets (SOX11, EMX2OS, DDIT4L).

### **Inmunomodulación**

* > 170 genes inmunes alterados en senescencia replicativa.
* PD-L1 e IDO1 consistentemente downregulated.
* Enriquecimiento de rutas inflamatorias (SASP, NF-κB, complement).

### **Reguladores candidatos**

* Módulos TF/miRNA con regulación distribuida.
* Candidatos para rescate: IRF1, STAT1/3, miR-146a, miR-155.

### **Panel MSC-ImmunoScore**

* Panel preliminar de 33 genes.
* Diferencia clara entre hMSC funcionales vs senescentes.

### **Machine Learning**

* **Accuracy CV:** 92.7 %
* **AUC:** 0.978
* **Validación externa (GSE35958):** 78 % accuracy

---

# 🖼️ 8. Figuras Principales (Thumbnails)
```
Fase1-Discovery/figures/PCA_GSE39035_GSE7888.png
Fase1-Discovery/figures/Volcano_Edad.png
Fase2-Enrichment/figures/GO_Edad.png
Fase2-Enrichment/figures/Red_Regulatoria.png
Fase3 - Modelo de clasificación Machine Learning/feature_importance.png
```

