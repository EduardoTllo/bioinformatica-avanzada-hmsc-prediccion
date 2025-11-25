# Bioinformática Avanzada: Predicción y Rescate de la Función Inmunomoduladora en hMSC mediante Análisis Transcriptómico

## 1. Temática General del Proyecto

El proyecto busca predecir y entender la pérdida de la función inmunomoduladora en células madre mesenquimales humanas (hMSC) a partir de datos transcriptómicos. La idea central es identificar biomarcadores de senescencia (por edad del donante y por número de pasajes in vitro) que expliquen por qué algunas hMSC dejan de suprimir respuestas inmunes y, a partir de eso, proponer reguladores potenciales para “rescatar” esa función.

## 2. Introducción y Contexto

Las hMSC se usan en terapia celular por tres propiedades clave: **inmunomodulación, autorrenovación y potencial de diferenciación**.

En bioprocesos industriales el gran reto es la **escalabilidad**: hay que expandir las hMSC en cultivo para obtener dosis terapéuticas suficientes. Sin embargo, hoy no existen protocolos estandarizados sobre la edad máxima del donante ni el número máximo de pasajes aceptables, lo que lleva a una alta variabilidad en la calidad inmunomoduladora entre lotes de hMSC.

El proyecto se plantea entonces como un esfuerzo por definir **criterios transcriptómicos objetivos** que permitan anticipar cuándo una hMSC ha perdido su “buena” inmunomodulación.

## 3. Estado del Arte y Vacío Identificado

### Qué se sabe
La senescencia cronológica (donantes ancianos) y la senescencia replicativa (muchos pasajes) se asocian a una **disminución de 30–50 % en la capacidad de suprimir linfocitos T** respecto a hMSC de donantes jóvenes.

Esta pérdida funcional se vincula con:
*   **Aumento de SASP** (IL-6, IL-8, quimiocinas) y daño en el ADN.
*   **Activación crónica de NF-κB y rutas MAPK**.
*   Cambios metabólicos y en fosfolípidos que modulan la respuesta inmune.

A nivel molecular, se ha visto que en hMSC envejecidas disminuyen **PD-L1 e IDO1**, dos reguladores clave de inmunosupresión, en parte por la reducción del factor de transcripción **GATA2**. La sobreexpresión de GATA2 puede rescatar parcialmente la inmunomodulación.

### Qué falta
Existen estudios de transcriptómica de senescencia, pipelines de microarrays y hasta modelos de ML tipo SenPred, pero **no hay un panel transcriptómico estándar** específicamente diseñado para predecir la pérdida de inmunomodulación en hMSC, ni criterios claros de edad/pasaje aceptables para uso terapéutico.

### Pregunta Central
> **¿Qué biomarcadores transcriptómicos predicen la pérdida de inmunomodulación en hMSC y qué reguladores podrían rescatarla?**

---

## 4. Fase 1 – Descubrimiento de Firma “Core” de Senescencia

### Objetivo Específico
Identificar firmas transcriptómicas diferenciales asociadas a edad y pasaje en hMSC y derivar una **firma core de senescencia** compartida por ambas formas (cronológica y replicativa).

### Datos y Diseño
Se integran dos datasets de microarrays:
1.  **GSE39035**: Diseño factorial (donantes jóvenes vs ancianos, distintos pasajes).
2.  **GSE7888**: Serie de pasajes (early, mid, late).

### Pipeline Bioinformático
1.  Carga de datos y metadata (`GEOQuery`).
2.  Anotación de probes a genes y colapso por mediana.
3.  Transformación log2 + normalización cuantílica para unificar plataformas.
4.  Filtrado de baja expresión.
5.  Corrección de efectos de batch/donante usando `removeBatchEffect` y `duplicateCorrelation`.
6.  Modelado lineal con `limma` con varios contrastes:
    *   Edad: Old vs Young
    *   Pasaje: High vs Low, Late vs Early
    *   Efecto aditivo Edad+Pasaje
    *   Interacción Edad×Pasaje

### Resultados Clave
*   Tras la corrección de batch, la PCA muestra que emergen claramente los patrones biológicos de Edad y Pasaje.
*   El número de DEGs es mayor para Pasaje que para Edad, indicando que la senescencia replicativa tiene un impacto transcriptómico más profundo.
*   La interacción Edad×Pasaje ≈ 0, sugiriendo efectos mayormente aditivos.
*   **Firma Core**: Se identificaron genes como **EMX2OS, SOX11 y DDIT4L** (downregulated en ambos contextos), indicando una pérdida funcional estable.

---

## 5. Fase 2 – Senescencia e Inmunomodulación

### Objetivo Específico
Identificar genes y vías inmunomoduladoras alteradas durante la senescencia de hMSC y proponer reguladores upstream (TF/miRNA) y un panel preliminar **MSC-ImmunoScore**.

### Enfoque General
1.  Construir un “catálogo inmune” usando anotaciones GO y MSigDB.
2.  Filtrar los DEGs de Fase 1 a DEGs inmunomoduladores.
3.  Realizar enriquecimiento funcional GO/KEGG.
4.  Explorar reguladores upstream (TF y miRNAs).
5.  Integrar resultados en un panel preliminar.

### Resultados Clave
*   **Catálogo Inmune**: En Pasaje High vs Low, >170 genes inmunes son diferenciales (muchos downregulated), reflejando un apagamiento de funciones inmunes.
*   **Enriquecimiento Funcional**:
    *   **GO BP**: Regulación de linfocitos/leucocitos, respuesta inflamatoria, producción de citoquinas.
    *   **KEGG**: Complemento y cascadas de coagulación.
*   **Reguladores Upstream**: Regulación distribuida (módulos miRNA-target) más que un único "master regulator".
*   **Panel MSC-ImmunoScore**: Genes candidatos priorizados (e.g., **SOX11, EMX2OS, RBP4, NTF3, ND1N, DPPA3, RRAGD, BST1, TNFRSF11B**) que capturan señal de senescencia e impacto inmune.

---

## 6. Fase 3 – Modelo de Clasificación Machine Learning

Esta fase utiliza el panel de genes identificado para entrenar un modelo predictivo capaz de clasificar nuevas muestras.

### 📊 Datasets

#### Datos de Entrenamiento
- **Samples**: 28 muestras de hMSC (13 funcionales, 15 senescentes).
- **Features**: 33 genes seleccionados del panel MSC-ImmunoScore.
- **Source**: Datos de expresión génica normalizados (log2).

#### Datos de Validación Externa
- **Dataset**: GSE35958
- **Samples**: 9 muestras de donantes ancianos (79-94 años).
- **Grupos**: Controles ancianos y pacientes con osteoporosis.

### 🧬 Genes Biomarcadores Identificados (Top Features)

| Rank | Gene | Importancia | Función Biológica |
|------|------|-------------|-------------------|
| 1 | **SCN9A** | 0.1853 | Canal de sodio, asociado con senescencia |
| 2 | **HDAC9** | 0.0979 | Histona deacetilasa, regulación epigenética |
| 3 | **KCTD16** | 0.0897 | Regulación de la degradación proteica |
| 4 | **CD55** | 0.0627 | Proteína reguladora del complemento |
| 5 | **EPHA5** | 0.0545 | Receptor tirosina quinasa |

### 🤖 Resultados del Modelo (Random Forest)

#### Entrenamiento (5-Fold Cross-Validation)
| Métrica | Valor |
|---------|-------|
| **Accuracy** | 92.7% |
| **ROC-AUC** | **0.978** |

#### Validación Externa (GSE35958)
| Métrica | Valor |
|---------|-------|
| **Accuracy** | 78% |
| **Recall** | 78% |
| **Precision** | 100% |

**Interpretación**: El modelo identifica correctamente el 78% de las muestras senescentes en un dataset independiente de donantes ancianos.

---

## 📁 Estructura del Proyecto

```
bioinformatica-avanzada-hmsc-prediccion/
├── README.md                          # Documentación principal
└── Fase 3 - Modelo de clasificación Machine Learning/
    ├── train_model.py                 # Script principal de entrenamiento
    ├── validate_model.py              # Script de validación
    ├── X_features_matrix.csv          # Matriz de características (genes)
    ├── y_labels.csv                   # Etiquetas de clase
    ├── metadata_samples.csv           # Metadata de las muestras
    ├── best_model_msc_senescence.pkl  # Modelo entrenado
    ├── scaler_msc_senescence.pkl      # Escalador StandardScaler
    ├── feature_importance.png         # Gráfico de importancia de features
    ├── .gitignore                     # Archivos ignorados por Git
    ├── ControlNegativo/               # Datos de control negativo
    ├── ControlPositivo/               # Datos de control positivo
    └── validation/                    # Datos de validación GSE35958
        ├── X_test_GSE35958 (1).csv
        ├── y_test_GSE35958 (1) copy.csv
        └── metadata_GSE35958 (1).csv
```

## 🚀 Uso (Fase 3)

### 1. Entrenamiento del Modelo
```bash
cd "Fase 3 - Modelo de clasificación Machine Learning"
python train_model.py
```

### 2. Validación con Datos Externos
```bash
cd "Fase 3 - Modelo de clasificación Machine Learning"
python validate_model.py
```

## 👤 Autor
**Eduardo**
- Proyecto: Bioinformática Avanzada
- Año: 2025

## 📄 Licencia
Este proyecto está bajo la Licencia MIT.
