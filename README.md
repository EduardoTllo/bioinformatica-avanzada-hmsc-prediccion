# Bioinformática Avanzada: Predicción y Rescate de la Función Inmunomoduladora en hMSC mediante Análisis Transcriptómico

## 📋 Descripción del Proyecto

Este proyecto utiliza técnicas de Machine Learning para predecir la senescencia en células madre mesenquimales humanas (hMSC) basándose en perfiles de expresión génica. El objetivo es identificar biomarcadores transcriptómicos que permitan distinguir entre hMSC funcionales y senescentes, facilitando estrategias de rescate de la función inmunomoduladora.

## 🎯 Objetivos

- **Predicción de Senescencia**: Desarrollar modelos de clasificación para identificar hMSC senescentes vs. funcionales
- **Identificación de Biomarcadores**: Descubrir genes clave asociados con la senescencia en hMSC
- **Validación Externa**: Evaluar el rendimiento del modelo en datasets independientes (GSE35958)

## 📊 Datasets

### Datos de Entrenamiento
- **Samples**: 28 muestras de hMSC
  - 13 muestras funcionales (jóvenes)
  - 15 muestras senescentes (ancianas)
- **Features**: 33 genes seleccionados
- **Source**: Datos de expresión génica normalizados (log2)

### Datos de Validación
- **Dataset**: GSE35958
- **Samples**: 9 muestras de donantes ancianos (79-94 años)
- **Grupos**: Controles ancianos y pacientes con osteoporosis

## 🧬 Genes Biomarcadores Identificados

### Top 10 Features (por importancia)

| Rank | Gene | Importancia | Función Biológica |
|------|------|-------------|-------------------|
| 1 | **SCN9A** | 0.1853 | Canal de sodio, asociado con senescencia |
| 2 | **HDAC9** | 0.0979 | Histona deacetilasa, regulación epigenética |
| 3 | **KCTD16** | 0.0897 | Regulación de la degradación proteica |
| 4 | **CD55** | 0.0627 | Proteína reguladora del complemento |
| 5 | **EPHA5** | 0.0545 | Receptor tirosina quinasa |
| 6 | **SOX11** | 0.0495 | Factor de transcripción |
| 7 | **C8orf34** | 0.0479 | Función desconocida |
| 8 | **FGD4** | 0.0478 | Activador de GTPasa |
| 9 | **EMX2** | 0.0438 | Desarrollo y diferenciación |
| 10 | **RBM24** | 0.0361 | Proteína de unión a ARN |

## 🤖 Modelos Implementados

### Algoritmos Evaluados
1. **Random Forest** ⭐ (Mejor modelo)
2. Logistic Regression
3. Support Vector Machine (Linear, Poly, RBF)
4. Decision Tree

### Resultados del Modelo Final (Random Forest)

#### Entrenamiento (5-Fold Cross-Validation)
| Métrica | Valor |
|---------|-------|
| **Accuracy** | 92.7% |
| **Precision** | 93.3% |
| **Recall** | 93.3% |
| **F1-Score** | 92.0% |
| **ROC-AUC** | **0.978** |

#### Validación Externa (GSE35958)
| Métrica | Valor |
|---------|-------|
| **Accuracy** | 78% |
| **Recall** | 78% |
| **Precision** | 100% |
| True Positives | 7/9 |
| False Negatives | 2/9 |

### Matriz de Confusión (Entrenamiento)

```
                Predicho Senescente (0)  Predicho Funcional (1)
Actual Senescente (0)        13                    1
Actual Funcional (1)          1                   13
```

**Interpretación**: Solo 2 muestras mal clasificadas de 28 totales.

## 🛠️ Tecnologías y Herramientas

- **Lenguaje**: Python 3.13
- **Machine Learning**: scikit-learn
- **Análisis de Datos**: pandas, numpy
- **Visualización**: matplotlib, seaborn
- **Validación**: Stratified K-Fold Cross-Validation

## 📁 Estructura del Proyecto

```
bioinformatica-avanzada-hmsc-prediccion/
├── README.md                          # Documentación principal
└── Modelo Machine Learning/
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

## 🚀 Uso

### 1. Entrenamiento del Modelo

```bash
cd "Modelo Machine Learning"
python train_model.py
```

Este script:
- Carga y preprocesa los datos
- Entrena múltiples modelos (Random Forest, SVM, Logistic Regression, Decision Tree)
- Evalúa con 5-Fold Cross-Validation
- Guarda el mejor modelo
- Genera visualizaciones (matriz de confusión, feature importance)

### 2. Validación con Datos Externos

```bash
cd "Modelo Machine Learning"
python validate_model.py
```

Este script:
- Carga el modelo entrenado
- Procesa datos de validación (GSE35958)
- Alinea features faltantes
- Calcula métricas de rendimiento

### 3. Uso del Modelo para Predicción

```python
import joblib
import pandas as pd
import os

# Cambiar al directorio del modelo
os.chdir("Modelo Machine Learning")

# Cargar modelo y scaler
model = joblib.load('best_model_msc_senescence.pkl')
scaler = joblib.load('scaler_msc_senescence.pkl')

# Cargar nuevos datos
X_new = pd.read_csv('nuevas_muestras.csv')

# Preprocesar
X_scaled = scaler.transform(X_new)

# Predecir
predictions = model.predict(X_scaled)
probabilities = model.predict_proba(X_scaled)

print(f"Predicción: {predictions}")
print(f"Probabilidades: {probabilities}")
```

## 📈 Análisis de Resultados

### Hallazgos Clave

1. **Alto Rendimiento en Entrenamiento**: El modelo Random Forest alcanza 97.8% de ROC-AUC en validación cruzada.

2. **Validación Externa Robusta**: 78% de accuracy en muestras de donantes ancianos (GSE35958), demostrando buena generalización.

3. **SCN9A como Biomarcador Principal**: El gen **SCN9A** (canal de sodio) es el predictor más importante (18.5% de importancia), sugiriendo un rol crucial en la senescencia de hMSC.

4. **Perfil de Error**: El modelo tiene alta precisión (100%) cuando predice senescencia, pero puede tener falsos negativos (2/9 muestras), indicando que algunos donantes ancianos conservan perfiles de expresión "jóvenes".

## 🔬 Implicaciones Biológicas

- **SCN9A** y **HDAC9** emergen como potenciales dianas terapéuticas para rescatar la función inmunomoduladora
- La expresión génica puede ser un mejor indicador de senescencia funcional que la edad cronológica
- Algunos donantes ancianos mantienen perfiles transcriptómicos juveniles, sugiriendo heterogeneidad en el envejecimiento

## 📚 Referencias

- Dataset de Validación: GSE35958 (Osteoporosis vs Control, donantes ancianos 79-94 años)
- Metodología: Random Forest con Stratified 5-Fold Cross-Validation

## 👤 Autor

**Eduardo**
- Proyecto: Bioinformática Avanzada
- Institución: [Tu Institución]
- Año: 2025

## 📄 Licencia

Este proyecto está bajo la Licencia MIT - ver el archivo LICENSE para más detalles.

## 🙏 Agradecimientos

- GEO (Gene Expression Omnibus) por los datos públicos de validación
- Comunidad de scikit-learn por las herramientas de ML

---

**Nota**: Este proyecto es parte de un estudio de investigación en bioinformática aplicada a células madre mesenquimales y senescencia celular.
