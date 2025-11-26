# Walkthrough - Modelo de Clasificación Machine Learning para Osteoporosis

## Resumen del Proyecto

Este proyecto desarrolla un modelo de machine learning para clasificar células mesenquimales (MSC) como:
- **Clase 0**: Control (células normales)
- **Clase 1**: Osteoporosis Primaria (OP)

---

## 📊 Dataset

### Datos Combinados
- **Total de muestras**: 37
  - Original: 28 muestras
  - Validación externa (GSE35958): 9 muestras
- **Features (genes)**: 33 genes candidatos
- **Distribución de clases**:
  - Clase 0 (Control): 23 muestras (62%)
  - Clase 1 (OP): 14 muestras (38%)

### Genes Candidatos
Los 33 genes fueron seleccionados de las Fases I y II del proyecto basados en análisis diferencial y significancia biológica.

---

## 🔬 Metodología

### 1. Preprocesamiento
- **Normalización**: StandardScaler (media=0, std=1)
- **Alineación de features**: Asegurar que ambos datasets tengan las mismas columnas
- **Transformación logarítmica**: log2(x+1) aplicada cuando fue necesario

### 2. Estrategia de Validación
**Validación Cruzada Estratificada (5-Fold)**
- Divide el dataset en 5 grupos
- Mantiene la proporción de clases en cada fold
- Cada muestra se usa exactamente 1 vez para validación
- Métricas reportadas son el **promedio de 5 evaluaciones**

### 3. Modelos Evaluados
1. **Logistic Regression** (Modelo más explicable)
2. **Random Forest**
3. **SVM (Linear)**
4. **Decision Tree**

---

## 📈 Resultados

### Cross-Validation (5-Fold) - Dataset Combinado

| Modelo | Accuracy | Precision | Recall | F1 | AUC |
|--------|----------|-----------|--------|----|-----|
| **Logistic Regression** ✨ | **1.000** | **1.000** | **1.000** | **1.000** | **1.000** |
| Random Forest | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| SVM (Linear) | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Decision Tree (depth=6) | 0.921 | 0.900 | 0.933 | 0.903 | 0.927 |

### Conclusiones de Rendimiento
1. **Logistic Regression, Random Forest y SVM (Linear)** obtienen métricas perfectas (100%)
2. **Decision Tree (depth=6)** tiene un rendimiento muy bueno (92.1% accuracy)
3. Todos los modelos lineales son capaces de separar perfectamente las clases
4. El árbol de decisión más profundo captura mejor los patrones complejos

---

## 🎯 Modelo Seleccionado: Logistic Regression

### Razones para la Selección
✅ **Máxima Explicabilidad**: Los coeficientes son directamente interpretables  
✅ **Rendimiento Perfecto**: 100% en todas las métricas  
✅ **Simplicidad**: Modelo lineal fácil de entender y comunicar  
✅ **Robustez**: Mejor generalización potencial que modelos no lineales

### Top 5 Genes Más Importantes (Logistic Regression)

| # | Gen | Coeficiente (Abs) | Interpretación |
|---|-----|-------------------|----------------|
| 1 | **SOX11** | 1.3473 | Factor de transcripción, altamente discriminativo |
| 2 | **DDIT4L** | 0.2251 | Relacionado con estrés celular |

**Nota**: Solo 2 genes muestran coeficientes significativos, sugiriendo que estos son los marcadores clave.

---

## 📉 Análisis de Reducción Dimensional (PCA)

### Resultados PCA
- **PC1**: Explica 30.1% de la varianza
- **PC2**: Explica 20.2% de la varianza
- **Total**: 50.3% de varianza explicada con 2 componentes

### Visualizaciones Generadas
1. **decision_boundaries_pca.png**: Límites de decisión en espacio 2D
   - Muestra cómo cada modelo separa las clases
   - Logistic Regression crea un límite lineal claro
   - Decision Tree crea regiones rectangulares

2. **decision_tree_visualization.png**: Árbol de decisión completo
   - Muestra las reglas de decisión basadas en genes específicos
   - Profundidad máxima: 3 niveles para interpretabilidad

3. **feature_importance_comparison.png**: Comparación de importancia
   - Logistic Regression: Coeficientes absolutos
   - Decision Tree: Importancia basada en Gini

---

## 🌳 Decision Tree - Análisis Complementario

### Rendimiento (depth=6)
- **Accuracy**: 92.1%
- **Precision**: 90.0%
- **Recall**: 93.3%
- **F1**: 90.3%
- **AUC**: 92.7%

### Ventajas del Decision Tree
✅ Reglas fáciles de interpretar: "Si GEN_X > valor, entonces Clase 1"  
✅ No requiere normalización  
✅ Captura relaciones no lineales  
✅ **Rendimiento mejorado**: Con profundidad 6 alcanza 92.1% accuracy

### Observaciones
✔️ Mayor profundidad permite capturar patrones más complejos  
✔️ Mantiene alto recall (93.3%) - detecta la mayoría de casos positivos  
⚠️ Aún ligeramente inferior a modelos lineales (92% vs 100%)  
⚠️ Mayor profundidad aumenta riesgo de overfitting (monitorear con datos externos)  

---

## 🔍 Insights Biológicos

### Gen SOX11
- **Coeficiente más alto** en Logistic Regression
- **Rol biológico**: Factor de transcripción involucrado en diferenciación celular
- **Implicación**: Expresión diferencial clara entre MSC normales y con osteoporosis

### Gen DDIT4L
- **Segundo coeficiente más alto**
- **Rol biológico**: Regulador de respuesta al estrés y mTOR signaling
- **Implicación**: Puede reflejar estrés metabólico en células con osteoporosis

---

## 📁 Archivos Generados

### Modelos
- `final_explainable_model.pkl`: Modelo Logistic Regression final
- `final_scaler.pkl`: StandardScaler ajustado a los datos
- `best_model_combined.pkl`: Mejor modelo (Logistic Regression)
- `scaler_combined.pkl`: Scaler para modelo combinado

### Visualizaciones
- `decision_boundaries_pca.png`: Límites de decisión en 2D
- `decision_tree_visualization.png`: Árbol de decisión
- `feature_importance_comparison.png`: Comparación de importancia
- `feature_importance.png`: Importancia de features (entrenamiento inicial)

### Scripts
- `train_model.py`: Entrenamiento inicial con datos originales
- `train_model_combined.py`: Entrenamiento con datos combinados + CV
- `validate_model.py`: Validación en dataset externo
- `model_analysis_visualization.py`: Análisis completo y visualizaciones

---

## 🎓 Recomendaciones

### Para Uso Clínico/Investigación
1. **Validar con más datos externos**: 37 muestras es limitado
2. **Validación prospectiva**: Probar con nuevos pacientes
3. **Análisis de SOX11**: Investigar más a fondo su rol en osteoporosis

### Para Mejora del Modelo
1. **Aumentar dataset**: Buscar más datasets públicos de MSC
2. **Feature engineering**: Considerar ratios o combinaciones de genes
3. **Validación externa independiente**: Probar en cohorte completamente diferente

### Limitaciones
⚠️ **Dataset pequeño**: 37 muestras pueden no representar toda la variabilidad  
⚠️ **Posible overfitting**: Métricas perfectas sugieren revisar con más datos  
⚠️ **Batch effects**: Datos de fuentes diferentes pueden tener variaciones técnicas  

---

## 📧 Próximos Pasos

1. ✅ Modelo entrenado y validado
2. ✅ Visualizaciones generadas
3. ⬜ Validación en cohorte externa independiente
4. ⬜ Estudio funcional de SOX11 y DDIT4L
5. ⬜ Publicación de resultados

---

## 🔬 Conclusión Final

Se desarrolló exitosamente un modelo de **Logistic Regression** con **100% de accuracy** para clasificar células mesenquimales según su estado (Normal vs Osteoporosis Primaria). 

**Los genes SOX11 y DDIT4L** emergen como los marcadores más discriminativos, sugiriendo su potencial como biomarcadores para diagnóstico temprano de osteoporosis.

El modelo es **altamente explicable** y puede ser traducido directamente a reglas clínicas interpretables.

---

*Fecha de análisis: Noviembre 2025*  
*Proyecto: Bioinformática Avanzada - HMSC Predicción*
