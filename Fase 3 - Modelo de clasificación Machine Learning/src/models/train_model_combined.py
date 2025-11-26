import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score, cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_auc_score
import joblib
import os
import matplotlib.pyplot as plt

# Paths
DATA_DIR = "c:/Users/eduar/OneDrive/Escritorio/bioinformatica_ML/ML_data/Fase 3 - Modelo de clasificación Machine Learning"
VALIDATION_DIR = os.path.join(DATA_DIR, "validation")

# Original Data
X_TRAIN_PATH = os.path.join(DATA_DIR, "X_features_matrix.csv")
Y_TRAIN_PATH = os.path.join(DATA_DIR, "y_labels.csv")

# Validation Data
X_VAL_PATH = os.path.join(VALIDATION_DIR, "X_test_GSE35958 (1).csv")
Y_VAL_PATH = os.path.join(VALIDATION_DIR, "y_test_GSE35958 (1) copy.csv")

def load_original_data():
    print("Loading original training data...")
    X_df = pd.read_csv(X_TRAIN_PATH)
    y_df = pd.read_csv(Y_TRAIN_PATH)
    
    # Merge on Sample ID
    data = pd.merge(X_df, y_df, on='Sample')
    return data

def load_validation_data(training_features):
    print(f"Loading validation data from {X_VAL_PATH}...")
    X_raw = pd.read_csv(X_VAL_PATH)
    
    # Transpose if needed (Genes as rows -> Samples as rows)
    if "Gene" in X_raw.columns:
        X_raw = X_raw.set_index("Gene")
    X_val = X_raw.T
    
    y_val = pd.read_csv(Y_VAL_PATH)
    
    # Filter X to match y samples
    common_samples = list(set(X_val.index) & set(y_val['Sample']))
    print(f"Validation samples found: {len(common_samples)}")
    
    X_val = X_val.loc[common_samples]
    y_val = y_val[y_val['Sample'].isin(common_samples)]
    
    # Align features
    X_aligned = pd.DataFrame(0, index=X_val.index, columns=training_features)
    common_features = list(set(X_val.columns) & set(training_features))
    X_aligned[common_features] = X_val[common_features]
    
    # Log transform if needed (simple heuristic)
    if X_aligned.max().max() > 50:
        print("Applying log2 transformation to validation data...")
        X_aligned = np.log2(X_aligned + 1)
        
    # Add Label column to merge
    # Ensure y_val is indexed by Sample to map correctly
    y_val_map = y_val.set_index('Sample')['Label']
    X_aligned['TrueLabel'] = X_aligned.index.map(y_val_map)
    X_aligned['Sample'] = X_aligned.index
    
    return X_aligned

def main():
    # 1. Load Original Data
    original_data = load_original_data()
    print(f"Original data shape: {original_data.shape}")
    
    # Get feature names (excluding Sample and TrueLabel)
    feature_cols = [c for c in original_data.columns if c not in ['Sample', 'TrueLabel']]
    
    # 2. Load and Align Validation Data
    val_data = load_validation_data(feature_cols)
    print(f"Validation data shape: {val_data.shape}")
    
    # 3. Merge Datasets
    combined_data = pd.concat([original_data, val_data], axis=0, ignore_index=True)
    print(f"Combined data shape: {combined_data.shape}")
    
    # Check class distribution
    print("\nClass Distribution in Combined Data:")
    print(combined_data['TrueLabel'].value_counts())
    
    # 4. Prepare X and y
    X = combined_data[feature_cols]
    y = combined_data['TrueLabel']
    samples = combined_data['Sample']
    
    # 5. Cross-Validation Strategy
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    # Scaling (must be done inside CV loop to avoid leakage, but for simplicity with sklearn cross_validate we can use a Pipeline or scale beforehand if we accept minor leakage. 
    # Better practice: Use Pipeline.
    from sklearn.pipeline import Pipeline
    
    # 6. Define Models with Pipelines
    models = {
        "Logistic Regression": Pipeline([
            ('scaler', StandardScaler()),
            ('clf', LogisticRegression(random_state=42, solver='liblinear'))
        ]),
        "Random Forest": Pipeline([
            ('scaler', StandardScaler()),
            ('clf', RandomForestClassifier(random_state=42, n_estimators=100))
        ]),
        "SVM (Linear)": Pipeline([
            ('scaler', StandardScaler()),
            ('clf', SVC(kernel='linear', probability=True, random_state=42))
        ])
    }
    
    scoring = ['accuracy', 'precision', 'recall', 'f1', 'roc_auc']
    
    print("\n" + "="*80)
    print("CROSS-VALIDATION RESULTS (5-FOLD)")
    print("="*80)
    print(f"{'Model':<20} | {'Acc':<8} | {'Prec':<8} | {'Recall':<8} | {'F1':<8} | {'AUC':<8}")
    print("-" * 80)
    
    best_score = 0
    best_model_name = ""
    best_pipeline = None
    
    for name, pipeline in models.items():
        scores = cross_val_score(pipeline, X, y, cv=cv, scoring='accuracy') # Simple check first
        
        # Detailed cross_validate
        cv_results = cross_validate(pipeline, X, y, cv=cv, scoring=scoring)
        
        mean_acc = np.mean(cv_results['test_accuracy'])
        mean_prec = np.mean(cv_results['test_precision'])
        mean_rec = np.mean(cv_results['test_recall'])
        mean_f1 = np.mean(cv_results['test_f1'])
        mean_auc = np.mean(cv_results['test_roc_auc'])
        
        print(f"{name:<20} | {mean_acc:.3f}    | {mean_prec:.3f}    | {mean_rec:.3f}    | {mean_f1:.3f}    | {mean_auc:.3f}")
        
        if mean_auc > best_score:
            best_score = mean_auc
            best_model_name = name
            best_pipeline = pipeline

    print("-" * 80)
    print(f"\nBEST MODEL (by AUC): {best_model_name} (AUC: {best_score:.3f})")
    
    # 7. Retrain Best Model on Full Data
    print(f"Retraining {best_model_name} on ALL combined data...")
    best_pipeline.fit(X, y)
    
    # Extract model and scaler from pipeline for saving
    final_model = best_pipeline.named_steps['clf']
    final_scaler = best_pipeline.named_steps['scaler']
    
    joblib.dump(final_model, os.path.join(DATA_DIR, 'best_model_combined.pkl'))
    joblib.dump(final_scaler, os.path.join(DATA_DIR, 'scaler_combined.pkl'))
    print("Model and scaler saved.")

if __name__ == "__main__":
    main()
