import pandas as pd
import joblib
import os
import numpy as np
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, roc_auc_score

# Paths
DATA_DIR = "c:/Users/eduar/OneDrive/Escritorio/bioinformatica_ML/ML_data/Fase 3 - Modelo de clasificación Machine Learning"
VALIDATION_DIR = os.path.join(DATA_DIR, "validation")
MODEL_PATH = os.path.join(DATA_DIR, "best_model_msc_senescence.pkl")
SCALER_PATH = os.path.join(DATA_DIR, "scaler_msc_senescence.pkl")
TRAIN_FEATURES_PATH = os.path.join(DATA_DIR, "X_features_matrix.csv")

# New Validation Files
VALIDATION_X = os.path.join(VALIDATION_DIR, "X_test_GSE35958 (1).csv")
VALIDATION_Y = os.path.join(VALIDATION_DIR, "y_test_GSE35958 (1) copy.csv")# es la data de validación positiva

def load_model_and_scaler():
    print("Loading model and scaler...")
    model = joblib.load(MODEL_PATH)
    scaler = joblib.load(SCALER_PATH)
    return model, scaler

def get_training_features():
    print("Loading training feature names...")
    # Read only the header
    df = pd.read_csv(TRAIN_FEATURES_PATH, nrows=0)
    features = df.columns.tolist()
    if "Sample" in features:
        features.remove("Sample")
    if "Unnamed: 0" in features:
        features.remove("Unnamed: 0")
    return features

def align_features(df, training_features):
    print("Aligning features...")
    # Create a DataFrame with the training features, filled with 0
    df_aligned = pd.DataFrame(0, index=df.index, columns=training_features)
    
    # Identify common features
    common_features = list(set(df.columns) & set(training_features))
    print(f"Number of common features: {len(common_features)}")
    
    if len(common_features) < len(training_features):
        missing = list(set(training_features) - set(df.columns))
        print(f"Missing features (first 10): {missing[:10]}")
        print(f"Total missing: {len(missing)}")
    
    # Update the aligned DataFrame with values from the input DataFrame
    df_aligned[common_features] = df[common_features]
    
    return df_aligned

def validate_new_data(model, scaler, training_features):
    print("\n--- Validating on New Data (GSE35958) ---")
    try:
        # Load data
        print(f"Loading X from {VALIDATION_X}")
        X_raw = pd.read_csv(VALIDATION_X)
        
        # Transpose: Genes are rows, Samples are columns in the raw file
        # We need Samples as rows, Genes as columns
        # Assuming first column is 'Gene'
        if "Gene" in X_raw.columns:
            X_raw = X_raw.set_index("Gene")
        X = X_raw.T # Transpose
        
        print(f"Loading y from {VALIDATION_Y}")
        y = pd.read_csv(VALIDATION_Y)
        
        # Ensure samples match
        # X index should be sample IDs (GSM...), y 'Sample' column has IDs
        # Filter X to match y samples and order
        # Check intersection
        common_samples = list(set(X.index) & set(y['Sample']))
        print(f"Common samples: {len(common_samples)}")
        
        X = X.loc[y['Sample']]
        y_true = y['Label']
        
        print(f"Loaded {X.shape[0]} samples.")
        
        # Align features
        X_aligned = align_features(X, training_features)
        
        print("Data stats before scaling:")
        print(X_aligned.describe().loc[['min', 'max', 'mean']].T.head())

        # Check if log transformation is needed
        # If max value is > 100, it's likely raw counts. If < 20, likely log2.
        if X_aligned.max().max() > 50:
             print("Max value > 50, applying log2(x + 1) transformation...")
             X_aligned = np.log2(X_aligned + 1)
             print("Data stats after log2 transformation:")
             print(X_aligned.describe().loc[['min', 'max', 'mean']].T.head())
        else:
            print("Max value <= 50, assuming data is already log-transformed.")

        # Scale features
        X_scaled = scaler.transform(X_aligned)
        
        # Predict
        y_pred = model.predict(X_scaled)
        y_prob = model.predict_proba(X_scaled)[:, 1]
        
        # Create detailed results DataFrame
        results_df = pd.DataFrame({
            'Sample': y['Sample'].values,
            'True_Label': y_true.values,
            'Predicted_Label': y_pred,
            'Prob_Class_0': 1 - y_prob,
            'Prob_Class_1': y_prob,
            'Correct': y_true.values == y_pred
        })
        
        # Evaluate
        print("\n" + "="*80)
        print("VALIDATION RESULTS")
        print("="*80)
        print(f"Accuracy: {accuracy_score(y_true, y_pred):.4f}")
        print(f"\nConfusion Matrix:")
        print(confusion_matrix(y_true, y_pred))
        print(f"\nClassification Report:")
        print(classification_report(y_true, y_pred))
        
        if len(y_true.unique()) > 1:
             print(f"ROC AUC: {roc_auc_score(y_true, y_prob):.4f}")
        else:
            print("ROC AUC: Not defined (only one class in y_true)")
        
        # Show all predictions
        print("\n" + "="*80)
        print("DETAILED PREDICTIONS FOR ALL SAMPLES")
        print("="*80)
        print(results_df.to_string(index=False))
        
        # Show misclassified samples
        misclassified = results_df[~results_df['Correct']]
        if len(misclassified) > 0:
            print("\n" + "="*80)
            print(f"[X] MISCLASSIFIED SAMPLES ({len(misclassified)}/{len(results_df)})")
            print("="*80)
            print(misclassified.to_string(index=False))
        else:
            print("\n" + "="*80)
            print("[OK] ALL SAMPLES CORRECTLY CLASSIFIED!")
            print("="*80)

    except Exception as e:
        print(f"Error validating new data: {e}")
        import traceback
        traceback.print_exc()

def main():
    model, scaler = load_model_and_scaler()
    training_features = get_training_features()
    
    validate_new_data(model, scaler, training_features)

if __name__ == "__main__":
    main()
