import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold, cross_validate, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, confusion_matrix
import joblib
import os

# 1. Load Data
print("Loading data...")
X_path = 'X_features_matrix.csv'
y_path = 'y_labels.csv'

if not os.path.exists(X_path) or not os.path.exists(y_path):
    print(f"Error: Files not found in {os.getcwd()}")
    exit()

X_df = pd.read_csv(X_path)
y_df = pd.read_csv(y_path)

# Merge on Sample ID to ensure alignment
data = pd.merge(X_df, y_df, on='Sample')
print(f"Data merged. Shape: {data.shape}")

# Separate Features and Target
X = data.drop(columns=['Sample', 'TrueLabel'])
y = data['TrueLabel']
sample_ids = data['Sample']

# 2. Preprocessing
print("Scaling features...")
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
X_scaled = pd.DataFrame(X_scaled, columns=X.columns)

# 3. Define Models
models = {
    "Logistic Regression": LogisticRegression(random_state=42, solver='liblinear'),
    "Random Forest": RandomForestClassifier(random_state=42, n_estimators=100),
    "SVM (RBF)": SVC(probability=True, random_state=42),
    "SVM (Linear)": SVC(kernel='linear', probability=True, random_state=42),
    "SVM (Poly)": SVC(kernel='poly', degree=3, probability=True, random_state=42,C=1),
    "Decision Tree": DecisionTreeClassifier(random_state=42, max_depth=5)
}

# 4. Evaluation Strategy (Stratified K-Fold)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

scoring = {
    'accuracy': 'accuracy',
    'precision': 'precision',
    'recall': 'recall',
    'f1': 'f1',
    'roc_auc': 'roc_auc'
}

results = []

print("\nStarting Cross-Validation (5-Fold)...")
print("-" * 80)
print(f"{'Model':<20} | {'Acc':<8} | {'Prec':<8} | {'Recall':<8} | {'F1':<8} | {'AUC':<8}")
print("-" * 80)

best_model_name = ""
best_score = 0
best_model_obj = None

for name, model in models.items():
    scores = cross_validate(model, X_scaled, y, cv=cv, scoring=scoring)
    
    mean_acc = np.mean(scores['test_accuracy'])
    mean_prec = np.mean(scores['test_precision'])
    mean_rec = np.mean(scores['test_recall'])
    mean_f1 = np.mean(scores['test_f1'])
    mean_auc = np.mean(scores['test_roc_auc'])
    
    results.append({
        "Model": name,
        "Accuracy": mean_acc,
        "Precision": mean_prec,
        "Recall": mean_rec,
        "F1": mean_f1,
        "AUC": mean_auc
    })
    
    print(f"{name:<20} | {mean_acc:.3f}    | {mean_prec:.3f}    | {mean_rec:.3f}    | {mean_f1:.3f}    | {mean_auc:.3f}")

    # Track best model based on AUC (or F1)
    if mean_auc > best_score:
        best_score = mean_auc
        best_model_name = name
        best_model_obj = model

print("-" * 80)

# 5. Save Best Model (Retrained on full data)
print(f"\nBest Model: {best_model_name} (AUC: {best_score:.3f})")

# Confusion Matrix
print(f"\nConfusion Matrix for {best_model_name} (Aggregated via CV):")
y_pred_cv = cross_val_predict(best_model_obj, X_scaled, y, cv=cv)
cm = confusion_matrix(y, y_pred_cv)
print(cm)
print(f"TN: {cm[0][0]}, FP: {cm[0][1]}")
print(f"FN: {cm[1][0]}, TP: {cm[1][1]}")

print(f"Retraining {best_model_name} on full dataset and saving...")

best_model_obj.fit(X_scaled, y)
joblib.dump(best_model_obj, 'best_model_msc_senescence.pkl')
joblib.dump(scaler, 'scaler_msc_senescence.pkl')
print("Model and Scaler saved successfully.")

import matplotlib.pyplot as plt

# 6. Feature Importance (if applicable)
if best_model_name == "Random Forest":
    importances = best_model_obj.feature_importances_
    indices = np.argsort(importances)[::-1]
    
    top_n = 10
    top_indices = indices[:top_n]
    top_importances = importances[top_indices]
    top_features = X.columns[top_indices]

    print(f"\nTop {top_n} Feature Importances:")
    for f in range(top_n):
        print(f"{top_features[f]}: {top_importances[f]:.4f}")
    
    # Plotting
    plt.figure(figsize=(10, 6))
    plt.barh(range(top_n), top_importances, align='center', color='skyblue')
    plt.yticks(range(top_n), top_features)
    plt.xlabel('Feature Importance')
    plt.title('Top 10 Genes - Random Forest Importance')
    plt.gca().invert_yaxis()  # Highest importance at top
    plt.tight_layout()
    plt.savefig('feature_importance.png')
    print("Feature importance plot saved as 'feature_importance.png'")

elif best_model_name == "Logistic Regression":
    importances = np.abs(best_model_obj.coef_[0])
    indices = np.argsort(importances)[::-1]
    
    top_n = 10
    top_indices = indices[:top_n]
    top_importances = importances[top_indices]
    top_features = X.columns[top_indices]

    print(f"\nTop {top_n} Coefficients (Absolute Value):")
    for f in range(top_n):
        print(f"{top_features[f]}: {top_importances[f]:.4f}")

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.barh(range(top_n), top_importances, align='center', color='salmon')
    plt.yticks(range(top_n), top_features)
    plt.xlabel('Absolute Coefficient Value')
    plt.title('Top 10 Genes - Logistic Regression Coefficients')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig('feature_importance.png')
    print("Feature importance plot saved as 'feature_importance.png'")

elif best_model_name == "Decision Tree":
    importances = best_model_obj.feature_importances_
    indices = np.argsort(importances)[::-1]
    print("\nTop 10 Feature Importances:")
    for f in range(10):
        print(f"{X.columns[indices[f]]}: {importances[indices[f]]:.4f}")

