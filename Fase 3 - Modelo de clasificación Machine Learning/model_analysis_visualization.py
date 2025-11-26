import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier, plot_tree
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
import joblib
import os

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
    data = pd.merge(X_df, y_df, on='Sample')
    return data

def load_validation_data(training_features):
    print(f"Loading validation data from {X_VAL_PATH}...")
    X_raw = pd.read_csv(X_VAL_PATH)
    
    if "Gene" in X_raw.columns:
        X_raw = X_raw.set_index("Gene")
    X_val = X_raw.T
    
    y_val = pd.read_csv(Y_VAL_PATH)
    
    common_samples = list(set(X_val.index) & set(y_val['Sample']))
    print(f"Validation samples found: {len(common_samples)}")
    
    X_val = X_val.loc[common_samples]
    y_val = y_val[y_val['Sample'].isin(common_samples)]
    
    X_aligned = pd.DataFrame(0, index=X_val.index, columns=training_features)
    common_features = list(set(X_val.columns) & set(training_features))
    X_aligned[common_features] = X_val[common_features]
    
    if X_aligned.max().max() > 50:
        print("Applying log2 transformation to validation data...")
        X_aligned = np.log2(X_aligned + 1)
        
    y_val_map = y_val.set_index('Sample')['Label']
    X_aligned['TrueLabel'] = X_aligned.index.map(y_val_map)
    X_aligned['Sample'] = X_aligned.index
    
    return X_aligned

def plot_decision_boundaries_pca(X, y, models_dict, feature_names):
    """Reduce to 2D with PCA and plot decision boundaries"""
    print("\n" + "="*80)
    print("CREATING DECISION BOUNDARY VISUALIZATIONS (PCA)")
    print("="*80)
    
    # Apply PCA to reduce to 2D
    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X)
    
    explained_var = pca.explained_variance_ratio_
    print(f"PCA explained variance: PC1={explained_var[0]:.3f}, PC2={explained_var[1]:.3f}")
    print(f"Total variance explained: {sum(explained_var):.3f}")
    
    # Create mesh grid
    h = 0.02
    x_min, x_max = X_pca[:, 0].min() - 1, X_pca[:, 0].max() + 1
    y_min, y_max = X_pca[:, 1].min() - 1, X_pca[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                         np.arange(y_min, y_max, h))
    
    # Plot for each model
    fig, axes = plt.subplots(2, 1, figsize=(15, 12))
    axes = axes.ravel()
    
    for idx, (name, model) in enumerate(models_dict.items()):
        # Train model on PCA-transformed data
        model.fit(X_pca, y)
        
        # Predict on mesh
        Z = model.predict(np.c_[xx.ravel(), yy.ravel()])
        Z = Z.reshape(xx.shape)
        
        # Plot
        axes[idx].contourf(xx, yy, Z, alpha=0.3, cmap='RdYlBu')
        scatter = axes[idx].scatter(X_pca[:, 0], X_pca[:, 1], c=y, 
                                   cmap='RdYlBu', edgecolors='black', s=100)
        axes[idx].set_xlabel(f'PC1 ({explained_var[0]:.1%} var)', fontsize=12)
        axes[idx].set_ylabel(f'PC2 ({explained_var[1]:.1%} var)', fontsize=12)
        axes[idx].set_title(f'{name}', fontsize=14, fontweight='bold')
        axes[idx].legend(*scatter.legend_elements(), title="Classes", loc='upper right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(DATA_DIR, 'decision_boundaries_pca.png'), dpi=300)
    print("Decision boundary plot saved as 'decision_boundaries_pca.png'")
    plt.close()

def visualize_decision_tree(tree_model, feature_names, max_depth=3):
    """Visualize decision tree"""
    print("\n" + "="*80)
    print("CREATING DECISION TREE VISUALIZATION")
    print("="*80)
    
    plt.figure(figsize=(20, 10))
    plot_tree(tree_model, 
              feature_names=feature_names,
              class_names=['Class 0 (Control)', 'Class 1 (OP)'],
              filled=True,
              rounded=True,
              fontsize=10,
              max_depth=max_depth)
    plt.title('Decision Tree Classifier (Simplified)', fontsize=16, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(DATA_DIR, 'decision_tree_visualization.png'), dpi=300, bbox_inches='tight')
    print("Decision tree visualization saved as 'decision_tree_visualization.png'")
    plt.close()

def plot_feature_importance_comparison(models_dict, feature_names):
    """Compare feature importance across models"""
    print("\n" + "="*80)
    print("CREATING FEATURE IMPORTANCE COMPARISON")
    print("="*80)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Logistic Regression coefficients
    lr_model = models_dict['Logistic Regression']
    lr_coef = np.abs(lr_model.coef_[0])
    n_features = min(10, len(lr_coef))
    lr_top_idx = np.argsort(lr_coef)[-n_features:]
    
    axes[0].barh(range(n_features), lr_coef[lr_top_idx], color='salmon')
    axes[0].set_yticks(range(n_features))
    axes[0].set_yticklabels([feature_names[i] for i in lr_top_idx])
    axes[0].set_xlabel('Absolute Coefficient Value', fontsize=12)
    axes[0].set_title(f'Logistic Regression - Top {n_features} Features', fontsize=14, fontweight='bold')
    axes[0].invert_yaxis()
    
    # Decision Tree feature importance
    dt_model = models_dict['Decision Tree']
    dt_imp = dt_model.feature_importances_
    dt_top_idx = np.argsort(dt_imp)[-n_features:]
    
    axes[1].barh(range(n_features), dt_imp[dt_top_idx], color='skyblue')
    axes[1].set_yticks(range(n_features))
    axes[1].set_yticklabels([feature_names[i] for i in dt_top_idx])
    axes[1].set_xlabel('Feature Importance', fontsize=12)
    axes[1].set_title(f'Decision Tree - Top {n_features} Features', fontsize=14, fontweight='bold')
    axes[1].invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(os.path.join(DATA_DIR, 'feature_importance_comparison.png'), dpi=300)
    print("Feature importance comparison saved as 'feature_importance_comparison.png'")
    plt.close()

def main():
    # 1. Load and combine data
    original_data = load_original_data()
    print(f"Original data shape: {original_data.shape}")
    
    feature_cols = [c for c in original_data.columns if c not in ['Sample', 'TrueLabel']]
    
    val_data = load_validation_data(feature_cols)
    print(f"Validation data shape: {val_data.shape}")
    
    combined_data = pd.concat([original_data, val_data], axis=0, ignore_index=True)
    print(f"Combined data shape: {combined_data.shape}")
    
    print("\nClass Distribution in Combined Data:")
    print(combined_data['TrueLabel'].value_counts())
    
    X = combined_data[feature_cols]
    y = combined_data['TrueLabel']
    
    # 2. Cross-Validation with Decision Tree added
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    models = {
        "Logistic Regression": Pipeline([
            ('scaler', StandardScaler()),
            ('clf', LogisticRegression(random_state=42, solver='liblinear'))
        ]),
        "Decision Tree": Pipeline([
            ('scaler', StandardScaler()),
            ('clf', DecisionTreeClassifier(random_state=42, max_depth=6, min_samples_split=3))
        ])
    }
    
    scoring = ['accuracy', 'precision', 'recall', 'f1', 'roc_auc']
    
    print("\n" + "="*80)
    print("CROSS-VALIDATION RESULTS (5-FOLD) - INCLUDING DECISION TREE")
    print("="*80)
    print(f"{'Model':<25} | {'Acc':<8} | {'Prec':<8} | {'Recall':<8} | {'F1':<8} | {'AUC':<8}")
    print("-" * 80)
    
    results = {}
    
    for name, pipeline in models.items():
        cv_results = cross_validate(pipeline, X, y, cv=cv, scoring=scoring)
        
        mean_acc = np.mean(cv_results['test_accuracy'])
        mean_prec = np.mean(cv_results['test_precision'])
        mean_rec = np.mean(cv_results['test_recall'])
        mean_f1 = np.mean(cv_results['test_f1'])
        mean_auc = np.mean(cv_results['test_roc_auc'])
        
        results[name] = {
            'accuracy': mean_acc,
            'precision': mean_prec,
            'recall': mean_rec,
            'f1': mean_f1,
            'auc': mean_auc
        }
        
        print(f"{name:<25} | {mean_acc:.3f}    | {mean_prec:.3f}    | {mean_rec:.3f}    | {mean_f1:.3f}    | {mean_auc:.3f}")
    
    print("-" * 80)
    
    # 3. Train models on full data for visualization
    print("\nTraining models on full dataset for visualization...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    models_for_viz = {
        'Logistic Regression': LogisticRegression(random_state=42, solver='liblinear'),
        'Decision Tree (depth=6)': DecisionTreeClassifier(random_state=42, max_depth=6),
    }
    
    for name, model in models_for_viz.items():
        model.fit(X_scaled, y)
    
    # 4. Create visualizations
    plot_decision_boundaries_pca(X_scaled, y, models_for_viz, feature_cols)
    
    # Visualize decision tree
    visualize_decision_tree(models_for_viz['Decision Tree (depth=6)'], feature_cols, max_depth=6)
    
    # Feature importance comparison
    comparison_models = {
        'Logistic Regression': models_for_viz['Logistic Regression'],
        'Decision Tree': models_for_viz['Decision Tree (depth=6)']
    }
    plot_feature_importance_comparison(comparison_models, feature_cols)
    
    # 5. Save summary statistics
    print("\n" + "="*80)
    print("SUMMARY - MODEL SELECTION")
    print("="*80)
    print("\nMost Explainable Model: Logistic Regression")
    print("Reason: Linear coefficients are directly interpretable")
    print("\nTop 5 Most Important Genes (Logistic Regression):")
    lr_coef = np.abs(models_for_viz['Logistic Regression'].coef_[0])
    top_5_idx = np.argsort(lr_coef)[-5:][::-1]
    for i, idx in enumerate(top_5_idx, 1):
        print(f"{i}. {feature_cols[idx]}: {lr_coef[idx]:.4f}")
    
    # 6. Save the explainable model
    joblib.dump(models_for_viz['Logistic Regression'], os.path.join(DATA_DIR, 'final_explainable_model.pkl'))
    joblib.dump(scaler, os.path.join(DATA_DIR, 'final_scaler.pkl'))
    print("\nFinal explainable model saved as 'final_explainable_model.pkl'")
    
    return results, feature_cols

if __name__ == "__main__":
    results, features = main()
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
