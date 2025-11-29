import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans

# -----------------------------------------
# STEP 1: Load the cleaned dataset
# -----------------------------------------

df = pd.read_csv("clean_GSE64456.csv")
print("Original shape:", df.shape)

# -----------------------------------------
# STEP 2: Clean gene names and set index
# -----------------------------------------

df = df.dropna(subset=[df.columns[0]])   # remove missing gene names
df = df.rename(columns={df.columns[0]: "Gene"})
df = df.set_index("Gene")

print("\nAfter setting index:", df.shape)
print(df.head())

# -----------------------------------------
# STEP 3: Remove NaNs
# -----------------------------------------

df = df.fillna(0)
print("\nAfter removing NaN:", df.shape)

# -----------------------------------------
# STEP 4: PCA Visualization (SAVED)
# -----------------------------------------

print("\nRunning PCA...")

pca = PCA(n_components=2)
pca_result = pca.fit_transform(df.T)

plt.figure(figsize=(8, 6))
plt.scatter(pca_result[:, 0], pca_result[:, 1])
plt.title("PCA of Samples")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.grid(True)
plt.savefig("PCA_plot.png", dpi=300)
plt.close()

print("PCA plot saved as PCA_plot.png")

# -----------------------------------------
# STEP 5: Feature Selection - Top 100 variable genes
# -----------------------------------------

print("\nSelecting top variable genes...")

gene_variances = df.var(axis=1)
top_genes = gene_variances.sort_values(ascending=False).head(100)

print("\nTop 20 variable genes:")
print(top_genes.head(20))

# -----------------------------------------
# STEP 6: Heatmap of Top 50 Genes (SAVED)
# -----------------------------------------

print("\nPreparing top 50 gene matrix...")

top50_gene_names = list(top_genes.index[:50])
df_top50 = df.loc[top50_gene_names]

print("Top 50 gene matrix shape:", df_top50.shape)

print("\nGenerating heatmap...")

plt.figure(figsize=(12, 10))
sns.heatmap(df_top50, cmap='viridis')
plt.title("Heatmap of Top 50 Highly Variable Genes")
plt.savefig("Heatmap_top50.png", dpi=300)
plt.close()

print("Heatmap saved as Heatmap_top50.png")

# -----------------------------------------
# STEP 7: KMeans Clustering
# -----------------------------------------

print("\nRunning KMeans clustering...")

kmeans = KMeans(n_clusters=2, random_state=42)
clusters = kmeans.fit_predict(df.T)

print("\nCluster labels for samples:")
print(clusters)

# -----------------------------------------
# STEP 8: MACHINE LEARNING MODEL (Random Forest)
# -----------------------------------------

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.model_selection import train_test_split

print("\nTraining Random Forest model...")

# X = gene expression matrix (samples as rows)
X = df.T

# y = cluster labels from KMeans
y = clusters

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Random Forest classifier
rf = RandomForestClassifier(n_estimators=200, random_state=42)
rf.fit(X_train, y_train)

# Predictions
y_pred = rf.predict(X_test)

# Evaluation
print("\nRandom Forest Accuracy:", accuracy_score(y_test, y_pred))
print("\nClassification Report:\n", classification_report(y_test, y_pred))
print("\nConfusion Matrix:\n", confusion_matrix(y_test, y_pred))

# -----------------------------------------
# STEP 9: Feature Importance from Random Forest
# -----------------------------------------

importances = pd.Series(rf.feature_importances_, index=X.columns)
top20_rf_genes = importances.sort_values(ascending=False).head(20)

print("\nTop 20 predictive genes:")
print(top20_rf_genes)


import numpy as np
np.save("cluster_labels.npy", clusters)
