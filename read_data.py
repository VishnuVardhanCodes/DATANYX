import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from fpdf import FPDF
import io

# -----------------------------------------------------
# Load cleaned dataset + clusters
# -----------------------------------------------------

@st.cache_data
def load_data():
    df = pd.read_csv("clean_GSE64456.csv")
    df = df.rename(columns={df.columns[0]: "Gene"}).set_index("Gene")
    df = df.fillna(0)
    return df

df = load_data()
X = df.T
kmeans_labels = np.load("cluster_labels.npy")

# Train RF model
rf = RandomForestClassifier(n_estimators=200, random_state=42)
rf.fit(X, kmeans_labels)

# Pre-calc top variable genes
variances = df.var(axis=1)
top_genes = variances.sort_values(ascending=False).head(50)
top50_gene_names = list(top_genes.index)

# -----------------------------------------------------
# Streamlit UI
# -----------------------------------------------------

st.title("🧬 Immunodeficiency Gene Expression Analyzer")
uploaded = st.file_uploader("📤 Upload Patient Gene Expression CSV", type=["csv"])

# -----------------------------------------------------
# Process Uploaded File (SAFE)
# -----------------------------------------------------

user_df = None
if uploaded:
    try:
        temp_df = pd.read_csv(uploaded)

        if temp_df.empty:
            st.error("❌ Empty CSV uploaded.")
            st.stop()

        first_col = temp_df.columns[0]
        temp_df = temp_df.rename(columns={first_col: "Gene"}).set_index("Gene")
        temp_df = temp_df.reindex(df.index).fillna(0)

        user_df = temp_df.copy()

    except Exception:
        st.error("❌ Invalid CSV format.")
        st.stop()

# -----------------------------------------------------
# PCA Plot
# -----------------------------------------------------

if user_df is not None and st.checkbox("📊 Show PCA Plot (Patient Highlighted)"):

    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(df.T)
    patient_pca = pca.transform(user_df.T)

    plt.figure(figsize=(8, 6))
    plt.scatter(pca_result[:, 0], pca_result[:, 1], c=kmeans_labels,
                cmap="coolwarm", alpha=0.5)
    plt.scatter(patient_pca[0][0], patient_pca[0][1], s=200, color="yellow",
                marker="X", edgecolors="black", linewidths=1.5)
    plt.title("PCA Plot (Patient Highlighted)")
    st.pyplot(plt)

# -----------------------------------------------------
# Dataset Heatmap
# -----------------------------------------------------

if st.checkbox("🧪 Show Dataset Heatmap (Top 50 Genes)"):

    df_top50 = df.loc[top50_gene_names]
    plt.figure(figsize=(12, 10))
    sns.heatmap(df_top50, cmap="viridis")
    plt.title("Top 50 Most Variable Genes (Dataset)")
    st.pyplot(plt)

# -----------------------------------------------------
# Prediction Section
# -----------------------------------------------------

if user_df is not None:

    prediction = rf.predict(user_df.T)[0]

    st.subheader("🔮 Prediction Result")
    if prediction == 0:
        st.success("🟢 Healthy Gene Expression Pattern")
    else:
        st.error("🔴 Possible Immunodeficiency Pattern")

    risk = rf.predict_proba(user_df.T)[0][1] * 100
    st.metric("Immune Dysfunction Risk", f"{risk:.2f}%")

    if risk < 30:
        st.success("🟢 Low Risk")
    elif risk < 60:
        st.warning("🟡 Moderate Risk")
    else:
        st.error("🔴 High Risk")

    # Similar Patient Search
    st.subheader("🔍 Closest Matching Patient")
    similarities = cosine_similarity(user_df.T, df.T)[0]
    idx = similarities.argmax()
    st.write("Closest Sample:", df.columns[idx])
    st.write("Similarity Score:", round(similarities[idx], 3))

# -----------------------------------------------------
# Patient Heatmap
# -----------------------------------------------------

if user_df is not None and st.checkbox("🧬 Show Patient Heatmap (Top 50 Genes)"):

    df_user_top50 = user_df.loc[top50_gene_names]

    plt.figure(figsize=(8, 6))
    sns.heatmap(df_user_top50, cmap="coolwarm")
    plt.title("Patient Gene Heatmap (Top 50 Genes)")
    st.pyplot(plt)

# -----------------------------------------------------
# Risk Timeline Simulation
# -----------------------------------------------------

if user_df is not None and st.checkbox("📈 Show Patient Risk Timeline (Simulated)"):

    timeline = np.linspace(max(5, risk - 20), risk, 12)

    plt.figure(figsize=(8, 5))
    plt.plot(timeline, marker="o", linewidth=2)
    plt.title("Patient Risk Timeline (Simulated 12 Months)")
    plt.xlabel("Month")
    plt.ylabel("Risk Score %")
    st.pyplot(plt)

# -----------------------------------------------------
# PDF Report Generator
# -----------------------------------------------------

if user_df is not None and st.checkbox("📄 Download PDF Report"):

    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)

    pdf.cell(0, 10, "Immunodeficiency Gene Expression Report", ln=True)
    pdf.ln(5)
    pdf.cell(0, 10, f"Prediction: {prediction}", ln=True)
    pdf.cell(0, 10, f"Risk Score: {risk:.2f}%", ln=True)
    pdf.cell(0, 10, f"Closest Patient Match: {df.columns[idx]}", ln=True)

    buffer = io.BytesIO()
    pdf.output(buffer)
    buffer.seek(0)

    st.download_button(
        label="⬇ Download PDF",
        data=buffer,
        file_name="Patient_Report.pdf",
        mime="application/pdf",
    )

    st.success("PDF Generated Successfully!")
