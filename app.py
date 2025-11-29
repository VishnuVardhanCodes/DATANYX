import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import shap
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

# Train RF using saved clusters
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
# Process Uploaded File
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

    except:
        st.error("❌ Invalid CSV format.")
        st.stop()

# -----------------------------------------------------
# PCA Plot
# -----------------------------------------------------

if user_df is not None and st.checkbox("📊 Show PCA Plot"):
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(df.T)
    patient_pca = pca.transform(user_df.T)

    plt.figure(figsize=(8, 6))
    plt.scatter(pca_result[:, 0], pca_result[:, 1], c=kmeans_labels, cmap="coolwarm", alpha=0.5)
    plt.scatter(patient_pca[0][0], patient_pca[0][1],
                s=180, color="yellow", marker="X", edgecolors="black")
    plt.title("PCA: Patient Highlighted")
    st.pyplot(plt)

# -----------------------------------------------------
# Dataset Heatmap
# -----------------------------------------------------

if st.checkbox("🧪 Show Dataset Heatmap (Top 50 Genes)"):
    df_top50 = df.loc[top50_gene_names]
    plt.figure(figsize=(12, 10))
    sns.heatmap(df_top50, cmap="viridis")
    st.pyplot(plt)

# -----------------------------------------------------
# Prediction
# -----------------------------------------------------

if user_df is not None:

    prediction = rf.predict(user_df.T)[0]

    st.subheader("🔮 Prediction Result")

    if prediction == 0:
        st.success("🟢 Healthy Expression Pattern")
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

    # Nearest patient
    similarities = cosine_similarity(user_df.T, df.T)[0]
    idx = similarities.argmax()
    st.subheader("🔍 Closest Patient in Dataset")
    st.write("Closest Sample:", df.columns[idx])
    st.write("Similarity Score:", round(similarities[idx], 3))

# -----------------------------------------------------
# Patient Heatmap
# -----------------------------------------------------

if user_df is not None and st.checkbox("🧬 Patient Gene Heatmap (Top 50)"):
    df_user_top50 = user_df.loc[top50_gene_names]
    plt.figure(figsize=(9, 6))
    sns.heatmap(df_user_top50, cmap="coolwarm")
    st.pyplot(plt)

# -----------------------------------------------------
# SHAP Explainability
# -----------------------------------------------------

if user_df is not None and st.checkbox("🧠 Show SHAP Explainability"):

    st.write("Calculating SHAP values (may take 10 sec)...")

    explainer = shap.TreeExplainer(rf)
    shap_values = explainer.shap_values(user_df.T)

    # Plot SHAP summary
    fig, ax = plt.subplots(figsize=(8, 5))
    shap.summary_plot(shap_values[1], user_df.T, plot_type="bar", show=False)
    st.pyplot(fig)

    st.success("Top SHAP genes explain why ML predicted this risk.")

# -----------------------------------------------------
# Patient Risk Trend Animation (Fake Clinical Timeline)
# -----------------------------------------------------

if user_df is not None and st.checkbox("📈 Generate Patient Risk Timeline Animation"):

    timeline = np.linspace(max(5, risk - 20), risk, 12)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(timeline, marker="o")
    ax.set_title("Patient Risk Trend (Simulated Clinical Follow-up)")
    ax.set_xlabel("Month")
    ax.set_ylabel("Risk Score %")
    st.pyplot(fig)

    st.info("This simulates how risk may evolve across 12 months.")

# -----------------------------------------------------
# Doctor PDF Report Generator
# -----------------------------------------------------

if user_df is not None and st.checkbox("📄 Download Doctor-Friendly PDF Report"):

    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)

    pdf.cell(0, 10, "Immunodeficiency Gene Expression Report", ln=True)
    pdf.ln(5)

    pdf.cell(0, 10, f"Predicted Cluster: {prediction}", ln=True)
    pdf.cell(0, 10, f"Risk Score: {risk:.2f}%", ln=True)
    pdf.ln(5)

    pdf.cell(0, 10, f"Closest Matching Sample: {df.columns[idx]}", ln=True)

    # Save PDF to memory
    buffer = io.BytesIO()
    pdf.output(buffer)
    buffer.seek(0)

    st.download_button(
        label="⬇ Download PDF Report",
        data=buffer,
        file_name="Patient_Report.pdf",
        mime="application/pdf",
    )
    st.success("PDF Generated Successfully!")
