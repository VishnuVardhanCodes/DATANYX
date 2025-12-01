🧬 **Immunodeficiency Gene Expression Analyzer**

A powerful Streamlit-based AI tool that analyzes patient gene expression profiles, identifies immunodeficiency patterns, visualizes top genes, compares against known patient clusters, and generates medical-style PDF reports.

This project uses:

Machine Learning (Random Forest)

Dimensionality Reduction (PCA)

Gene Variability Analysis

Heatmaps for Gene Visualization

Risk Scoring

Closest-Patient Matching

PDF Report Generation

and an Interactive Web Dashboard (Streamlit)

**🚀 Features**

📤 Upload patient gene expression files

🔮 Predict immune dysfunction using ML

📊 PCA plot with patient highlighted

🧪 Dataset heatmap (Top 50 variable genes)

🧬 Patient-specific gene expression heatmap

📈 Simulated patient risk timeline

📄 Downloadable doctor-style PDF report

🔍 Similarity search using cosine similarity

⚡ Fast, accurate, and hackathon-ready

**🛠️ Technologies Used**

Python

Streamlit

Pandas

NumPy

Scikit-learn

Seaborn + Matplotlib

FPDF

Cosine Similarity ML techniques

Gene Expression Dataset (GSE64456 – preprocessed)

📂 **Project Structure**
DATANYX/
│── app.py                   # Main Streamlit application
│── analyze_data.py          # Dataset preprocessing + PCA + clustering
│── clean_GSE64456.csv       # Cleaned dataset (after preprocessing)
│── cluster_labels.npy       # Saved KMeans cluster labels
│── patient_sample.csv       # Example patient file
│── PCA_plot.png             # PCA visualization output
│── Heatmap_top50.png        # Heatmap visualization output
│── README.md                # Project documentation


**⭐ Why This Project Stands Out**

Full ML pipeline

Medical-grade visualizations

Interactive dashboard

Auto PDF generation

Works on real gene expression datasets

Professional, clean-looking UI

Perfect for hackathons & research portfolios

**🤝 Contributing**

Feel free to submit Pull Requests, report issues, or suggest new features!

**📜 License**

MIT License

