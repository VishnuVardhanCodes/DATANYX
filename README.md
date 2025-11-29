🧬 Immunodeficiency Gene Expression Analyzer

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

🚀 Features

📤 Upload patient gene expression files

🔮 Predict immune dysfunction using ML

📊 PCA plot with patient highlighted

🧪 Dataset heatmap (Top 50 variable genes)

🧬 Patient-specific gene expression heatmap

📈 Simulated patient risk timeline

📄 Downloadable doctor-style PDF report

🔍 Similarity search using cosine similarity

⚡ Fast, accurate, and hackathon-ready

🛠️ Technologies Used

Python

Streamlit

Pandas

NumPy

Scikit-learn

Seaborn + Matplotlib

FPDF

Cosine Similarity ML techniques

Gene Expression Dataset (GSE64456 – preprocessed)

📂 Project Structure
DATANYX/
│── app.py                   # Main Streamlit application
│── analyze_data.py          # Dataset preprocessing + PCA + clustering
│── clean_GSE64456.csv       # Cleaned dataset (after preprocessing)
│── cluster_labels.npy       # Saved KMeans cluster labels
│── patient_sample.csv       # Example patient file
│── PCA_plot.png             # PCA visualization output
│── Heatmap_top50.png        # Heatmap visualization output
│── README.md                # Project documentation

🧪 How to Run the Project
1. Clone the Repository
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name

2. Create Virtual Environment
python -m venv venv

3. Activate Virtual Environment

Windows:

venv\Scripts\activate


Mac/Linux:

source venv/bin/activate

4. Install Dependencies
pip install -r requirements.txt


(If you want I can generate the full requirements.txt for you.)

🎯 Running the Dashboard
streamlit run app.py


Open browser →
👉 http://localhost:8501/

📤 Uploading Gene Expression CSV

Your CSV must follow this format:

Gene	GSMxxxx1	GSMxxxx2	...
ILMN_12345	120	240	...
ILMN_54321	560	180	...

If needed, you can use the provided patient_sample.csv.

📄 PDF Report

The system generates a downloadable clinical-style PDF summarizing:

Predicted cluster

Immune dysfunction risk score

Closest database patient

Gene activity summary

Perfect for doctors, researchers, or hackathon demo judges.

⭐ Why This Project Stands Out

Full ML pipeline

Medical-grade visualizations

Interactive dashboard

Auto PDF generation

Works on real gene expression datasets

Professional, clean-looking UI

Perfect for hackathons & research portfolios

🤝 Contributing

Feel free to submit Pull Requests, report issues, or suggest new features!

📜 License

MIT License