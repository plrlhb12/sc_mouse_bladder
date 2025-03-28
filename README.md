
# 🧬 Seurat Analysis: WT vs HDAC11 KO Mouse Bladder

This repository contains the single-cell RNA-seq analysis pipeline for comparing wild-type (WT) and HDAC11 knockout (KO) mouse bladder samples using [Seurat](https://satijalab.org/seurat/).

## 📁 Project Structure

```
seurat_HDAC11_mouse_bladder/
├── data/               # Raw 10X Genomics output (filtered_feature_bc_matrix)
├── scripts/            # R scripts for each analysis step
├── qc/                 # Quality control results and filtered objects
├── results/            # Plots, tables, and final DE results
├── notebooks/          # RMarkdown/Jupyter notebooks (optional)
├── README.md           # Project overview
├── .gitignore          # Ignore rules
```

## 🧪 Dataset Overview

- **Organism:** Mouse
- **Tissue:** Bladder
- **Conditions:** WT vs HDAC11 KO
- **Replicates:** 2 biological replicates per condition
- **Platform:** 10x Genomics (v3 chemistry)

## 🧭 Analysis Workflow

1. **Preprocessing and Quality Control**  
   - Filtering cells based on gene count, UMI count, and mitochondrial percentage  
   - Violin plots and scatterplots for quality inspection

2. **Normalization & Feature Selection**  
   - SCTransform or LogNormalize (configurable)  
   - Selection of variable features

3. **Dimensionality Reduction & Clustering**  
   - PCA, UMAP  
   - Clustering using shared nearest neighbor (SNN)

4. **Cell Type Annotation**  
   - Marker-based or automated annotation with reference datasets

5. **Differential Gene Expression (DGE)**  
   - Between WT and KO for each major cell type

6. **Visualization**  
   - UMAPs, violin plots, dot plots, heatmaps

## 🔧 Setup Instructions

1. Clone the repo:
   ```bash
   git clone https://github.com/YOUR_USERNAME/seurat-hdac11-mouse-bladder.git
   cd seurat-hdac11-mouse-bladder
   ```

2. Load R (preferably ≥4.2) and install dependencies:
   ```r
   install.packages("Seurat")
   # or use renv if version-controlled environment is preferred
   ```

3. Organize your 10X outputs:
   ```
   data/
   ├── WT_1/outs/filtered_feature_bc_matrix/
   ├── WT_2/outs/filtered_feature_bc_matrix/
   ├── KO_1/outs/filtered_feature_bc_matrix/
   ├── KO_2/outs/filtered_feature_bc_matrix/
   ```

4. Run scripts in order from `scripts/` folder.

## 🧾 License

MIT License.

## 👩‍🔬 Acknowledgments

- Developed by [Your Name]
- Based on the framework provided by the [Seurat team](https://satijalab.org/seurat/)
- Funded by

