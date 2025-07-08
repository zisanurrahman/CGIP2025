# CIMPLE-Seq: Exploring Chemical-Genetic Interactions in *Burkholderia cenocepacia* K56-2

This repository contains the data analysis scripts and a Shiny web application developed to explore the chemical-genetic interaction profiles of *Burkholderia cenocepacia* K56-2 essential gene knockdown mutants. The interaction profiles were generated using the **CIMPLE-Seq** (CRISPRi-based Interaction Mapping via Pooled Library Enrichment Sequencing) approach.

> 📌 Reference:  
> **Rahman ASMZ et al.** (2024). *Rationally designed pooled CRISPRi-seq uncovers an inhibitor of bacterial peptidyl-tRNA hydrolase*.  
> Cell Reports, Volume 43, Issue 11, 114967

---

## 🧪 Project Summary

- CRISPRi was used to systematically knock down essential genes in *B. cenocepacia* K56-2.
- Mutant pools were screened against ~5000 compounds.
- CIMPLE-Seq was used to quantify mutant abundance post-treatment and compute **chemical-genetic interaction profiles (CGIPs)**.
- These CGIPs reveal mechanistic insights into compound action and essential gene function.

---

## 📂 Repository Structure

- `data/`  
  - Processed CGIP matrices and metadata  
  - Annotations and reference gene/compound information  

- `main/`  
  - Scripts for CGIP normalization, clustering, dimensionality reduction (t-SNE), and enrichment analysis

- `main/`  
  - Interactive R Shiny app to explore CGIPs, visualize compound-gene associations, and search by compound or gene name

---

## 🚀 Shiny App Features

- Explore gene-compound interaction landscapes  
- Visualize CGIP similarity using t-SNE/UMAP plots  
- Filter by compound, gene target, or mechanism of action  
- Export custom plots and result tables

