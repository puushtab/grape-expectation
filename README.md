# 🍇 Grape Expectation

A complete statistical learning pipeline to classify grape varieties from morphometric image features — based on the Raisin dataset (Cinar et al., 2020). Developed as the final project for the **STA203 Statistical Learning course** at ENSTA Paris.

---

## 🔬 Project Overview

**Objective:** Predict whether a grape belongs to the *Besni* or *Kecimen* variety based on 7 image-derived features (area, perimeter, convex area, etc.).

We apply a **hybrid approach** combining **unsupervised learning** (PCA, clustering) and **supervised learning** (logistic regression, Lasso, SVM, LDA/QDA) to:

* Understand structure in the data
* Reduce dimensionality
* Train and compare multiple classifiers

---

## ⚖️ Unsupervised Learning

### 1. Exploratory Analysis

* Summary stats, boxplots, and pairwise correlations
* Highlights strong linear relationships among size-related variables (e.g., Area ≈ ConvexArea)

### 2. PCA (Principal Component Analysis)

* PCA performed both manually and using `FactoMineR`
* PC1 explains \~88% of the variance and aligns with grape size
* PC2 captures shape-related variance
* Biplots and 2D projections show good visual separation between classes

### 3. Clustering

* Hierarchical clustering (Ward method) using original features and PCA scores
* Validated with classification error, ARI (Adjusted Rand Index), and silhouette analysis
* Clustering aligns well with known labels when using PC space

---

## 🎓 Supervised Learning

### 1. Logistic Regression

* Models trained with raw variables and PCs
* Best tradeoff found using **PC1 + PC2** as predictors
* Stepwise AIC and Lasso (`glmnet`) further refine feature selection

### 2. SVM (Support Vector Machine)

* Trained with linear and polynomial kernels using `e1071::svm`
* Polynomial kernel underperforms compared to linear
* Grid search used to tune hyperparameters

### 3. ROC and AUC

* All models evaluated with ROC curves and AUC on test set
* Computed with `pROC` and `ROCR`
* Best models (logit-2PC and SVM-linear) reach **AUC \~0.932**

### 4. Discriminant Analysis

* LDA and QDA implemented both manually and using `MASS`
* LDA directions consistent with PCA structure

---

## 🏆 Performance Summary

| Model            | AUC   | Test Error |
| ---------------- | ----- | ---------- |
| Logistic (2 PCs) | 0.932 | 12.8%      |
| SVM (Linear)     | 0.931 | 12.1%      |
| Lasso            | 0.923 | 13.2%      |
| Stepwise AIC     | 0.925 | 12.5%      |
| LDA              | 0.918 | 13.0%      |
| QDA              | 0.912 | 13.7%      |
| SVM (Poly)       | 0.911 | 15.2%      |

---

## 📊 Key Insights

* PCA greatly enhances performance and interpretability
* The **logistic regression on PC1 + PC2** achieves strong results with high transparency
* **SVM with linear kernel** slightly outperforms in error but at higher model complexity
* Discriminant methods (LDA, QDA) show consistent trends with the structure revealed by PCA

---

## 📁 Files

| File          | Description                                                         |
| ------------- | ------------------------------------------------------------------- |
| `exo1.R`      | Unsupervised learning: EDA, PCA, clustering                         |
| `exo2.R`      | Supervised learning: logistic regression, Lasso, SVM, ROC, LDA, QDA |
| `grape.R`     | Combined script with modular functions and full visualizations      |
| `Raisin.xlsx` | Source dataset (Cinar et al., 2020)                                 |

---

## ✏️ Authors

* Project completed by **Gabriel Dupuis** and **Hussein Rammal** — MSc in Maths & CS at Institut Polytechnique de Paris - ENSTA Paris

---

## 📚 Reference

> Cinar, I., Koklu, M., & Tasdemir, S. (2020). Classification of Raisin Varieties Using Machine Learning Algorithms. *Measurement*, 168, 108419.

---

## 🚀 Manual vs. Library Code

* PCA, logistic regression, clustering, and LDA/QDA all implemented both **from scratch** and using R packages
* Ensures transparency, reproducibility, and deeper understanding of model mechanics

---

*Let the data ferment into insights!* 🥃
