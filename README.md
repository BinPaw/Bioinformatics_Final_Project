# Molecular Classification of Cancer by Gene Expression

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.12-brightgreen.svg)
![scikit-learn](https://img.shields.io/badge/scikit--learn-1.x-orange.svg)

## Abstract

მოლეკულური მონაცემებიდან კიბოს ტიპებისა და ქვეტიპების იდენტიფიცირება გამოთვლითი ბიოლოგიის ფუნდამენტური პრობლემაა. 1999 წელს გოლუბმა და სხვებმა (Golub et al.) ჟურნალ Science-ში გამოაქვეყნეს საეტაპო ნაშრომი, სადაც აჩვენეს, რომ გენთა ექსპრესიის მიკრომატრიცულ მონაცემებზე მანქანური სწავლების გამოყენებით შესაძლებელი იყო მწვავე ლიმფობლასტური ლეიკემიის (ALL) და მწვავე მიელოიდური ლეიკემიის (AML) კლასიფიკაცია 85.3%-იანი სიზუსტით დამოუკიდებელ სატესტო ნაკრებზე. აღნიშნულმა ნაშრომმა, რომელიც დღეისთვის 20,000-ზე მეტჯერაა ციტირებული, საფუძველი ჩაუყარა ექსპრესიაზე დაფუძნებული კიბოს დიაგნოსტიკის პარადიგმას.
წინამდებარე კვლევაში, ლეიკემიის ორიგინალურ მონაცემებზე (72 ნიმუში, 7,129 გენი) და თანამედროვე Python/scikit-learn ინსტრუმენტების გამოყენებით, ჩვენ აღვადგენთ გოლუბის კლასიფიკაციის სქემას, ხოლო შემდეგ ვამოწმებთ, რამდენად განზოგადდება იგივე მეთოდოლოგია არსებითად უფრო რთულ ამოცანაზე: ძუძუს კიბოს ექვსი მოლეკულური ქვეტიპის კლასიფიცირებაზე დამოუკიდებელი მიკრომატრიცული მონაცემთა ნაკრებიდან (GSE45827, 151 ნიმუში, 54,675 გენი). მახასიათებელთა შესარჩევად ვიყენებთ გოლუბის „სიგნალი-ხმაურის თანაფარდობას“ (S2N) და ვადარებთ მას ANOVA F-ტესტს.
შეფასებულია რვა კლასიფიკატორი, მათ შორის: უახლოესი ცენტროიდი (Golub-ის ანალოგი), წრფივი SVM, RBF SVM, შემთხვევითი ტყე (Random Forest), KNN, ლოგისტიკური რეგრესია და გრადიენტული ბუსტინგი. ჩვენმა „უახლოესი ცენტროიდის“ კლასიფიკატორმა ლეიკემიის სატესტო ნაკრებზე 82.4% აჩვენა, რაც ახლოსაა გოლუბის 85.3%-იან შედეგთან. თანამედროვე კლასიფიკატორებმა 1999 წლის მიდგომას აჯობეს: წრფივმა SVM-მა მიაღწია 91.2%-ს, ხოლო პარამეტრებზე ოპტიმიზებულმა (grid-search) SVM-მა — 94.1%-ს. ძუძუს კიბოს შემთხვევაში, იგივე სქემით, ოპტიმიზებული SVM-ის გამოყენებით, 98.7% სიზუსტე იქნა მიღწეული. წრფივმა SVM-მა ორივე მონაცემთა ნაკრებზე ჯვარედინი ვალიდაციის (CV) 96.0% სიზუსტე აჩვენა, რაც ადასტურებს, რომ გენთა ექსპრესიაზე დაფუძნებული მანქანური სწავლების კლასიფიკაციის პარადიგმა არის როგორც რეპროდუცირებადი, ასევე განზოგადებადი სხვადასხვა ტიპის კიბოს, მიკრომატრიცული პლატფორმებისა და კლასიფიკაციის სირთულის მიმართ.

## Introduction

Cancer is fundamentally a disease of aberrant gene expression. Different cancer types and subtypes exhibit distinct transcriptomic profiles driven by particular oncogenic mutations, epigenetic alterations, and microenvironmental factors. Leukemias arise from hematopoietic precursors, with ALL originating from lymphoid progenitors and AML from myeloid progenitors. These distinct cellular origins produce markedly different gene expression signatures. Similarly, breast cancer is classified into molecular subtypes — basal-like, HER2-enriched, luminal A, luminal B, and normal-like — each with different prognosis, treatment response, and transcriptomic identity.

In 1999, Golub et al. published "Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring" in *Science* (Vol. 286, pp. 531–537). Using Affymetrix Hgu6800 microarrays measuring 7,129 gene probes across 72 bone marrow samples (47 ALL, 25 AML), they demonstrated that (1) unsupervised methods (self-organizing maps) could discover cancer classes without prior labeling, and (2) supervised methods (weighted voting with signal-to-noise gene selection) could predict class membership with 85.3% accuracy (29/34) on a held-out test set. This paper established the entire field of expression-based cancer diagnostics and remains one of the most cited works in computational biology.

Our central research question is: does Golub's gene expression-based classification paradigm generalize beyond binary leukemia to a more complex, multi-class problem? We test three linked hypotheses: (1) Golub's classification results on ALL vs. AML leukemia can be reproduced using modern Python/scikit-learn implementations, (2) the same ML pipeline accurately classifies molecular subtypes of a completely different cancer (breast cancer, 6 classes), and (3) modern classifiers (SVM, Random Forest, Gradient Boosting) outperform the 1999-era approach on both datasets. Our approach is to reproduce Golub's analysis on the original data, then apply the identical pipeline to the GSE45827 breast cancer gene expression dataset.

## Materials and Methods

### Datasets

We used two publicly available gene expression microarray datasets:

**Phase 1 — Golub Leukemia Dataset:** 72 bone marrow samples profiled on Affymetrix Hgu6800 arrays (7,129 probes). Pre-split into 38 training and 34 test samples. Class distribution: 47 ALL (27 train, 20 test) and 25 AML (11 train, 14 test). Obtained as CSV files from Kaggle (crawford/gene-expression). Data required transposition (genes as rows to samples as rows) and removal of Affymetrix detection call columns.

**Phase 2 — GSE45827 Breast Cancer Dataset:** 151 samples profiled on Affymetrix HG-U133 Plus 2.0 arrays (54,675 probes). Six molecular classes: basal-like (41), HER2-enriched (30), luminal B (30), luminal A (29), cell line (14), and normal tissue (7). Obtained through CuMiDa (Curated Microarray Database) from Kaggle (brunogrisci/breast-cancer-gene-expression-cumida). CuMiDa applied background correction, RMA normalization, and quality assessment, providing a machine-learning-ready CSV file.

### Preprocessing

Both datasets were standardized using scikit-learn's StandardScaler (zero mean, unit variance per feature). For the breast cancer data, a variance threshold filter (median variance) was applied to remove low-information probes, reducing the feature space from 54,675 to 11,715 genes. Missing values were checked and none were found in either dataset.

### Feature Selection

For the binary leukemia problem, we implemented Golub's signal-to-noise ratio (S2N): for each gene *g*, S2N(*g*) = (μ_ALL − μ_AML) / (σ_ALL + σ_AML). Genes were ranked by |S2N| and the top 50 were selected (25 ALL-correlated, 25 AML-correlated). For comparison, we also applied scikit-learn's SelectKBest with ANOVA F-test (k=50). The two methods agreed on 31/50 genes (62% overlap), confirming that both approaches capture similar discriminative information.

For the multi-class breast cancer problem, we used the ANOVA F-test exclusively (S2N is only defined for binary classification). Feature sets of 50, 100, 200, and 500 genes were tested. Grid search identified 200 as optimal.

### Model Architectures

We evaluated eight classifiers on both datasets:

1. **Nearest Centroid:** The closest analog to Golub's weighted voting classifier. Assigns samples to the class whose centroid is nearest. Serves as our reproduction baseline.
2. **Linear SVM:** Support vector machine with linear kernel. Well-suited for high-dimensional, low-sample-size settings where p >> n.
3. **RBF SVM:** SVM with radial basis function kernel for nonlinear decision boundaries.
4. **Random Forest:** Ensemble of 200 decision trees with bootstrap aggregation.
5. **KNN (k=3, k=5):** K-nearest neighbors for instance-based classification.
6. **Logistic Regression (L1/L2):** Regularized logistic regression with L1 (sparsity) and L2 (ridge) penalties.
7. **Gradient Boosting:** Sequential ensemble of 100 boosted decision trees.

Grid search with cross-validation was used to optimize hyperparameters for SVM (C, kernel, number of features) and Logistic Regression (C, penalty, PCA components).

### Evaluation Metrics

We evaluated our models using the following metrics:
1. **Accuracy:** Fraction of correctly classified samples (test set and cross-validation).
2. **Macro F1-Score:** Harmonic mean of precision and recall, averaged across classes.
3. **Confusion Matrix:** Per-class classification performance visualization.
4. **Adjusted Rand Index (ARI):** Agreement between K-Means clustering and true labels for unsupervised evaluation.
5. **Stratified 5-Fold CV:** Robust performance estimation preserving class proportions in each fold.

### Comparison with Golub et al. (1999)

We compared our results directly with Golub's original weighted voting classifier. Both approaches were evaluated on the same 38/34 train–test split using the same dataset to ensure a fair comparison. Additionally, we ran stratified 5-fold cross-validation on all 72 samples to provide more robust performance estimates.

### Implementation

All analyses were conducted in Python 3.12 on Kaggle Notebooks. We used the scikit-learn library for building and evaluating classifiers, PCA, feature selection, and grid search. Additional packages included pandas, numpy, matplotlib, seaborn, and scipy (hierarchical clustering). All experiments were conducted on Kaggle's cloud infrastructure. No specialized bioinformatics software was required.

## Results

### Phase 1: Leukemia Classification

| Classifier | Test Accuracy | Correct | CV Mean | CV Std | CV F1 |
|---|---|---|---|---|---|
| Golub (1999) Original | 85.3% | 29/34 | N/A | N/A | N/A |
| Nearest Centroid | 82.4% | 28/34 | 0.834 | 0.176 | 0.800 |
| **Linear SVM** | **91.2%** | **31/34** | **0.960** | 0.080 | 0.910 |
| **Gradient Boosting** | **91.2%** | **31/34** | 0.846 | 0.056 | 0.910 |
| Logistic Reg. (L1) | 88.2% | 30/34 | 0.863 | 0.146 | 0.880 |
| Logistic Reg. (L2) | 88.2% | 30/34 | 0.905 | 0.123 | 0.880 |
| Random Forest | 73.5% | 25/34 | 0.959 | 0.054 | 0.670 |
| KNN (k=5) | 70.6% | 24/34 | 0.820 | 0.090 | 0.650 |
| RBF SVM | 61.8% | 21/34 | 0.860 | 0.079 | 0.440 |
| **Grid Search SVM** | **94.1%** | **32/34** | — | — | 0.940 |

Grid Search SVM best parameters: C=1, RBF kernel, top 25 S2N genes.

Our Nearest Centroid classifier (82.4%) closely reproduced Golub's original result (85.3%), with the small discrepancy expected given the tiny test set (n=34) where a single misclassification changes accuracy by ~3%. Modern classifiers clearly outperformed the 1999 approach: Linear SVM and Gradient Boosting both reached 91.2% (+5.9% over Golub), and grid-search SVM achieved 94.1%.

### Phase 2: Breast Cancer Classification

| Classifier | CV Accuracy | CV Std | Macro F1 |
|---|---|---|---|
| **Linear SVM** | **96.0%** | ±0.024 | 0.960 |
| Random Forest | 95.4% | ±0.040 | 0.954 |
| RBF SVM | 94.7% | ±0.040 | 0.947 |
| Nearest Centroid | 92.7% | ±0.039 | 0.927 |
| Logistic Regression | 92.7% | ±0.049 | 0.927 |
| KNN (k=5) | 86.7% | ±0.043 | 0.867 |
| Gradient Boosting | 83.4% | ±0.099 | 0.834 |
| **Grid Search SVM** | **98.7%** | — | 0.990 |

Grid Search SVM best parameters: C=0.1, linear kernel, top 200 ANOVA genes.

The optimized SVM classified all six breast cancer subtypes with near-perfect accuracy. Per-class F1 scores ranged from 0.966 (HER2) to 1.000 (luminal A, cell line, normal). HER2 was the only subtype with any misclassifications (2 samples).

### Cross-Dataset Comparison

| Metric | Leukemia | Breast Cancer |
|---|---|---|
| Classes | 2 | 6 |
| Samples | 72 | 151 |
| Genes | 7,129 | 54,675 |
| PCA (PC1+PC2) | 27.0% | 28.1% |
| K-Means ARI | 0.260 | 0.726 |
| Best CV Accuracy | 96.0% (Linear SVM) | 96.0% (Linear SVM) |
| Grid Search Best | 94.1% | 98.7% |

Linear SVM achieved identical 96.0% CV accuracy on both datasets, despite the breast cancer problem having 3× more classes and 8× more features. This demonstrates remarkable consistency of the expression-based classification paradigm.

## Conclusion

This study demonstrates that the gene expression-based ML classification paradigm established by Golub et al. in 1999 is both reproducible and generalizable. Our Nearest Centroid classifier (82.4%, 28/34) closely reproduced Golub's original result (85.3%, 29/34), with the small discrepancy expected given the tiny test set where a single misclassification changes accuracy by ~3%.

Modern classifiers clearly outperformed the 1999 approach. Linear SVM and Gradient Boosting both achieved 91.2% (31/34) on the original test split, improving on Golub by +5.9 percentage points. Grid search SVM reached 94.1%.

Most importantly, the identical pipeline generalized to breast cancer: the optimized SVM classified six breast cancer subtypes at 98.7% accuracy — despite 3× more classes and 8× more features. Three subtypes (luminal A, cell line, normal) achieved perfect classification (F1 = 1.000). HER2 was the only subtype with any misclassifications (2 samples, F1 = 0.966).

Linear SVM was the most consistent classifier, achieving 96.0% CV accuracy on both datasets regardless of cancer type, platform, or classification complexity. Feature selection proved critical: reducing 54,675 genes to 200 via ANOVA improved both accuracy and computational efficiency (270× reduction).

## Individual Contributions

#### Badri — Data Lead
Loaded and cleaned both datasets (Golub CSV transposition, call column removal, CuMiDa preprocessing). Performed exploratory data analysis including class distributions, expression statistics, and variance analysis. Computed summary statistics and verified data integrity for both phases.

#### Dimitri — Feature Selection Specialist
Implemented Golub’s Signal-to-Noise Ratio from scratch in Python. Applied ANOVA F-test (SelectKBest) for multi-class feature selection on GSE45827. Compared S2N vs. ANOVA overlap (31/50 genes, 62%). Created feature heatmaps and gene expression boxplots for both datasets.

#### Giorgi — Classifier Group 1
Trained and evaluated Nearest Centroid, Linear SVM, RBF SVM, and Logistic Regression (L1/L2) on both datasets. Performed Grid Search optimization for SVM and Logistic Regression pipelines. Generated confusion matrices and per-class classification reports.

#### Lasha — Classifier Group 2
Trained and evaluated Random Forest, KNN (k=3, k=5), and Gradient Boosting on both datasets. Ran stratified 5-fold cross-validation across all classifiers. Extracted Random Forest feature importances and analyzed learning curves.

#### Anano — Visualization & Synthesis
Performed all PCA analyses (2D and 3D) and unsupervised clustering (hierarchical, K-Means) for both datasets. Created all publication-quality figures across three phases (20+ figures). Designed the presentation, compiled the comparative analysis (Phase 3), and led report assembly.

## Future Work

1. Application to a third cancer type (e.g., lung cancer) to further validate generalizability
2. Use of RNA-seq data from TCGA for modern platform comparison
3. Investigation of deep learning classifiers (neural networks, autoencoders)
4. Pathway enrichment analysis on top discriminating genes to understand biological mechanisms
5. Development of an interactive web application (Streamlit) for real-time classification

## References

1. Golub, T. R., Slonim, D. K., Tamayo, P., et al. (1999). Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring. *Science*, 286(5439), 531–537.

2. Gruosso, T., Mieulet, V., Carber, M., et al. (2016). Chronic oxidative stress promotes H-Ras-induced superinvasion. *PNAS*. GEO Accession: GSE45827.

3. Sørlie, T., Perou, C. M., Tibshirani, R., et al. (2001). Gene expression patterns of breast carcinomas distinguish tumor subclasses with clinical implications. *PNAS*, 98(19), 10869–10874.

4. Parker, J. S., Mullins, M., Cheang, M. C. U., et al. (2009). Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes. *J Clin Oncol*, 27(8), 1160–1167.

5. Pedregosa, F., Varoquaux, G., Gramfort, A., et al. (2011). Scikit-learn: Machine Learning in Python. *JMLR*, 12, 2825–2830.

6. Grisci, B. I., Krause, M. J., & Dorn, M. (2019). CuMiDa: An Extensively Curated Microarray Database for Benchmarking and Testing Machine Learning Approaches. *J Comput Biol*, 26(4), 376–386.

## Contact

For questions or collaborations, please open an issue in this repository.
