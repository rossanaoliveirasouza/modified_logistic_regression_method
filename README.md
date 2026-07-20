>> # Revisiting Logistic Regression for High-Dimensional Gene Expression Data

This repository contains the source code, analysis scripts, and supporting materials for the study:

> **Revisiting Logistic Regression for High-Dimensional Gene Expression Data**

The study proposes a logistic-regression-inspired linear framework for feature selection and classification in high-dimensional gene expression datasets, particularly in settings where the number of attributes is substantially larger than the number of samples.

The method was evaluated using three publicly available sepsis transcriptomic datasets:

- GSE12624
- GSE13205
- GSE69063

## Overview

High-dimensional gene expression datasets commonly present the \(p \gg n\) problem, in which thousands of molecular attributes are available for a relatively small number of samples.

The proposed framework:

1. transforms binary class labels using a log-odds representation;
2. estimates feature coefficients through a regularized linear system;
3. ranks genes according to the absolute magnitude of their coefficients;
4. selects a compact subset of relevant genes;
5. evaluates classification performance using stratified cross-validation.

The proposed method was compared with:

- LASSO logistic regression;
- elastic net logistic regression.

The evaluation considered:

- area under the ROC curve;
- accuracy;
- sensitivity;
- specificity;
- F1-score;
- feature-selection stability;
- number of selected genes;
- computational runtime.

## Mathematical formulation

Let \(A \in \mathbb{R}^{N \times p}\) represent the gene expression matrix, where \(N\) is the number of samples and \(p\) is the number of attributes.

Binary class labels are transformed into a log-odds response vector \(b\). The coefficient vector \(\alpha\) is estimated by minimizing:

\[
\|A\alpha-b\|^2+\|\alpha\|^2.
\]

The corresponding normal equation is:

\[
(I+A^TA)\alpha=A^Tb.
\]

To avoid explicitly forming \(A^TA\), the problem is represented by the augmented linear system:

\[
\begin{pmatrix}
I_N & A \\
-A^T & I_p
\end{pmatrix}
\begin{pmatrix}
r \\
\alpha
\end{pmatrix}
=
\begin{pmatrix}
b \\
0
\end{pmatrix}.
\]

Features are ranked according to:

\[
|\alpha_i|.
\]

The attributes with the largest absolute coefficient values are retained for classification.

## Datasets

The datasets analyzed in this study are publicly available from the NCBI Gene Expression Omnibus.

| Dataset | Samples | Control/non-sepsis | Sepsis | GEO accession |
|---|---:|---:|---:|---|
| GSE12624 | 70 | 36 | 34 | GSE12624 |
| GSE13205 | 21 | 8 | 13 | GSE13205 |
| GSE69063 | 90 | 33 | 57 | GSE69063 |

No new primary biological data were generated for this study.

