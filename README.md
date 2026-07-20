>> # Revisiting Logistic Regression for High-Dimensional Gene Expression Data

This repository contains the source code, analysis scripts, and supporting materials for the study:

> **Revisiting Logistic Regression for High-Dimensional Gene Expression Data**

The study proposes a logistic-regression-inspired linear framework for feature selection and classification in high-dimensional gene expression datasets, particularly in settings where the number of attributes is substantially larger than the number of samples.

The method was evaluated using three publicly available sepsis transcriptomic datasets:

- GSE12624
- GSE13205
- GSE69063


## Running the analyses

### 1. Prepare the datasets

Run the preprocessing script in MATLAB:

```matlab
run('main.m')
```

## Datasets

The datasets analyzed in this study are publicly available from the NCBI Gene Expression Omnibus.

| Dataset | Samples | Control/non-sepsis | Sepsis | GEO accession |
|---|---:|---:|---:|---|
| GSE12624 | 70 | 36 | 34 | GSE12624 |
| GSE13205 | 21 | 8 | 13 | GSE13205 |
| GSE69063 | 90 | 33 | 57 | GSE69063 |

No new primary biological data were generated for this study.

## Authors

- Rossana O. Souza
- Wellington Francisco Rodrigues
- Bráulio R. G. M. Couto
- Marcos A. dos Santos

Corresponding author:

**Rossana O. Souza**  
Email: rossanaoliveirasouza@gmail.com

## Ethical statement

This study used exclusively publicly available, de-identified datasets obtained from the NCBI Gene Expression Omnibus.

No new participants were recruited, no new biological samples were collected, and no identifiable personal data were accessed.

## Contact

Questions about the implementation or the study may be directed to:

**Rossana O. Souza**  
rossanaoliveirasouza@gmail.com

