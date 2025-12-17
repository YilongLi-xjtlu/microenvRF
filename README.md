README
================

# microenvRF

`microenvRF` is an R package for **tumor microenvironment (TME)**
classification using a Random Forest model. It provides:

- **Cell type classification**: Cancer / T_Cell / Fibroblast  
- **Disease status classification**: Tumor vs Healthy_Control

> Note (coursework requirement): the exported prediction function
> returns **class labels** (not probabilities).

## 1. Installation (GitHub)

``` r
# install.packages("remotes")
remotes::install_github("YilongLi-xjtlu/microenvRF")
```

## 2. Usage demo

### 2.1 After installation

``` r
library(microenvRF)

example_path <- system.file("extdata", "example_input.csv", package = "microenvRF")
newdata <- read.csv(example_path)

predict_microenv(newdata, task = "cell_type")
predict_microenv(newdata, task = "disease_status")
```

### 2.2 Development mode (for knitting this repository)

``` r
devtools::load_all()

example_path <- system.file("extdata", "example_input.csv", package = "microenvRF")
if (example_path == "") example_path <- file.path("inst", "extdata", "example_input.csv")
newdata <- read.csv(example_path)

predict_microenv(newdata, task = "cell_type")
predict_microenv(newdata, task = "disease_status")
```

## 3. Model performance (Performance)

The following figures are taken from the project report and summarize
model performance:

### 3.1 Cell type (Cancer / T_Cell / Fibroblast)

![](man/figures/roc_pr_celltype.png)

### 3.2 Disease status (Tumor vs Healthy_Control)

![](man/figures/roc_pr_disease.png)
