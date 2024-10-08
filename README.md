# agevarycoef


`agevarycoef` is an R package that models time-varying coefficients in a variety of contexts, specifically designed for age-varying coefficients. This package allows you to fit models where the relationship between covariates and outcomes changes smoothly with age.

## Features

- **Time-varying effect models (TVEM):** Allows for modeling coefficients that change as a function of age.
- **Cross-validation support:** Perform model validation using K-fold cross-validation to ensure robust results.
- **Elastic net support:** Apply elastic net regularization to improve model prediction and avoid overfitting.

Also support for fitting cubic splines, xgboost, and Gaussian processes.


## Installation

To install the development version of `agevarycoef`, run the following command:

```r
# First, install devtools if not installed
# install.packages("devtools")

# Then, install the package from GitHub
devtools::install_github("bengriffin3/agevarycoef")
