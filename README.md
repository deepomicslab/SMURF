# SCEnd
A matrix factorization method for single-cell

## Pre-requirements
* python3
* numpy
* pandas
* scipy
* scikit-learn

## Installation

### Installation with pip
To install with pip, run the following from a terminal:
```Bash
pip install scend
```

## Usage

### Basic use
```Python
import scend
import pandas as pd

# read your data, the rows in the data represent genes, and the columns represent cells
data = pd.read_csv("data.csv", header=0, index_col=0)

# create a SCEnd object which only return the imputed data
operator = scend.SCEnd(n_features=10, estimate_only=True)

# impute
data_imputed = operator.scend_impute(data)

# create a SCEnd object
operator = scend.SCEnd(n_features=10, estimate_only=False)

# impute
res = operator.scend_impute(data)

# get the results
data_imputed = res["estimate"]

gene_matrix = res["gene latent factor matrix"]

cell_matrix = res["cell latent factor matrix"]

```

### Parameters
```Python
SCEnd(n_features=20, steps=10, alpha=1e-5, eps=10, noise_model="Fano", normalize=True, estimate_only=False)
```
Parameters

* n_features : int, optional, default: 20

    The number of features during the matrix factorizaiton.

* steps : int, optional, default: 0.5

    The max number of iteration.

* alpha : float, optional, default: 1e-5

    gradient update step size. It can be so different with different dataset, please try more for a better result.
  
* eps : float, optional, default: 10
    
    The threshold at which the objective function stops updating
  
* noise_model: boolean, optional, default: "Fano"
    
    Our hypothetical noise model. We offer three options:
    * CV : constant variance
    * Fano : Fano factor
    * CCV : constant coefficient of variation
    
    We found that generally the fano model is the most stable.
    
* normalize : boolean, optional, default: True

    By default, SCEnd takes in an unnormalized matrix and performs library size normalization during the denoising step. However, if your data is already normalized or normalization is not desired, you can set normalize=False.

* estimate_only : boolean, optional, default: False

    Generally, the SCEnd returns a dictionary which contains the imputed matrix and gene latent factor matrix and cell latent factor matrix. If you have no need of the latent factor matrix, you can set estimate_only=True.


