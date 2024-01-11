# SMURF
A matrix factorization method for single-cell

## Pre-requirements
* python3
* numpy
* pandas
* scipy
* scikit-learn
* umap-learn

## Installation

### Installation with pip
To install with pip, run the following from a terminal:
```Bash
pip install smurf-imputation
```

## Usage

### Basic use
```Python
import smurf
import pandas as pd

# read your data, the rows in the data represent genes, and the columns represent cells
# demo file https://github.com/deepomicslab/SMURF/blob/master/data.csv
data = pd.read_csv("data.csv", header=0, index_col=0)

# create a SMURF object which only return the imputed data
operator = smurf.SMURF(n_features=10, estimate_only=True)

# impute
data_imputed = operator.smurf_impute(data)

# create a SMURF object
operator = smurf.SMURF(n_features=10, estimate_only=False)

# impute
res = operator.smurf_impute(data)

# get the results
data_imputed = res["estimate"]

gene_matrix = res["gene latent factor matrix"]

cell_matrix = res["cell latent factor matrix"]


# get cell-circle
cell_circle = operator.smurf_cell_circle(
  n_neighbors=20, min_dist=0.01, major_axis=3, minor_axis=2, k=0.2
)

# or get cell-cirecle directly from you own data
mapper = smurf.SMURF()
cell_circle = mapper.smurf_cell_circle(cells_data=your_own_data)

# get result in different coordinate
angle = cell_circle["angle"]
plane_embedding = cell_circle["plane_embedding"]



```

### Parameters
```Python
SMURF(n_features=20, steps=10, alpha=1e-5, eps=10,lambda2=0.1, noise_model="Fano", normalize=True, estimate_only=False)
```
Parameters

* n_features : int, optional, default: 20

    The number of features during the matrix factorizaiton.

* steps : int, optional, default: 10

    The max number of iteration.

* alpha : float, optional, default: 1e-5

    gradient update step size. It can be so different with different dataset, please try more for a better result.
  
* eps : float, optional, default: 10
    
    The threshold at which the objective function stops updating

* lambda2 : float, optional, default: 0.1
    
    The coefficient of L2 regularization
  
* noise_model: boolean, optional, default: "Fano"
    
    Our hypothetical noise model. We offer three options:
    * CV : constant variance
    * Fano : Fano factor
    * CCV : constant coefficient of variation
    
    We found that generally the fano model is the most stable.
    
* normalize : boolean, optional, default: True

    By default, SMURF takes in an unnormalized matrix and performs library size normalization during the denoising step. However, if your data is already normalized or normalization is not desired, you can set normalize=False.

* estimate_only : boolean, optional, default: False

    Generally, the SMURF returns a dictionary which contains the imputed matrix and gene latent factor matrix and cell latent factor matrix. If you have no need of the latent factor matrix, you can set estimate_only=True.

```Python
smurf_cell_circle(cells_data=None, n_neighbors=20, min_dist=0.01, major_axis=3, minor_axis=2, k=0.2)
```
* cells_data : array of 2D, optional, default: None
  
    Cells data to be processed. If it's not None, the model will process your own data, or please use SMURF process the original data and the model will calculate the cell circle from the cell latent factor matrix of the feedback.
  
* n_neighbors : int, optional, default: 20
  
    The parameter controls how our model balances local versus global structure in the data.
  
* min_dist : float, optional, default: 0.01

    This parameter controls how tightly SMURF is allowed to pack points together

* major_axis : float, optional, default: 3

    Major axis length of the oval.

* minor_axis : float, optional, default: 2
    
    Minor axis length of the oval.

* k : float, optional, default: 0.2

    Deformation parameter of the oval.




