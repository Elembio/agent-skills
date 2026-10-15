# Code Generation — Library Reference

Available libraries for AVITI24 analysis scripts, organized by category. Import only what the task requires. See [SKILL.md](./SKILL.md) for annotation conventions and [IMPLEMENTATION.md](./IMPLEMENTATION.md) for usage patterns.

---

## Data Handling

```python
import pandas as pd           # DataFrames
import numpy as np            # Numerical operations
import scipy                  # Scientific computing
import scipy.sparse as sp     # Sparse matrices
```

## Single-Cell Analysis

```python
import anndata as ad          # AnnData objects
import scanpy as sc           # Single-cell analysis
import squidpy as sq          # Spatial analysis
import spatialdata as sd      # SpatialData objects
```

## Element Biosciences Specific

Use the `pp` / `qc` / `tl` / `vis` / `load` namespaces.

```python
from elembio_spatialdata_analysis.load import Loader
import elembio_spatialdata_analysis as ebsd

```

## Visualization

```python
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px    # Interactive plots
```

## Machine Learning

```python
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
import umap                    # UMAP dimensionality reduction
```

## Statistical Analysis

```python
from scipy import stats
import statsmodels.api as sm
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
```

## File Formats

```python
import zarr                    # Zarr storage
import h5py                    # HDF5 files
import json                    # JSON I/O
```
