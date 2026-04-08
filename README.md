# TLS Finder

**An algorithm for identifying tertiary lymphoid structures using spatial immune cell organization**

---

## 🔬 Overview

TLS Finder is a computational framework for detecting **tertiary lymphoid structures (TLS)** from spatial single-cell or multiplex imaging data.

It is designed to bridge **biological interpretability** and **computational rigor**, enabling:

* Spatial identification of immune niches
* Quantification of immune cell organization
* Hypothesis generation for tumor-immune interactions

---

## 🚀 Key Features

* Detection of TLS from spatial coordinates
* Graph-based spatial analysis (neighborhood modeling)
* Flexible input (CODEX, mIF, spatial transcriptomics)
* Scalable to large datasets
* Designed for translational and clinical research

---

## 🧠 Method Summary

TLS Finder models spatial relationships between immune cells using distance-based neighborhood graphs and clustering logic to identify structured immune aggregates consistent with TLS biology.

---

## 📂 Input Requirements

* Single-cell spatial coordinates (X, Y)
* Cell type annotations (e.g., CD3+, B cells, etc.)
* Optional: region labels (tumor vs stroma)

---

## 📊 Output

* TLS regions / clusters
* Spatial metrics
* Visualization-ready outputs

---

## 📖 Citation

If you use TLS Finder in your research, please cite:

Koksoy, A.A. (2024).
*TLS_Finder: An algorithm for identifying tertiary lymphoid structures using immune cell spatial coordinates.*
bioRxiv
https://doi.org/10.1101/2024.12.26.630405

---

## 🤝 Contributing

Contributions, improvements, and discussions are welcome.
Please open an issue or submit a pull request.

---

## 📜 License

This project is licensed under the Apache License 2.0.

---

## 👩‍🔬 Author

**Ayse A. Koksoy, MD, PhD**
Computational scientist in spatial omics and multiplex imaging

"I build AI-assisted pipelines that integrate spatial, transcriptomic, and imaging data into interpretable biological insights."

---

## 💡 Notes

This tool is intended for research use and aims to support reproducible and interpretable spatial analysis workflows.

---

# Specifics of TLS Finder

## Description
The TLS Finder is an algorithm written in Python, designed to classify clusters of cells into Tertiary Lymphoid Structures (TLS) categories using spatial data from immunofluorescence or immunohistochemistry datasets. The code uses a combination of KDTree for efficient spatial querying and NetworkX for graph-based cluster analysis.

## Dependencies
- pandas
- numpy
- scipy
- networkx

## Usage
Import the script and use the `measure_TLS` function by passing a DataFrame that contains the cell coordinates and phenotypes. The DataFrame should have the columns `X`, `Y`, and `Phenotype`.

## Function Details
### `measure_TLS(df, radius)`
This function takes in a DataFrame and a radius for neighborhood analysis to classify cellular aggregates that resemble TLS based on their spatial proximity and phenotype.

#### Parameters:
- `df` : pandas DataFrame containing the columns `X`, `Y`, and `Phenotype`.
  - `X` and `Y` are coordinates of the cells.
  - `Phenotype` indicates the type of cell, important for identifying relevant clusters.
- `radius` : float. The radius within which to consider other cells as neighbors.

#### Returns:
- A DataFrame containing the indices, coordinates of the centers of identified TLS, and their classifications.

### Process Overview:
1. **KDTree Construction**: A KDTree is built for the entire set of cells to facilitate rapid spatial queries.
2. **Graph Construction**: A NetworkX graph is initialized, and nodes are added for each cell, with attributes for phenotype and coordinates.
3. **Edge Creation**: Edges are added between nodes that are within the specified radius of each other, creating a connected graph based on spatial proximity.
4. **TLS Classification**:
   - For each B cell (`CD20`) or T cell (`CD8`) that hasn't been processed, the function searches for neighboring cells within the given radius.
   - It then checks the number of neighboring B and T cells. If the count meets a specific threshold (e.g., 50), the cluster is classified as a "Lymphoid Aggregate".
   - Nodes in a classified cluster are marked to prevent reclassification in subsequent iterations.
5. **Result Compilation**: The function compiles the results into a new DataFrame, listing the indices of the aggregate centers, their coordinates, and the classification.
6. **Local App**: The function runs as an app on a local computer that asks for 2 variable entries; and produces the image with LA locations.
7. **Web App**: The function runs as a webb app,that asks for 2 variable entries; and produces the image with LA locations.

## Example

```python
import pandas as pd

# Example data
data = {
    'X': [1, 2, 3, 4, 5],
    'Y': [5, 4, 3, 2, 1],
    'Phenotype': ['CD20', 'CD8', 'CD20', 'CD8', 'CD20']
}

df = pd.DataFrame(data)
```

