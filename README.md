<img width="1254" height="1254" alt="image" src="https://github.com/user-attachments/assets/1da54564-1b3f-4dfa-978e-78d598527fac" />

# TLS Finder

## Project Information

| Item | Details |
|------|---------|
| **Published** | December 2024 (bioRxiv) |
| **DOI** | https://doi.org/10.1101/2024.12.26.630405 |
| **First Public Release** | December 26, 2024 |
| **Purpose** | Automated detection of tertiary lymphoid structures from spatial immune-cell coordinates |
| **Status** | Stable |

**TLS Finder is now part of the MOSAIC ecosystem for graph-native spatial biology.

# TLS Finder

**Detection of tertiary lymphoid structures from spatial single-cell data**

---

## 🔬 Overview

TLS Finder is a Python-based computational method for identifying **tertiary lymphoid structures (TLS)** from spatial imaging datasets such as multiplex immunofluorescence (mIF), CODEX, and immunohistochemistry.

The method combines:

* efficient spatial indexing
* graph-based modeling
* biologically informed clustering

to detect organized immune aggregates relevant to tumor immunology and clinical outcomes.

---

## 🚀 Key Features

* TLS detection from spatial coordinates
* Graph-based modeling using NetworkX
* Fast spatial queries via KDTree
* Compatible with mIF, CODEX, and spatial transcriptomics
* Scalable to large datasets

---

## ⚙️ Method Overview

TLS Finder identifies structured immune aggregates through the following steps:

1. **KDTree Construction**
   Efficient spatial indexing of all cells for rapid neighbor lookup.

2. **Graph Construction**
   Cells are represented as nodes in a graph with phenotype and spatial attributes.

3. **Edge Formation**
   Edges are created between cells within a defined radius, forming spatial neighborhoods.

4. **Cluster Identification**
   Connected components are evaluated based on immune cell composition.

5. **TLS Classification**
   Clusters meeting biologically informed thresholds (e.g., B-cell and T-cell density) are classified as:

   * Lymphoid Aggregates (LA)
   * (extendable to more refined TLS categories)

---

## 📂 Input Requirements

A pandas DataFrame with:

| Column    | Description                 |
| --------- | --------------------------- |
| X         | Cell X-coordinate           |
| Y         | Cell Y-coordinate           |
| Phenotype | Cell type (e.g., CD20, CD8) |

---

## 📊 Output

* DataFrame containing:

  * TLS / aggregate centers
  * coordinates
  * classification labels

* Visualization-ready outputs (for plotting TLS locations)

---

## 🧪 Core Function

### `measure_TLS(df, radius)`

#### Parameters

* **df**: pandas DataFrame
  Contains X, Y, and Phenotype columns

* **radius**: float
  Distance threshold for defining cellular neighborhoods

#### Returns

* DataFrame with:

  * cluster indices
  * TLS center coordinates
  * classification labels

---

## 💻 Example Usage

```python
import pandas as pd
from tls_finder import measure_TLS

data = {
    'X': [1, 2, 3, 4, 5],
    'Y': [5, 4, 3, 2, 1],
    'Phenotype': ['CD20', 'CD8', 'CD20', 'CD8', 'CD20']
}

df = pd.DataFrame(data)

result = measure_TLS(df, radius=50)
print(result)
```

---

## 📦 Dependencies

* pandas
* numpy
* scipy
* networkx

---

## 🖥️ Applications

* **Local analysis**: run on single datasets for TLS detection
* **Web app integration**: interactive parameter selection and visualization
* **Pipeline integration**: plug into spatial omics workflows

---

## 📖 Citation

If you use TLS Finder in your research, please cite:

Koksoy, A.A. (2024).
*TLS_Finder: An algorithm for identifying tertiary lymphoid structures using immune cell spatial coordinates.*
bioRxiv
https://doi.org/10.1101/2024.12.26.630405

---

## 🤝 Contributing

Contributions and improvements are welcome.
Please open an issue or submit a pull request.

---

## 📜 License

© 2026 Ayse A. Koksoy

This project is licensed under the Apache License 2.0.
See the LICENSE file for details.

---

## 👩‍🔬 Author

**Ayse A. Koksoy, MD, PhD**
Computational scientist in spatial omics and multiplex imaging

"I build AI-assisted pipelines that integrate spatial, transcriptomic, and imaging data into interpretable biological insights."

---
