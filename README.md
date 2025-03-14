# TLS-Finder
An algorithm in Python that locates the TLS and lymphoid aggregates in tissues using immune cell coordinates.
If you find this code useful please cite as__ TLS_Finder: An algorithm for Identifying Tertiary Lymphoid Structures Using Immune Cell Spatial Coordinates; Ayse A Koksoy, Maria Esther Salvatierra, Patient Mosaic Team, Luisa Maren Solis Soto, Cara Haymaker; 
### doi: https://doi.org/10.1101/2024.12.26.630405

# TLS Finder

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



### JavaScript Snippet for Tracking

```markdown
## Tracking Globe

Below is an HTML script for embedding a ClustrMaps globe to visually track visits to your site. Ensure this fits your project's privacy and functionality requirements before using.

```html
<script type="text/javascript" id="clstr_globe" src="//clustrmaps.com/globe.js?d=DmltlICmBQZ3Y


