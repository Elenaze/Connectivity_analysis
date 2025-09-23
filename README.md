
# Connectivity analysis for Molecular Graph Construction

A computational tool for determining optimal connectivity parameters (α-values) in molecular graph construction from atomic coordinates. This analysis establishes the foundation for molecular graph encoding by defining precise connectivity criteria based on pair distribution functions (PDFs).

## Purpose & Context

Molecular graphs are fundamental representations in computational chemistry, where atoms are nodes and bonds are edges. However, determining when two atoms should be considered "connected" from coordinate data alone is non-trivial. Traditional approaches using fixed cutoff distances fail for flexible molecules or diverse chemical environments.

This tool solves this problem by:
- Analyzing the statistical distribution of atom-pair distances in your dataset
- Determining optimal scaling factors (α-values) that account for bond flexibility
- Providing dataset-specific connectivity criteria rather than rigid cutoffs

**Applications:**
- Molecular dynamics trajectory analysis
- Conformational sampling studies  
- Chemical reaction network construction

## Overview

This tool was developed to assess molecular connectivity for any type of dataset. The connectivity criterion is defined as:

$$
   \text{if} \quad \lvert \vec{r}_i - \vec{r}_j \rvert \le \alpha \cdot (R_i + R_j) \to \quad \text{i and j atoms are connected}
$$

where $R_i$ and $R_j$ are the atomic radii of atom i and j respectively and $\alpha$ is an optimized parameter to determine connectivity which also accounts for bond flexibility.


## Features

- **XYZ File Processing**: Load and analyze molecular structure data
- **PDF Analysis**: Calculate pair distribution functions for all atom pairs invariant by permutation.
- **Parameter Determination**:
  - Atomic radii from PDF first local maximum, iteration and comparison with predefined atomic radii from ASE library.
  - Alpha values extracted from PDF first local minimum 
- **Visualization**:
  - PDF plots with marked maxima/minima

## Installation

```bash
# Clone the repository
git clone https://github.com/Elenaze/Connectivity_analysis
cd Connectivity_analysis 

# Install dependencies
pip install -r requirements.txt
```

## Usage

### Basic Usage
```bash
python main.py
```

This will:
1. Automatically find and load all `.xyz` files from the `data/` directory (supports both single-structure XYZ files and multi-structure files with multiple geometries vertically stacked together)
2. Calculate pair distribution functions for all atom types
3. Determine optimal atomic radii and α-values
4. Generate visualization plots
5. Save results to `results/` directory

### Sample Output
```
Loading systems...
Finding unique atoms...
Calculating pair distances...
Computing KDE...

Radius for C: 0.748 Å (ASE: 0.760 Å)
Radius for O: 0.474 Å
Pair: ('C', 'C'), Alpha: 1.256
Pair: ('N', 'C'), Alpha: 1.329

Results saved to: /path/to/results
```

## Output & Results

### Generated Files
```
results/
├── analysis_results.json    # Summary statistics and final α-values
├── kde_[ATOM]_[ATOM].csv   # Kernel density estimation data for each atom pair
├── full_results.pkl        # Complete analysis data (Python pickle format)
└── plots/
    └── pdf_[ATOM]_[ATOM].pdf  # Distribution plots with marked peaks/valleys
```

### Interpreting Results

**Key Outputs:**
- **Atomic Radii**: Optimized radii for each element type in your dataset
- **Alpha Values**: Scaling factors for each atom pair type (typically 1.0-3.0)
- **Distribution Statistics**: Count, range, and shape of distance distributions

**Example Alpha Interpretation:**
- α = 1.0: Standard covalent bonding (rigid bonds)
- α = 1.5: Flexible bonding allowing for conformational variation  
- α > 2.0: Very flexible connections or weak interactions

**Quality Indicators:**
- Well-defined PDF peaks indicate good structural data
- "Distance too large" warnings suggest outliers or very flexible regions
- α-values should be chemically reasonable (avoid α > 5.0)

## Method Details

### Atomic Radii Determination
1. Analysis of same-atom pair distributions
2. Refinement using different-atom pairs
3. Validation against ASE reference values

### Alpha Value Calculation
1. PDF maxima/minima identification
2. Alpha calculation: minimum position / sum of atomic radii
3. Physical constraints validation

## References

- ASE Covalent Radii: Cordero et al., Dalton Trans., 2008, 2832-2838 DOI:10.1039/B801115J
- Perez-Mellor, Ariel F., and Riccardo Spezia. "Determination of kinetic properties in unimolecular dissociation of complex systems from graph theory based analysis of an ensemble of reactive trajectories." The Journal of Chemical Physics 155.12 (2021).

## License

MIT License

