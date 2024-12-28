# NCBI Genotype-Phenotype Analysis

A Python-based tool for analyzing and visualizing relationships between phenotypes and genotypes involving the brain using NCBI's E-utilities API.

## Project Overview

This project aims to:
- Fetch gene expression data from NCBI databases
- Analyze relationships between brain-related phenotypes and genotypes
- Create visualizations of gene expression patterns
- Provide statistical analysis of gene-phenotype associations

## Installation

1. Clone the repository:
```bash
git clone https://github.com/mycoleb/NCBIGenotypePhenotype.git
cd NCBIGenotypePhenotype
```

2. Set up Python environment:
```bash
# Using pip
python -m venv brain-genetics-env
source brain-genetics-env/Scripts/activate  # On Windows
# OR
./brain-genetics-env/Scripts/activate  # On Unix/MacOS

# Install dependencies
pip install -r requirements.txt
```

## Project Structure

```
NCBIGenotypePhenotype/
├── data/               # Data storage (not tracked by git)
├── notebooks/         # Jupyter notebooks for analysis
├── src/              # Source code
│   └── gene_analysis.py
├── requirements.txt   # Project dependencies
└── README.md         # Project documentation
```

## Dependencies

- Python 3.10+
- pandas
- numpy
- matplotlib
- seaborn
- requests
- jupyter

## Usage

1. Start with the Jupyter notebooks for data exploration:
```bash
jupyter notebook notebooks/
```

2. Use the gene analysis module:
```python
from src.gene_analysis import GeneExpressionAnalyzer

analyzer = GeneExpressionAnalyzer()
genes = ['BDNF', 'DRD2', 'COMT']
data = analyzer.fetch_expression_data(genes)
analyzer.create_visualization(data)
```

## Features

- NCBI E-utilities API integration
- Gene expression data retrieval
- Statistical analysis tools
- Visualization capabilities
- Phenotype-genotype relationship analysis

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- NCBI for providing the E-utilities API
- Contributors and maintainers
