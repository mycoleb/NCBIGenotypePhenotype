import requests
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict
import time
import json

class GeneExpressionAnalyzer:
    def __init__(self):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        # Add delay between requests to comply with NCBI's guidelines
        self.request_delay = 0.34  # seconds
    
    def fetch_gene_data(self, gene_id: str) -> Dict:
        """
        Fetch gene data from NCBI E-utilities
        
        Args:
            gene_id (str): Gene symbol or ID
            
        Returns:
            dict: JSON response from the API
        """
        esearch_url = (
            f"{self.base_url}esearch.fcgi?"
            f"db=gene&term={gene_id}[Gene Name]+AND+homo+sapiens[Organism]&format=json"
        )
        
        try:
            response = requests.get(esearch_url)
            response.raise_for_status()
            time.sleep(self.request_delay)  # Be nice to the API
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data for gene {gene_id}: {e}")
            return None

    def fetch_expression_data(self, gene_ids: List[str]) -> pd.DataFrame:
        """
        Fetch expression data for multiple genes
        
        Args:
            gene_ids (List[str]): List of gene symbols
            
        Returns:
            pd.DataFrame: DataFrame containing gene expression data
        """
        expression_data = []
        
        for gene in gene_ids:
            data = self.fetch_gene_data(gene)
            if data and 'esearchresult' in data:
                # In a real implementation, we would parse the expression data
                # Here we're using dummy data for illustration
                expression_data.append({
                    'gene': gene,
                    'expression_level': float(hash(gene) % 100) / 10,  # dummy data
                    'phenotype': 'brain_development' if hash(gene) % 2 else 'neurotransmission'
                })
        
        return pd.DataFrame(expression_data)

    def create_visualization(self, data: pd.DataFrame, output_path: str = None):
        """
        Create visualizations for gene expression data
        
        Args:
            data (pd.DataFrame): DataFrame containing gene expression data
            output_path (str, optional): Path to save the visualization
        """
        # Set the style
        plt.style.use('seaborn')
        
        # Create figure with multiple subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Expression levels bar plot
        sns.barplot(data=data, x='gene', y='expression_level', hue='phenotype', ax=ax1)
        ax1.set_title('Gene Expression Levels by Phenotype')
        ax1.set_xlabel('Gene')
        ax1.set_ylabel('Expression Level (TPM)')
        
        # Phenotype distribution
        phenotype_counts = data['phenotype'].value_counts()
        sns.pieplot(phenotype_counts, ax=ax2)
        ax2.set_title('Distribution of Phenotypes')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path)
        plt.show()

def main():
    # Initialize the analyzer
    analyzer = GeneExpressionAnalyzer()
    
    # Define genes of interest
    genes_to_study = [
        'BDNF',    # Brain-derived neurotrophic factor
        'DRD2',    # Dopamine receptor D2
        'COMT',    # Catechol-O-methyltransferase
        'HTR2A',   # Serotonin receptor 2A
        'CACNA1C'  # Calcium voltage-gated channel
    ]
    
    # Fetch and analyze data
    expression_data = analyzer.fetch_expression_data(genes_to_study)
    
    # Create and save visualization
    analyzer.create_visualization(
        expression_data,
        output_path='gene_expression_analysis.png'
    )
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(expression_data.groupby('phenotype')['expression_level'].agg(['mean', 'std']))

if __name__ == "__main__":
    main()