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
    def create_function_network(self, df: pd.DataFrame):
        """Create a network visualization showing gene functions and their relationships"""
        G = nx.Graph()
        
        # Create nodes and edges
        for _, row in df.iterrows():
            gene_symbol = row['symbol']
            G.add_node(gene_symbol, node_type='gene')
            
            for function, confidence in row['molecular_functions']:
                if not G.has_node(function):
                    G.add_node(function, node_type='molecular')
                G.add_edge(gene_symbol, function, weight=confidence)
            
            for process, confidence in row['biological_processes']:
                if not G.has_node(process):
                    G.add_node(process, node_type='process')
                G.add_edge(gene_symbol, process, weight=confidence)
        
    def create_function_network(self, df: pd.DataFrame):
        """Create a network visualization showing gene functions and their relationships"""
        G = nx.Graph()
        
        # Create nodes and edges
        for _, row in df.iterrows():
            gene_symbol = row['symbol']
            G.add_node(gene_symbol, node_type='gene')
            
            for function, confidence in row['molecular_functions']:
                if not G.has_node(function):
                    G.add_node(function, node_type='molecular')
                G.add_edge(gene_symbol, function, weight=confidence)
            
            for process, confidence in row['biological_processes']:
                if not G.has_node(process):
                    G.add_node(process, node_type='process')
                G.add_edge(gene_symbol, process, weight=confidence)
        
        # Create visualization with adjusted legend
        plt.figure(figsize=(15, 15))
        pos = nx.spring_layout(G, k=1, iterations=50)
        
        # Draw different types of nodes with different colors
        genes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'gene']
        functions = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'molecular']
        processes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'process']
        
        # Create separate legend entries with adjusted spacing
        plt.plot([], [], 'o', color='lightblue', markersize=15, label='Genes')
        plt.plot([], [], 'o', color='lightgreen', markersize=15, label='Molecular Functions')
        plt.plot([], [], 'o', color='lightpink', markersize=15, label='Biological Processes')
        
        # Draw nodes
        nx.draw_networkx_nodes(G, pos, nodelist=genes, node_color='lightblue', node_size=2000)
        nx.draw_networkx_nodes(G, pos, nodelist=functions, node_color='lightgreen', node_size=1500)
        nx.draw_networkx_nodes(G, pos, nodelist=processes, node_color='lightpink', node_size=1500)
        
        # Draw edges with varying thickness based on confidence
        edge_weights = [G[u][v]['weight'] * 2 for u, v in G.edges()]
        nx.draw_networkx_edges(G, pos, width=edge_weights, alpha=0.5)
        
        # Add labels
        nx.draw_networkx_labels(G, pos)
        
        plt.title("Gene Function Network\n(Edge thickness indicates confidence)", size=16, pad=20)
        
        # Adjust legend position and spacing
        plt.legend(bbox_to_anchor=(1.15, 1), 
                loc='upper right', 
                borderaxespad=0,
                frameon=True,
                labelspacing=1.5)
        
        plt.axis('off')
        plt.tight_layout()
        plt.savefig('visualizations/gene_function_network.png', 
                    bbox_inches='tight',
                    dpi=300)
        plt.close()

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