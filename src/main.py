import requests
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import time
from typing import Dict, List, Set
from pathlib import Path

class GeneAnalyzer:
    def __init__(self):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.delay = 0.34
        Path("data").mkdir(exist_ok=True)
        Path("visualizations").mkdir(exist_ok=True)
        sns.set_theme(style="whitegrid")
        
        # Define molecular functions and biological processes
        self.gene_functions = {
            'BDNF': {
                'molecular_functions': [
                    ('growth_factor_activity', 0.9),
                    ('protein_binding', 0.8),
                    ('receptor_binding', 0.7)
                ],
                'biological_processes': [
                    ('neuron_development', 0.9),
                    ('synaptic_plasticity', 0.8),
                    ('signal_transduction', 0.7)
                ]
            },
            'DRD2': {
                'molecular_functions': [
                    ('g_protein_coupled_receptor', 0.9),
                    ('dopamine_binding', 0.9),
                    ('signal_transducer', 0.8)
                ],
                'biological_processes': [
                    ('dopamine_signaling', 0.9),
                    ('synaptic_transmission', 0.8),
                    ('behavior_regulation', 0.7)
                ]
            }
        }

    def fetch_gene_data(self, gene_symbol: str) -> Dict:
        """Fetch gene information from NCBI"""
        try:
            # First, search for the gene
            search_url = f"{self.base_url}esearch.fcgi?db=gene&term={gene_symbol}[Gene Name]+AND+homo+sapiens[Organism]&format=json"
            response = requests.get(search_url)
            response.raise_for_status()
            search_data = response.json()
            
            if (search_data.get('esearchresult', {}).get('count', '0') != '0' and 
                search_data['esearchresult']['idlist']):
                
                gene_id = search_data['esearchresult']['idlist'][0]
                print(f"Found gene ID for {gene_symbol}: {gene_id}")
                
                # Get detailed information
                time.sleep(self.delay)
                summary_url = f"{self.base_url}esummary.fcgi?db=gene&id={gene_id}&format=json"
                summary_response = requests.get(summary_url)
                summary_response.raise_for_status()
                summary_data = summary_response.json()
                
                gene_info = summary_data['result'][gene_id]
                summary = gene_info.get('summary', '')
                
                # Get functions with confidence scores
                functions = self.get_gene_functions(gene_symbol, summary)
                
                return {
                    'symbol': gene_symbol,
                    'gene_id': gene_id,
                    'name': gene_info.get('description', 'N/A'),
                    'chromosome': gene_info.get('chromosome', 'N/A'),
                    'summary': summary,
                    'molecular_functions': functions.get('molecular_functions', []),
                    'biological_processes': functions.get('biological_processes', [])
                }
            
            print(f"No results found for {gene_symbol}")
            return None
                
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data for {gene_symbol}: {e}")
            return None
    
    def get_gene_functions(self, gene_symbol: str, summary: str) -> Dict:
        """Get detailed gene functions with confidence scores"""
        base_functions = self.gene_functions.get(gene_symbol, {})
        if not base_functions:
            return self.categorize_from_summary(summary)
        return base_functions

    def categorize_from_summary(self, summary: str) -> Dict:
        """Categorize gene functions from summary text with confidence scores"""
        summary = summary.lower()
        
        molecular_keywords = {
            'receptor': ('receptor_activity', 0.8),
            'enzyme': ('enzyme_activity', 0.8),
            'channel': ('ion_channel_activity', 0.8),
            'transcription': ('transcription_regulation', 0.8),
            'binding': ('protein_binding', 0.7)
        }
        
        process_keywords = {
            'development': ('neural_development', 0.8),
            'signaling': ('signal_transduction', 0.8),
            'synaptic': ('synaptic_function', 0.8),
            'plasticity': ('neural_plasticity', 0.8),
            'transmission': ('neurotransmission', 0.8)
        }
        
        molecular_functions = []
        biological_processes = []
        
        for keyword, (function, base_score) in molecular_keywords.items():
            if keyword in summary:
                confidence = base_score * (1 + summary.count(keyword) * 0.1)
                molecular_functions.append((function, min(confidence, 0.9)))
        
        for keyword, (process, base_score) in process_keywords.items():
            if keyword in summary:
                confidence = base_score * (1 + summary.count(keyword) * 0.1)
                biological_processes.append((process, min(confidence, 0.9)))
        
        return {
            'molecular_functions': molecular_functions,
            'biological_processes': biological_processes
        }

            # Create visualization with adjusted legend
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
    try:
        print("Starting gene analysis...")
        analyzer = GeneAnalyzer()
        
        genes = [
            'BDNF', 'DRD2', 'COMT', 'HTR2A', 'FOXP2',
            'SLC6A4', 'APOE', 'SNAP25', 'CACNA1C', 'NRG1'
        ]
        
        results = []
        for gene in genes:
            print(f"\nProcessing {gene}...")
            data = analyzer.fetch_gene_data(gene)
            if data:
                results.append(data)
                print(f"\nFunctions for {gene}:")
                print("Molecular functions:")
                for func, conf in data['molecular_functions']:
                    print(f"  - {func} (confidence: {conf:.2f})")
                print("Biological processes:")
                for proc, conf in data['biological_processes']:
                    print(f"  - {proc} (confidence: {conf:.2f})")
        
        if results:
            df = pd.DataFrame(results)
            analyzer.create_function_network(df)
            df.to_csv('data/gene_analysis_with_functions.csv', index=False)
            print("\nResults and visualizations saved.")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        import traceback
        print(traceback.format_exc())

if __name__ == "__main__":
    main()