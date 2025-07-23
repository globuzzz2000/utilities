import math
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display, clear_output

def prob_only_i(lambdas_per_droplet, i):
    """Calculate probability of droplets only positive for target i"""
    lambda_i = lambdas_per_droplet[i]
    others = lambdas_per_droplet[:i] + lambdas_per_droplet[i+1:]
    product_exp = math.prod([math.exp(-l) for l in others])
    return (1 - math.exp(-lambda_i)) * product_exp

def median(lst):
    sorted_lst = sorted(lst)
    n = len(lst)
    mid = n // 2
    return (sorted_lst[mid] if n % 2 else (sorted_lst[mid - 1] + sorted_lst[mid]) / 2)

def calculate_results(total_droplets, copies_per_sample):
    """Calculate all results for given parameters"""
    # Convert to per-droplet λᵢ values
    lambdas = [copies / total_droplets for copies in copies_per_sample]
    
    # Calculate droplet counts
    only_counts = []
    for i in range(5):
        p_only = prob_only_i(lambdas, i)
        count = total_droplets * p_only
        only_counts.append(count)
    
    # Normalization relative to average within 15% of the median
    median_count = median(only_counts)
    allowed_range = (0.85 * median_count, 1.15 * median_count)
    in_range = [c for c in only_counts if allowed_range[0] <= c <= allowed_range[1]]
    reference_mean = sum(in_range) / len(in_range) if in_range else float('nan')
    
    # Calculate relative abundances
    relative_abundances = [count / reference_mean if reference_mean else float('nan') 
                          for count in only_counts]
    
    return {
        'lambdas': lambdas,
        'only_counts': only_counts,
        'median_count': median_count,
        'allowed_range': allowed_range,
        'reference_mean': reference_mean,
        'relative_abundances': relative_abundances
    }

def print_results(results, copies_per_sample):
    """Print formatted results"""
    print("--- Droplet Counts ---")
    for i in range(5):
        print(f"Chromosome {i+1}: {results['only_counts'][i]:.1f} droplets")
    
    print(f"\nMedian count: {results['median_count']:.1f}")
    print(f"Allowed range: {results['allowed_range'][0]:.1f} - {results['allowed_range'][1]:.1f}")
    print(f"Reference mean: {results['reference_mean']:.1f}")
    
    print("\n--- Relative Abundance ---")
    for i, rel in enumerate(results['relative_abundances']):
        in_range = (results['allowed_range'][0] <= results['only_counts'][i] <= results['allowed_range'][1])
        status = "✓" if in_range else "✗"
        print(f"Chromosome {i+1}: {rel:.3f} {status}")

def create_chr4_abundance_plot():
    """Create the specific plot showing chromosome 4 abundance vs average copies"""
    
    # Parameters
    total_droplets = 19000
    avg_copies_range = np.linspace(1000, 10000, 100)
    
    # Two scenarios
    ratio_125 = [1, 1, 1, 1.25, 1]  # Chr4 = 1.25x
    ratio_075 = [1, 1, 1, 0.75, 1]  # Chr4 = 0.75x
    
    chr4_abundance_125 = []
    chr4_abundance_075 = []
    
    for avg_copies in avg_copies_range:
        # Scenario 1: Chr4 = 1.25x
        copies_125 = [avg_copies * r for r in ratio_125]
        results_125 = calculate_results(total_droplets, copies_125)
        chr4_abundance_125.append(results_125['relative_abundances'][3])  # Chr4 is index 3
        
        # Scenario 2: Chr4 = 0.75x
        copies_075 = [avg_copies * r for r in ratio_075]
        results_075 = calculate_results(total_droplets, copies_075)
        chr4_abundance_075.append(results_075['relative_abundances'][3])  # Chr4 is index 3
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 1, 1)
    plt.plot(avg_copies_range, chr4_abundance_125, 'b-', linewidth=2, label='Chr4 = 1.25x (excess)')
    plt.axhline(y=1.25, color='b', linestyle='--', alpha=0.7, label='Expected (1.25)')
    plt.xlabel('Average Copies per Sample')
    plt.ylabel('Chr4 Relative Abundance')
    plt.title('Chromosome 4 Abundance vs Average Copies (Chr4 = 1.25x scenario)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.subplot(2, 1, 2)
    plt.plot(avg_copies_range, chr4_abundance_075, 'r-', linewidth=2, label='Chr4 = 0.75x (deficit)')
    plt.axhline(y=0.75, color='r', linestyle='--', alpha=0.7, label='Expected (0.75)')
    plt.xlabel('Average Copies per Sample')
    plt.ylabel('Chr4 Relative Abundance')
    plt.title('Chromosome 4 Abundance vs Average Copies (Chr4 = 0.75x scenario)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plt.show()
    
    return avg_copies_range, chr4_abundance_125, chr4_abundance_075

def interactive_calculator():
    """Create interactive widgets for the calculator"""
    
    # Create widgets
    total_droplets_widget = widgets.IntSlider(
        value=19000, min=5000, max=50000, step=1000,
        description='Total Droplets:'
    )
    
    chr1_widget = widgets.IntSlider(value=5000, min=1000, max=10000, step=250, description='Chr1 Copies:')
    chr2_widget = widgets.IntSlider(value=5000, min=1000, max=10000, step=250, description='Chr2 Copies:')
    chr3_widget = widgets.IntSlider(value=5000, min=1000, max=10000, step=250, description='Chr3 Copies:')
    chr4_widget = widgets.IntSlider(value=6250, min=1000, max=10000, step=250, description='Chr4 Copies:')
    chr5_widget = widgets.IntSlider(value=5000, min=1000, max=10000, step=250, description='Chr5 Copies:')
    
    def update_results(total_droplets, chr1, chr2, chr3, chr4, chr5):
        clear_output(wait=True)
        
        copies_per_sample = [chr1, chr2, chr3, chr4, chr5]
        results = calculate_results(total_droplets, copies_per_sample)
        
        print_results(results, copies_per_sample)
        
        # Create a simple bar plot
        plt.figure(figsize=(12, 4))
        
        plt.subplot(1, 2, 1)
        chromosomes = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']
        colors = ['green' if results['allowed_range'][0] <= count <= results['allowed_range'][1] 
                 else 'red' for count in results['only_counts']]
        plt.bar(chromosomes, results['only_counts'], color=colors, alpha=0.7)
        plt.ylabel('Droplet Count')
        plt.title('Corrected Droplet Counts')
        plt.xticks(rotation=45)
        
        plt.subplot(1, 2, 2)
        plt.bar(chromosomes, results['relative_abundances'], color='steelblue', alpha=0.7)
        plt.axhline(y=1, color='black', linestyle='--', alpha=0.5)
        plt.ylabel('Relative Abundance')
        plt.title('Relative Abundance (normalized)')
        plt.xticks(rotation=45)
        
        plt.tight_layout()
        plt.show()
    
    # Create interactive output
    interactive_output = widgets.interactive_output(
        update_results,
        {
            'total_droplets': total_droplets_widget,
            'chr1': chr1_widget,
            'chr2': chr2_widget,
            'chr3': chr3_widget,
            'chr4': chr4_widget,
            'chr5': chr5_widget
        }
    )
    
    # Display widgets
    display(widgets.VBox([
        total_droplets_widget,
        widgets.HBox([chr1_widget, chr2_widget, chr3_widget]),
        widgets.HBox([chr4_widget, chr5_widget]),
        interactive_output
    ]))

# Main execution functions
def main():
    """Run the original calculator"""
    total_droplets = 19000
    copies_per_sample = [5000, 5000, 5000, 6250, 5000]
    
    results = calculate_results(total_droplets, copies_per_sample)
    print_results(results, copies_per_sample)

if __name__ == "__main__":
    print("=== Poisson Droplet Calculator ===\n")
    
    # Run original calculation
    print("1. Original calculation:")
    main()
    
    print("\n" + "="*50 + "\n")
    
    # Create the chromosome 4 abundance plot
    print("2. Creating Chromosome 4 abundance plot...")
    create_chr4_abundance_plot()
    
    print("\n" + "="*50 + "\n")
    
    # Instructions for interactive mode
    print("3. For interactive mode, run:")
    print("   interactive_calculator()")
    print("\nThis requires Jupyter notebook/lab with ipywidgets installed:")
    print("   pip install ipywidgets matplotlib numpy")