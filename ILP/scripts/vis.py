import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from enum import Enum

class LayoutType(Enum):
    """Enum for different graph layout types."""
    SPRING = "spring"
    CIRCULAR = "circular"
    SHELL = "shell"
    KAMADA_KAWAI = "kamada_kawai"

def generate_kmers(k):
    """Generate all possible binary strings of length k."""
    return [format(i, f'0{k}b') for i in range(2 ** k)]

def follows(kmer, other_kmer):
    """Check if the suffix of kmer (k-1) matches the prefix (k-1) of other_kmer."""
    return kmer[1:] == other_kmer[:-1]

def create_graph(k, order):
    """Create a directed graph with nodes representing k-mers and directed edges based on the follows() function."""
    G = nx.DiGraph()
    kmers = generate_kmers(k)

    for idx, kmer in enumerate(kmers):
        G.add_node(idx, label=kmer, value=order[idx])

    for i, kmer1 in enumerate(kmers):
        for j, kmer2 in enumerate(kmers):
            if follows(kmer1, kmer2):
                G.add_edge(i, j)

    return G

def get_text_color(color):
    """Calculate whether to use black or white text based on node color brightness."""
    r, g, b, _ = color
    brightness = 0.299 * r + 0.587 * g + 0.114 * b
    return 'black' if brightness > 0.5 else 'white'

def get_layout(G, layout_type):
    """Return the graph layout based on the selected layout type."""
    if layout_type == LayoutType.SPRING:
        return nx.spring_layout(G)
    elif layout_type == LayoutType.CIRCULAR:
        return nx.circular_layout(G)
    elif layout_type == LayoutType.SHELL:
        return nx.shell_layout(G)
    elif layout_type == LayoutType.KAMADA_KAWAI:
        return nx.kamada_kawai_layout(G)
    else:
        raise ValueError(f"Unsupported layout type: {layout_type}")

def draw_graph(G, layout_type=LayoutType.SPRING):
    """Draw the graph with labeled nodes and color based on their order values."""
    pos = get_layout(G, layout_type)
    values = np.array([G.nodes[node]['value'] for node in G.nodes])

    unique_values = np.unique(values)
    second_max_value = unique_values[-2] if len(unique_values) > 1 else unique_values[0]
    max_value = unique_values[-1]

    norm = plt.Normalize(vmin=values.min(), vmax=second_max_value)
    cmap = plt.get_cmap('viridis')

    node_colors = [
        cmap(norm(value)) if value < max_value else (1, 1, 1, 1)
        for value in values
    ]

    node_labels = {node: G.nodes[node]['label'] for node in G.nodes}
    font_colors = [get_text_color(color) for color in node_colors]

    # Create a larger figure
    plt.figure(figsize=(20, 15))  # Adjusted figure size

    # Adjust node and font sizes to ensure text fits inside nodes
    node_size = 1000  # Increased node size
    font_size = 8  # Increased font size for better readability

    # Draw nodes with adjusted size
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_size, cmap=cmap)

    # Draw labels with appropriate text color and font size
    for i, (node, label) in enumerate(node_labels.items()):
        nx.draw_networkx_labels(
            G, pos, labels={node: label}, font_color=font_colors[i], font_size=font_size,
            verticalalignment='center', horizontalalignment='center'
        )

    # Draw directed edges with larger arrows
    nx.draw_networkx_edges(G, pos, edge_color='black', arrows=True, arrowstyle='->', arrowsize=30)

    # Add a colorbar with larger label
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(values)
    plt.colorbar(sm, label="Order Value (Max Excluded)", fraction=0.03, pad=0.04)

    # Show the plot
    plt.show()



import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from enum import Enum
from scripts.scripts import prob_gc  # Assuming this is part of your local setup

class LayoutType(Enum):
    """Enum for different graph layout types."""
    SPRING = "spring"
    CIRCULAR = "circular"
    SHELL = "shell"
    KAMADA_KAWAI = "kamada_kawai"

def generate_kmers(k):
    """Generate all possible binary strings of length k."""
    return [format(i, f'0{k}b') for i in range(2 ** k)]

def follows(kmer, other_kmer):
    """Check if the suffix of kmer (k-1) matches the prefix (k-1) of other_kmer."""
    return kmer[1:] == other_kmer[:-1]

def create_graph(k, order):
    """Create a directed graph with nodes representing k-mers and directed edges based on the follows() function."""
    G = nx.DiGraph()
    kmers = generate_kmers(k)

    for idx, kmer in enumerate(kmers):
        G.add_node(idx, label=kmer, value=order[idx])

    for i, kmer1 in enumerate(kmers):
        for j, kmer2 in enumerate(kmers):
            if follows(kmer1, kmer2):
                G.add_edge(i, j)

    return G



def draw_filtered_graph(G, layout_type=LayoutType.SPRING):
    """Draw the graph, excluding nodes with the maximum unique value."""
    pos = get_layout(G, layout_type)
    values = np.array([G.nodes[node]['value'] for node in G.nodes])

    # Identify the maximum unique value
    unique_values = np.unique(values)
    max_value = unique_values[-1]

    # Filter nodes to exclude those with the max value
    filtered_nodes = [node for node in G.nodes if G.nodes[node]['value'] < max_value]

    # Filter edges to only include edges between the remaining nodes
    filtered_edges = [
        (u, v) for u, v in G.edges if u in filtered_nodes and v in filtered_nodes
    ]

    # Normalize the values for color mapping (without max value)
    norm = plt.Normalize(vmin=values.min(), vmax=unique_values[-2])
    cmap = plt.get_cmap('viridis')

    # Assign colors to the remaining nodes
    node_colors = [
        cmap(norm(G.nodes[node]['value'])) for node in filtered_nodes
    ]
    node_labels = {node: G.nodes[node]['label'] for node in filtered_nodes}
    font_colors = [get_text_color(color) for color in node_colors]

    # Create a larger figure
    plt.figure(figsize=(20, 15))

    # Draw only the filtered nodes and edges
    nx.draw_networkx_nodes(G, pos, nodelist=filtered_nodes, node_color=node_colors, 
                           node_size=1000, cmap=cmap)
    nx.draw_networkx_edges(G, pos, edgelist=filtered_edges, edge_color='black', 
                           arrows=True, arrowstyle='->', arrowsize=25)

    # Draw labels with appropriate font size
    for i, (node, label) in enumerate(node_labels.items()):
        nx.draw_networkx_labels(G, pos, labels={node: label}, 
                                font_color=font_colors[i], font_size=10,
                                verticalalignment='center', horizontalalignment='center')

    # Add a colorbar to the plot
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(values)
    plt.colorbar(sm, label="Order Value (Max Excluded)", fraction=0.03, pad=0.04)

    # Show the plot
    plt.show()

def get_text_color(color):
    """Calculate whether to use black or white text based on node color brightness."""
    r, g, b, _ = color
    brightness = 0.299 * r + 0.587 * g + 0.114 * b
    return 'black' if brightness > 0.5 else 'white'

def get_layout(G, layout_type):
    """Return the graph layout based on the selected layout type."""
    if layout_type == LayoutType.SPRING:
        return nx.spring_layout(G)
    elif layout_type == LayoutType.CIRCULAR:
        return nx.circular_layout(G)
    elif layout_type == LayoutType.SHELL:
        return nx.shell_layout(G)
    elif layout_type == LayoutType.KAMADA_KAWAI:
        return nx.kamada_kawai_layout(G)
    else:
        raise ValueError(f"Unsupported layout type: {layout_type}")
    

from scripts.scripts import prob_gc
def get_uhs_order(order,w ,k):
    # first validate the order

    # Step 1: Pair each element with its index
    value_index_list = list(enumerate(order))

    # Step 2: Sort the list based on values (ties broken by index)
    sorted_value_index_list = sorted(value_index_list, key=lambda x: x[1])

    # Step 3: Assign ranks and place them back into a list matching the original order
    ranks = [0] * len(order)
    for rank, (index, value) in enumerate(sorted_value_index_list, start=0):
        ranks[index] = rank

    #print(ranks)

    gc_count_initial = prob_gc(w,k ,order)

    for i in range(len(ranks)):
        # make copy of ranks
        ranks_copy = ranks.copy()
        # all cells containing more than i are set to 2^k
        for j in range(len(ranks_copy)):
            if ranks_copy[j] > i:
                ranks_copy[j] = 2 ** k
        gc_count_now = prob_gc(w,k ,ranks_copy)
        if gc_count_now == gc_count_initial:
            #print(i)
            return ranks_copy
