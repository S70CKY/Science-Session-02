import plotly.graph_objects as go

# Nodes
nodes = [
    "Exercise", "Diet", "Medication",
    "Genomic", "Epigenomic", "Transcriptomic", "Metabolomic", "Proteomic", "Metagenomic",
    "Muscle Function", "Neurological Recovery", "Metabolic Improvement", 
    "Immune Function", "Autonomic Function", "Altering Microbiome"
]

# create Links (1 Node, 2 Node, Value (weight))
links = [
    # Exercise -> Omics
    ("Exercise", "Epigenomic", 3),
    ("Exercise", "Transcriptomic", 11),
    ("Exercise", "Metagenomic", 1),

    # Diet -> Omics
    ("Diet", "Metagenomic", 2),
    ("Diet", "Metabolomic", 1),

    # Medication -> Omics
    ("Medication", "Proteomic", 2),
    ("Medication", "Metagenomic", 1),
    ("Medication", "Transcriptomic", 1),
    ("Medication", "Metabolomic", 1),
    
    # Omics -> Outcomes
    ("Metagenomic", "Immune Function", 1),
    ("Metagenomic", "Altering Microbiome", 2),
    ("Metagenomic", "Metabolic Improvement", 0.5),
    ("Metagenomic", "Neurological Recovery", 0.5),
    ("Epigenomic", "Muscle Function", 1.5),
    ("Epigenomic", "Neurological Recovery", 0.5),
    ("Epigenomic", "Metabolic Improvement", 1),
    ("Transcriptomic", "Muscle Function", 5),
    ("Transcriptomic", "Neurological Recovery", 1.37),
    ("Transcriptomic", "Metabolic Improvement", 3),
    ("Transcriptomic", "Immune Function", 1.30),
    ("Metabolomic", "Neurological Recovery", 1),
    ("Metabolomic", "Metabolic Improvement", 1),
    ("Proteomic", "Neurological Recovery", 1),
    ("Proteomic", "Metabolic Improvement", 0.5),
    ("Proteomic", "Immune Function", 0.5),
    ("Transcriptomic", "Autonomic Function", 0.68), 
]

# Nodes zu indices zuordnen
node_dict = {node: idx for idx, node in enumerate(nodes)}

# define colors and opacity -> 40%
node_colors = {
    "Exercise": "rgba(0, 194, 249, 0.4)",  
    "Diet": "rgba(0, 149, 3, 0.4)",      
    "Medication": "rgba(0, 145, 117, 0.4)",
    "Metagenomic": "rgba(120, 0, 75, 0.4)", 
    "Genomic": "rgba(90, 0, 15, 0.4)",   
    "Epigenomic": "rgba(205, 2, 45, 0.4)",  
    "Transcriptomic": "rgba(255, 172, 59, 0.4)",  
    "Metabolomic": "rgba(255, 99, 132, 0.4)",  
    "Proteomic": "rgba(132, 0, 205, 0.4)",  
    "Muscle Function": "rgba(75, 192, 192, 0.4)",  
    "Neurological Recovery": "rgba(75, 192, 192, 0.4)",  
    "Metabolic Improvement": "rgba(75, 192, 192, 0.4)",  
    "Immune Function": "rgba(75, 192, 192, 0.4)",  
    "Autonomic Function": "rgba(75, 192, 192, 0.4)", 
    "Altering Microbiome": "rgba(75, 192, 192, 0.4)"  
}

# Data arrays für Sankey
source = []
target = []
value = []
link_colors = []

# create Arrays for the links
for link in links:
    source.append(node_dict[link[0]])  # index of source node
    target.append(node_dict[link[1]])  # index of target node
    value.append(link[2])  # value 
    link_colors.append(node_colors[link[0]])  # link color to match node color

# Assign colors to nodes with opacity
node_borders = [node_colors[node] for node in nodes]

# Create the Sankey diagram
fig = go.Figure(go.Sankey(
    arrangement="snap", 
    node=dict(
        pad=50,  # spacing 
        thickness=200,  # thickness of nodes 
        line=dict(color=node_borders, width=7),  # border colors
        color=[node_colors[node] for node in nodes],  # Node fill colors with opacity
        label=["" for _ in nodes],  # Remove labels
    ),
    link=dict(
        source=source,
        target=target,
        value=value,
        color=link_colors,
        hovertemplate=" ",  
    )
))

# Adjust the layout
fig.update_layout(
    # title="",
    autosize=False,
    # DinA4 size 
    width=2970,  
    height=2100,  
    margin=dict(l=100, r=100, t=100, b=100),
    font=dict(size=12)
)

fig.show()