import plotly.graph_objects as go

# Flow data
flow_through = [
    ("CILL", "ET", "VAP", 1), ("CILL", "ND", "VAP", 3),
    ("Trauma", "ET", "VAP", 1), ("TraumaNN", "ND", "VAP", 1),
    ("SCI", "ET", "VAP", 1), ("TBI", "ET", "VAP", 1),
    ("Stroke", "ET", "VAP", 1), ("C19", "ET", "VAP", 1), ("C19", "ND", "VAP", 1),

    ("CILL", "ET", "MVD", 3), ("CILL", "ND", "MVD", 1),
    ("Trauma", "ET", "MVD", 1), ("TraumaNN", "ND", "MVD", 1),
    ("SCI", "ET", "MVD", 1), ("TBI", "ET", "MVD", 1),
    ("Stroke", "ET", "MVD", 1), ("C19", "ND", "MVD", 3),

    ("CILL", "ET", "STM", 2), ("CILL", "ND", "STM", 3),
    ("SCI", "ND", "STM", 1), ("C19", "ND", "STM", 2),

    ("Trauma", "ET", "ICULOS", 1), ("TraumaNN", "ND", "ICULOS", 1),
    ("SCI", "ET", "ICULOS", 1), ("TBI", "ET", "ICULOS", 1),
    ("Stroke", "ET", "ICULOS", 1), ("C19", "ND", "ICULOS", 2),
    ("CILL", "ET", "ICULOS", 4),

    ("CILL", "ND", "HLOS", 1), ("TraumaNN", "ND", "HLOS", 1),
    ("SCI", "ET", "HLOS", 1), ("TBI", "ET", "HLOS", 1),
    ("Stroke", "ET", "HLOS", 1), ("C19", "ND", "HLOS", 1),
]

# Source-based color mapping
source_colors = {
    "CILL": "rgba(27, 158, 119, 1)",
    "Trauma": "rgba(217, 95, 2, 1)",
    "TraumaNN": "rgba(117, 112, 179, 1)",
    "Stroke": "rgba(230, 171, 2, 1)",
    "TBI": "rgba(102, 166, 30, 1)",
    "C19": "rgba(166, 118, 29, 1)",
    "SCI": "rgba(231, 41, 138, 1)"
}

# Expand flow data
expanded_rows = []
for src, mid, dst, val in flow_through:
    for _ in range(val):
        expanded_rows.append((src, mid, dst))

sources = [r[0] for r in expanded_rows]
middles = [r[1] for r in expanded_rows]
targets = [r[2] for r in expanded_rows]
color_vals = [source_colors[src] for src in sources]

# Build the figure
fig = go.Figure(go.Parcats(
    dimensions=[
        dict(values=sources, label="Source"),
        dict(values=middles, label="Mid (ET/ND)"),
        dict(values=targets, label="Destination")
    ],
    line=dict(color=color_vals),
    hoveron='color',
    hoverinfo='count+probability',
    arrangement='freeform',
    bundlecolors=True,           # Smooth same-color lines
    sortpaths='forward',         # Order for smooth transitions
    labelfont=dict(size=18),
    tickfont=dict(size=16)
))

fig.update_layout(
    title="Smoothed Flow from Source → Mid (ET/ND) → Destination",
    font=dict(size=16)
)

fig.show()
