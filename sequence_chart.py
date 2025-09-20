"""
Protein sequence visualization utilities for the computational geometry dashboard.
Defines a function to generate a 2D scatter chart for protein sequences.
"""

import plotly.io as pio
import plotly.graph_objects as go

def custom_sequence_chart(
    seq: str,
    plt_title: str,
    highlight_indices: list = None
) -> go.Figure:
    """
    Generate a 2D scatter chart to display a protein sequence.

    Highlights specified residue indices in magenta and larger size.

    Args:
        seq (str): Protein sequence (one-letter codes).
        plt_title (str): Chart title.
        highlight_indices (list, optional): List of indices to highlight.

    Returns:
        plotly.graph_objects.Figure: The generated sequence chart.
    """
    if highlight_indices is None:
        highlight_indices = []

    spacing = 0.2
    x = [i * spacing for i in range(len(seq))]
    y = [0] * len(seq)
    labels = list(seq)
    seq_len = len(seq)

    # Marker settings
    colors = [
        'magenta' if i in highlight_indices else 'dodgerblue'
        for i in range(seq_len)
    ]
    sizes = [
        12 if i in highlight_indices else 10
        for i in range(seq_len)
    ]

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=x,
        y=y,
        mode='markers+text',
        text=labels,
        textposition="top center",
        hoverinfo='skip',
        marker=dict(
            size=sizes,
            color=colors,
            symbol='square',
            line=dict(width=1, color='white')
        )
    ))

    # X-axis tick labels every 5 characters
    tick_vals = [i * spacing for i in range(0, seq_len, 5)]
    tick_text = [str(i) for i in range(0, seq_len, 5)]

    fig.update_layout(
        plot_bgcolor='#242424',
        paper_bgcolor='#242424',
        title=plt_title,
        xaxis=dict(
            showgrid=False,
            showticklabels=True,
            tickmode='array',
            tickvals=tick_vals,
            ticktext=tick_text,
            zeroline=False,
            range=[-1, 25],
            fixedrange=False,
            tickfont=dict(color='white')
        ),
        yaxis=dict(
            visible=False,
            fixedrange=True
        ),
        font=dict(color='white'),
        height=100,  # Thin strip
        margin=dict(t=50, l=30, r=30, b=10),
        dragmode='pan' if seq_len > 100 else False
    )

    return fig