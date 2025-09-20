import math
import pandas as pd
import plotly.express as px
from Bio.PDB import PDBParser, PPBuilder

def custom_ramachandran_chart(pdb_id, pdb_path, highlight_indices=None):
    """
    Generate a Ramachandran Chart for Φ (Phi) and Ψ (Psi) angles in a protein chain.

    Only residues with both angles are displayed. Highlighted residues are colored and sized differently.
    If highlight_indices contains residues not present in the plot, a warning annotation is added.

    Args:
        pdb_id (str): PDB identifier.
        pdb_path (str): Path to PDB file.
        highlight_indices (list, optional): List of residue numbers to highlight.

    Returns:
        plotly.graph_objects.Figure: The generated Ramachandran chart.
    """
    if highlight_indices is None:
        highlight_indices = []

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_path)
    phi_psi_data = []

    for model in structure:
        for chain in model:
            ppb = PPBuilder()
            for pp in ppb.build_peptides(chain):
                phi_psi = pp.get_phi_psi_list()
                for res, (phi, psi) in zip(pp, phi_psi):
                    if phi and psi:
                        res_num = res.get_id()[1]
                        phi_deg = math.degrees(phi)
                        psi_deg = math.degrees(psi)
                        res_name = res.get_resname()
                        phi_psi_data.append({
                            "Residue": res_name,
                            "Residue Number": res_num,
                            "Phi": phi_deg,
                            "Psi": psi_deg,
                            "Chain": chain.id
                        })

    df_phi_psi = pd.DataFrame(phi_psi_data)

    # Set highlight category and sizes
    df_phi_psi["Selection"] = df_phi_psi["Residue Number"].apply(
        lambda x: "Selection" if x in highlight_indices else "Residues"
    )
    df_phi_psi["Size"] = df_phi_psi["Residue Number"].apply(
        lambda x: 12 if x in highlight_indices else 10
    )

    # Track missing residue numbers
    actual_residues = set(df_phi_psi["Residue Number"])
    missing = sorted(set(highlight_indices) - actual_residues)
    print(missing)

    # Build the chart
    fig = px.scatter(
        df_phi_psi,
        x='Phi',
        y='Psi',
        color="Selection",
        size="Size",
        hover_data=['Residue', "Residue Number", "Chain"],
        title="Ramachandran Chart",
        labels={'Phi': 'Φ (Phi)', 'Psi': 'Ψ (Psi)'},
        color_discrete_map={
            "Residues": "dodgerblue",
            "Selection": "magenta"
        }
    )

    fig.update_layout(
        margin=dict(l=50, r=50, t=50, b=50),
        font=dict(color='gray'),
        legend=dict(title=None),
        plot_bgcolor='#242424',
        paper_bgcolor="#242424",
        xaxis_range=[-180, 180],
        yaxis_range=[-180, 180],
        xaxis=dict(tickmode='array', tickvals=list(range(-180, 181, 60))),
        yaxis=dict(tickmode='array', tickvals=list(range(-180, 181, 60))),
        height=540,
    )

    fig.update_xaxes(gridcolor='grey') 
    fig.update_yaxes(gridcolor='grey') 

    # Add annotation if residues are missing
    if missing:
        fig.add_annotation(
            text="⚠️ First & Last residues not displayed.",
            xref="paper", yref="paper",
            x=0, y=-0.08,
            showarrow=False,
            font=dict(color="dodgerblue", size=12),
            align="left",
            bgcolor="rgba(255,255,255,0.1)",
            bordercolor="dodgerblue",
            borderwidth=1,
        )

    return fig