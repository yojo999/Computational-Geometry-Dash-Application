import numpy as np
import plotly.graph_objects as go
from pdb_parser import biopandas_pdb_parser
from scipy.spatial import ConvexHull, Delaunay
from visualization_helper import get_edges_tri, get_edges_tetra, arrays_of_selected

def della(coords_df, aa_seq_lst):
    """
    Perform Delaunay tessellation on coordinates and map results to amino acid sequence.

    Args:
        coords_df (np.ndarray): Array of coordinates.
        aa_seq_lst (list): List of amino acid sequence labels.

    Returns:
        tuple: (simplices, coplanar, aa labels, index labels, simplex coordinates)
    """
    della = Delaunay(coords_df, qhull_options='Qc QJ')
    sx_num_list = [list(i) for i in della.simplices]
    sx_aa_lst = ["-".join([aa_seq_lst[i] for i in sublist]) for sublist in sx_num_list]
    sx_num_lst_della = ["-".join([str(i) for i in sublist]) for sublist in sx_num_list]
    sx_coords_della = coords_df[della.simplices]
    return della.simplices, della.coplanar, sx_aa_lst, sx_num_lst_della, sx_coords_della

def cxhull(coords_df, aa_seq_lst):
    """
    Compute the convex hull of coordinates and map results to amino acid sequence.

    Args:
        coords_df (np.ndarray): Array of coordinates.
        aa_seq_lst (list): List of amino acid sequence labels.

    Returns:
        tuple: (area, volume, hull simplices, aa labels, index labels, simplex coordinates)
    """
    cx_hull = ConvexHull(coords_df, qhull_options='Qc QJ')
    tri_num_lst_cx = [list(i) for i in cx_hull.simplices]
    tri_aa_lst_cx = ["-".join([aa_seq_lst[i] for i in sublist]) for sublist in tri_num_lst_cx]
    tri_num_lst_cx_str = ["-".join([str(i) for i in sublist]) for sublist in tri_num_lst_cx]
    tri_coords_cx = coords_df[cx_hull.simplices]
    return cx_hull.area, cx_hull.volume, cx_hull.simplices, tri_aa_lst_cx, tri_num_lst_cx_str, tri_coords_cx

def cg_chart(pdb_id, pdb_ext, chain_id, pdb_loc, highlight_indices=None):
    """
    Generate a computational geometry chart for a 3D protein chain.

    Plots C-alpha trace, all atom residues, Delaunay tessellation, and convex hull.

    Args:
        pdb_id (str): PDB identifier.
        pdb_ext (str): PDB file extension.
        chain_id (str): Protein chain identifier.
        pdb_loc (str): Path to PDB file.
        highlight_indices (list, optional): Residue indices to highlight.

    Returns:
        plotly.graph_objects.Figure: The generated 3D chart.
    """
    hed, df_ca, df_aa, bf_ca, bf_aa, xyz_ca, xyz_aa, aa_lst_ca, aa_lst_aa, aa1_ca, aa1_aa = biopandas_pdb_parser(pdb_id, pdb_loc)
    sx_dca, cplnr_dca, sx_aa_dca, sx_num_dca, sx_xyz_dca = della(xyz_ca, aa_lst_ca)
    sx_daa, cplnr_daa, sx_aa_daa, sx_num_daa, sx_xyz_daa = della(xyz_aa, aa_lst_aa)

    x_min, x_max = xyz_aa[:, 0].min(), xyz_aa[:, 0].max()
    y_min, y_max = xyz_aa[:, 1].min(), xyz_aa[:, 1].max()
    z_min, z_max = xyz_aa[:, 2].min(), xyz_aa[:, 2].max()

    cx_area_ca, cx_vol_ca, cx_tri_ca, cx_tri_amino_ca, cx_tri_num_ca, cx_xyz_ca = cxhull(xyz_ca, aa_lst_ca)
    cx_area_aa, cx_vol_aa, cx_tri_aa, cx_tri_amino_aa, cx_tri_num_aa, cx_xyz_aa = cxhull(xyz_aa, aa_lst_aa)

    sx_xedge_ca, sx_yedge_ca, sx_zedge_ca = get_edges_tetra(sx_xyz_dca)
    sx_xedge_aa, sx_yedge_aa, sx_zedge_aa = get_edges_tetra(sx_xyz_daa)

    cx_xedge_ca, cx_yedge_ca, cx_zedge_ca = get_edges_tri(cx_xyz_ca)
    cx_xedge_aa, cx_yedge_aa, cx_zedge_aa = get_edges_tri(cx_xyz_aa)

    labels_ca = list(df_ca['label'])
    labels_aa = list(df_aa['label'])

    marker_sizes = [5] * len(df_ca)
    marker_colors = ['dodgerblue'] * len(df_ca)

    if highlight_indices:
        highlight_indices = [(i-1) for i in highlight_indices]
        sel_simplices = arrays_of_selected(sx_daa, highlight_indices)
        unique_select = np.unique(sel_simplices)
        sel_coords = xyz_aa[sel_simplices]
        sel_xedge_aa, sel_yedge_aa, sel_zedge_aa = get_edges_tetra(sel_coords)
        for res in unique_select:
            if 0 <= res < len(df_ca):
                marker_sizes[res] = 12
                marker_colors[res] = 'magenta'

    plt_title = f"{pdb_id[0:4]}  chain {chain_id}"
    plt_sub_title = (
        "▪ C-alpha residues: {}.<br>"
        "▪ Total residues: {}.<br>"
        "▪ {} c-alpha and {} all atom Convex Hull simplices.<br>"
        "▪ {} c-alpha and {} all atom Delaunay tetrahedra.<br>"
        "▪ all-atom convex hull surface area {} sq. A.<br>"
        "▪ all-atom convex hull volume {} cu. A."
        ).format(len(xyz_ca), len(xyz_aa), len(cx_tri_ca), len(cx_tri_aa),
                 len(sx_dca), len(sx_daa), round(cx_area_aa, 2), round(cx_vol_aa, 2))

    ca_trace = go.Scatter3d(
        x=df_ca['x_coord'], 
        y=df_ca['y_coord'], 
        z=df_ca['z_coord'],
        text=labels_ca,
        mode='markers+lines',
        name='C-alpha trace',
        line=dict(color=marker_colors, width=6),
        showlegend=True
    )

    all_atm_trace = go.Scatter3d(
        x=df_aa['x_coord'], 
        y=df_aa['y_coord'], 
        z=df_aa['z_coord'],
        text=labels_aa,
        mode='markers',
        name='all atom vertices',
        marker=dict(color='papayawhip', size=4, symbol='circle'),
        visible='legendonly',
        showlegend=True
    )

    fig = go.Figure(data=[ca_trace, all_atm_trace])

    fig.add_trace(go.Scatter3d(
        x=sx_xedge_ca, y=sx_yedge_ca, z=sx_zedge_ca,
        mode='lines', name='C-alpha edges',
        line=dict(color='whitesmoke', width=2),
    ))

    fig.add_trace(go.Scatter3d(
        x=sx_xedge_aa, y=sx_yedge_aa, z=sx_zedge_aa,
        mode='lines', name='all atom edges',
        line=dict(color='lightblue', width=2),
        visible='legendonly',
        showlegend=True
    ))

    for i, pt in enumerate(cx_xyz_aa):
        fig.add_trace(go.Scatter3d(
            x=pt[:, 0], y=pt[:, 1], z=pt[:, 2],
            mode='markers',
            name='points' if i == 0 else '',
            opacity=0.1,
            legendgroup='looped',
            showlegend=(i == 0),
            hoverinfo='skip',
            marker=dict(size=3, color='lightblue', opacity=1)
        ))

    for i, pt in enumerate(cx_xyz_aa):
        fig.add_trace(go.Mesh3d(
            x=pt[:, 0], y=pt[:, 1], z=pt[:, 2],
            name='surface' if i == 0 else '',
            alphahull=-1,
            color='lightyellow',
            flatshading=False,
            opacity=0.1,
            legendgroup='looped',
            legendgrouptitle=dict(text='Convex Hull') if i == 0 else None,
            showlegend=(i == 0),
            hoverinfo='skip'
        ))

    fig.update_layout(
        font=dict(color='gray'),
        plot_bgcolor='#242424',
        paper_bgcolor="#242424",
        margin=dict(l=20, r=20, t=50, b=20),
        uirevision='constant',
        legend=dict(
            groupclick="togglegroup",
            xanchor="right",
            yanchor="top",
            font=dict(size=11, color="white")
        ),
        annotations=[
            dict(
                text=plt_sub_title,
                showarrow=False,
                xref="paper",
                yref="paper",
                x=1,
                y=0,
                xanchor='right',
                yanchor='bottom',
                font=dict(size=12, color="gray")
            )
        ],
        height=540,
        scene=dict(
            xaxis=dict(range=[x_min, x_max], visible=False),
            yaxis=dict(range=[y_min, y_max], visible=False),
            zaxis=dict(range=[z_min, z_max], visible=False),
        ),
        title=dict(text=plt_title)
    )
    return fig