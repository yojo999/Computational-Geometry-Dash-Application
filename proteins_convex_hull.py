"""
Protein convex hull calculations and visualization for the computational geometry dashboard.

Defines ProteinsConvexHull class for computing convex hull features and
generating 3D charts for protein structures.
"""

import sys
import numpy as np
import pandas as pd
import plotly.io as pio
import plotly.graph_objects as go
from scipy.spatial import ConvexHull
from proteins_point_cloud import ProteinPointCloudData
from comp_geometry_helper import TriangleCalculations as tri_calc
from comp_geometry_helper import consecutive_among_three, get_flat_list

class ProteinsConvexHull:
    """
    Compute and visualize the convex hull of protein structures.

    Attributes:
        pdb_id (str): PDB file path or identifier.
        chain_id (str): Protein chain identifier.
        atomic_df (pd.DataFrame): Atomic dataframe.
        coords (np.ndarray): Coordinates array.
        aa_seq (list): Amino acid sequence.
        temp_factors (np.ndarray): Temperature factors.
        atom_type (str): Atom type ('calpha' or 'all_atom').
    """

    def __init__(self, pdb_id, chain_id, atomic_df, coords, aa_seq, temp_factors, atom_type):
        """
        Initialize ProteinsConvexHull.

        Args:
            pdb_id (str): PDB file path or identifier.
            chain_id (str): Protein chain identifier.
            atomic_df (pd.DataFrame): Atomic dataframe.
            coords (np.ndarray): Coordinates array.
            aa_seq (list): Amino acid sequence.
            temp_factors (np.ndarray): Temperature factors.
            atom_type (str): Atom type ('calpha' or 'all_atom').
        """
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.atomic_df = atomic_df
        self.coords = coords
        self.aa_seq = aa_seq
        self.temp_factors = temp_factors
        self.atom_type = atom_type

    def compute_convex_hull(self):
        """
        Compute convex hull features and save them to a CSV file.

        Calculates area, volume, triangle features, temperature factors,
        and saves results as a CSV file.
        """
        try:
            identifier = f"{self.pdb_id}_{self.chain_id}_cx_hull_{self.atom_type}"
            out_fname = f"{identifier}.csv"

            cx_hull = ConvexHull(self.coords, qhull_options='Qc')
            cx_hull_area = cx_hull.area
            cx_hull_volume = cx_hull.volume

            simplices = cx_hull.simplices
            simplices_number_list = ["-".join(map(str, row)) for row in simplices]
            simplices_coords = self.coords[simplices]
            simplices_seq_list = ['-'.join([self.aa_seq[i] for i in row]) for row in simplices]

            print(f"\t{self.pdb_id} chain {self.chain_id} contains {len(simplices)} simplices: {self.atom_type}")
            consecutive_list = [consecutive_among_three(list(row)) for row in simplices]

            # Temperature Factor for the Simplex residues
            tf_tri = self.temp_factors[simplices]
            tf_1_list = get_flat_list([row[0] for row in tf_tri])
            tf_2_list = get_flat_list([row[1] for row in tf_tri])
            tf_3_list = get_flat_list([row[2] for row in tf_tri])
            sum_tf_list = [round(np.sum(row), 2) for row in tf_tri]
            avg_tf_list = [round(np.mean(row), 2) for row in tf_tri]

            # Triangle calculations
            side_a_list, side_b_list, side_c_list, t_area_list = [], [], [], []
            for sim in simplices_coords:
                tc = tri_calc(sim)
                aa, bb, cc, ar = tc.triangle_properties()
                side_a_list.append(round(aa, 3))
                side_b_list.append(round(bb, 3))
                side_c_list.append(round(cc, 3))
                t_area_list.append(round(ar, 3))

            # Create and save DataFrame
            req_columns = list(zip(
                simplices_seq_list, simplices_number_list, consecutive_list,
                side_a_list, side_b_list, side_c_list, t_area_list,
                tf_1_list, tf_2_list, tf_3_list, sum_tf_list, avg_tf_list
            ))
            req_headers = [
                "amino_acid_residues", "simplex_indices", "simplex_type",
                "side_a", "side_b", "side_c", "triangle_area",
                "tf_1", "tf_2", "tf_3", "sum_tf", "avg_tf"
            ]
            cx_hull_df = pd.DataFrame(req_columns, columns=req_headers)
            cx_hull_df["cx_hull_area"] = cx_hull_area
            cx_hull_df["cx_hull_volume"] = cx_hull_volume
            cx_hull_df["pdb_identifier"] = identifier
            cx_hull_df_cols = list(cx_hull_df.columns)
            cx_hull_df_cols = [cx_hull_df_cols[-1]] + cx_hull_df_cols[:-1]
            cx_hull_df = cx_hull_df[cx_hull_df_cols]
            print(f"\tOutput file: {out_fname}")
            cx_hull_df.to_csv(out_fname, index=False)
        except Exception as e:
            print("\tUnknown exception occurred.", file=sys.stderr)

    def chart_convex_hull(self):
        """
        Generate and save a 3D Plotly chart of the protein convex hull.

        Saves the chart as an HTML file.
        """
        try:
            identifier = f"{self.pdb_id}_{self.chain_id}_cx_hull_{self.atom_type}"
            out_html = f"{identifier}.html"
            cx_hull = ConvexHull(self.coords, qhull_options='Qc')
            simplices_coords = self.coords[cx_hull.simplices]

            pio.templates.default = "plotly_dark"
            fig = go.Figure()

            fig.add_trace(go.Scatter3d(
                mode='markers',
                x=self.atomic_df['x_coord'],
                y=self.atomic_df['y_coord'],
                z=self.atomic_df['z_coord'],
                name=f"{self.pdb_id[:-4]} chain {self.chain_id}",
                ids=self.aa_seq,
                text=self.aa_seq,
                showlegend=False,
                marker=dict(color='lime', size=6, symbol='circle')
            ))

            for pt in simplices_coords:
                fig.add_trace(go.Mesh3d(
                    x=pt[:, 0], y=pt[:, 1], z=pt[:, 2],
                    alphahull=-1, color='mintcream',
                    flatshading=True, opacity=0.3,
                    delaunayaxis='y', showlegend=False,
                    hoverinfo='skip'
                ))

            fig.update_layout(
                scene=dict(
                    xaxis=dict(visible=False),
                    yaxis=dict(visible=False),
                    zaxis=dict(visible=False)
                )
            )

            print(f"\tSaving plotly chart to html file: {out_html}")
            fig.write_html(out_html)

        except Exception as e:
            print("\tConvex Hull Chart creation failed! Unknown Exception!", e)

if __name__ == "__main__":
    pdb_id = "./pdb_files/4MMK.pdb"
    chain_id = "A"
    atom_type = "calpha"
    ppcd = ProteinPointCloudData(pdb_id, chain_id)
    ca_df = ppcd.get_calpha_info()
    ca_df, coords, aa_seq, temp_factors = ppcd.get_point_cloud(ca_df)

    pcxh = ProteinsConvexHull(pdb_id, chain_id, ca_df, coords, aa_seq, temp_factors, atom_type)
    pcxh.compute_convex_hull()
    pcxh.chart_convex_hull()

    atom_type = "all_atom"
    aa_df = ppcd.get_all_atom_info()
    aa_df, coords, aa_seq, temp_factors = ppcd.get_point_cloud(aa_df)

    pcxh = ProteinsConvexHull(pdb_id, chain_id, aa_df, coords, aa_seq, temp_factors, atom_type)
    pcxh.compute_convex_hull()
    pcxh.chart_convex_hull()