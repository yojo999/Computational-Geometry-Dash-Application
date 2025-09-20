"""
Protein Delaunay tessellation calculations and visualization for the computational geometry dashboard.
Defines ProteinsDelaunayTessellation class for computing Delaunay features and
generating 3D charts for protein structures.
"""

import sys
import numpy as np
import pandas as pd
import plotly.io as pio
import plotly.graph_objects as go
from scipy.spatial import Delaunay, Voronoi
from proteins_point_cloud import ProteinPointCloudData
from comp_geometry_helper import TetrahedronCalculations as tet_calc
from comp_geometry_helper import consecutive_among_four, get_flat_list

class ProteinsDelaunayTessellation:
    """
    Compute and visualize the Delaunay tessellation of protein structures.

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
        Initialize ProteinsDelaunayTessellation.

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

    def compute_delaunay(self):
        """
        Compute Delaunay tessellation features and save them to a CSV file.

        Calculates tetrahedron features, temperature factors, and saves results as a CSV file.
        """
        try:
            identifier = f"{self.pdb_id}_{self.chain_id}_della_{self.atom_type}"
            out_fname = f"{identifier}.csv"

            della = Delaunay(self.coords, qhull_options='Qc', furthest_site=False)
            vorro = Voronoi(self.coords, qhull_options='Qc', furthest_site=False)
            num_vorro_regions = len(vorro.regions)
            print(f"\tGenerating Delaunay Tessellation for {self.pdb_id} chain {self.chain_id}...")
            print(f"\tNumber of input points : {len(della.points)}")

            simplices = della.simplices
            print("\tPerforming Calculations...")
            print(f"\t{self.pdb_id} chain {self.chain_id} contains {len(simplices)} simplices: {self.atom_type}")

            # Obtain simplex indices and amino acid residues
            simplices_number_list = ["-".join(map(str, row)) for row in simplices]
            simplices_coords = self.coords[della.simplices]
            simplices_seq_list = ['-'.join([self.aa_seq[i] for i in row]) for row in simplices]

            # Get simplex type based on number of consecutive residues
            simplices_list = [list(i) for i in simplices]
            consecutive_list_four = [consecutive_among_four(i) for i in simplices_list]

            # Temperature Factor for the Simplex residues
            tf_quad = self.temp_factors[della.simplices]
            tf_1_list = get_flat_list([list(i[0]) for i in tf_quad])
            tf_2_list = get_flat_list([list(i[1]) for i in tf_quad])
            tf_3_list = get_flat_list([list(i[2]) for i in tf_quad])
            tf_4_list = get_flat_list([list(i[3]) for i in tf_quad])
            sum_tf_list = [round(np.sum(i), 2) for i in tf_quad]
            avg_tf_list = [round(np.mean(i), 2) for i in tf_quad]

            # Calculate tetrahedron features
            Tab, Tac, Tad = [], [], []
            Tbc, Tbd, Tcd = [], [], []
            Tsuml, Tavgl, Tvol = [], [], []
            Tarea1, Tarea2, Tarea3, Tarea4 = [], [], [], []
            T_surf_area = []
            O1, O2, O3, O4 = [], [], [], []
            for sim in simplices_coords:
                tetc = tet_calc(sim)
                ab, ac, ad, bc, bd, cd, suml, avgl = tetc.tetra_edge_lengths()
                Tab.append(round(ab, 3))
                Tac.append(round(ac, 3))
                Tad.append(round(ad, 3))
                Tbc.append(round(bc, 3))
                Tbd.append(round(bd, 3))
                Tcd.append(round(cd, 3))
                Tsuml.append(round(suml, 3))
                Tavgl.append(round(avgl, 3))

                tvol = tetc.tetra_volume()
                Tvol.append(round(tvol, 3))

                t1, t2, t3, t4, t_area = tetc.tetra_surface_area()
                Tarea1.append(round(t1, 3))
                Tarea2.append(round(t2, 3))
                Tarea3.append(round(t3, 3))
                Tarea4.append(round(t4, 3))
                T_surf_area.append(round(t_area, 3))

                o1, o2, o3, o4 = tetc.tetra_solid_angles()
                O1.append(round(o1, 3))
                O2.append(round(o2, 3))
                O3.append(round(o3, 3))
                O4.append(round(o4, 3))

            req_columns = list(zip(
                simplices_seq_list, simplices_number_list, consecutive_list_four,
                Tab, Tac, Tad, Tbc, Tbd, Tcd, Tsuml, Tavgl, Tvol,
                Tarea1, Tarea2, Tarea3, Tarea4, T_surf_area, O1, O2, O3, O4
            ))

            req_headers = [
                "amino_acid_residues", "simplex_indices", "simplex_type",
                "L1", "L2", "L3", "L4", "L5", "L6", "SumL", "AvgL", "Tvol",
                "T1_area", "T2_area", "T3_area", "T4_area", "Tsurface_area",
                "O1", "O2", "O3", "O4"
            ]

            della_df = pd.DataFrame(req_columns, columns=req_headers)
            della_df["num_voronoi"] = num_vorro_regions
            della_df["pdb_identifier"] = identifier

            della_df_cols = list(della_df.columns)
            della_df_cols = [della_df_cols[-1]] + della_df_cols[:-1]
            della_df = della_df[della_df_cols]

            della_df.to_csv(out_fname, index=False)
            print(f"\tOutput file: {out_fname}")
            print("\tDone!")
        except Exception as e:
            print("Unknown exception occurred.", file=sys.stderr)

    def chart_delaunay(self):
        """
        Generate and save a 3D Plotly chart of the protein Delaunay tessellation.

        Saves the chart as an HTML file.
        """
        try:
            identifier = f"{self.pdb_id}_{self.chain_id}_della_{self.atom_type}"
            out_html = f"{identifier}.html"
            della = Delaunay(self.coords, qhull_options='Qc')
            simplices_coords = self.coords[della.simplices]

            pio.templates.default = "plotly_dark"
            fig = go.Figure()

            fig.add_trace(go.Scatter3d(
                mode='markers',
                x=self.atomic_df['x_coord'],
                y=self.atomic_df['y_coord'],
                z=self.atomic_df['z_coord'],
                name=f"{self.pdb_id[0:-4]} chain {self.chain_id}",
                ids=self.aa_seq,
                text=self.aa_seq,
                showlegend=True,
                marker=dict(
                    color='lime',
                    size=6,
                    symbol='circle'
                )
            ))

            for pt in simplices_coords:
                fig.add_trace(go.Mesh3d(
                    x=pt[:, 0],
                    y=pt[:, 1],
                    z=pt[:, 2],
                    alphahull=-1,
                    color='mintcream',
                    flatshading=False,
                    opacity=0.1,
                    delaunayaxis='y',
                    showlegend=False,
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
            print("\tDelaunay Tessellation Chart creation failed! Unknown Exception!", e)

if __name__ == "__main__":
    pdb_id = "./pdb_files/4MMK.pdb"
    chain_id = "A"
    atom_type = "calpha"
    ppcd = ProteinPointCloudData(pdb_id, chain_id)
    ca_df = ppcd.get_calpha_info()
    ca_df, coords, aa_seq, temp_factors = ppcd.get_point_cloud(ca_df)
    pdt = ProteinsDelaunayTessellation(pdb_id, chain_id, ca_df, coords, aa_seq, temp_factors, atom_type)
    pdt.compute_delaunay()
    pdt.chart_delaunay()

    pdb_id = "./pdb_files/4MMK.pdb"
    chain_id = "A"
    atom_type = "all_atom"
    ppcd = ProteinPointCloudData(pdb_id, chain_id)
    aa_df = ppcd.get_all_atom_info()
    aa_df, coords, aa_seq, temp_factors = ppcd.get_point_cloud(aa_df)
    pdt = ProteinsDelaunayTessellation(pdb_id, chain_id, aa_df, coords, aa_seq, temp_factors, atom_type)
    pdt.compute_delaunay()