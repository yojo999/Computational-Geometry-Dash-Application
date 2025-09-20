"""
Protein point cloud extraction utilities for the computational geometry dashboard.
Defines ProteinPointCloudData class for parsing PDB files and extracting
atomic and point cloud data for protein structures.
"""

import sys
from biopandas.pdb import PandasPdb
from comp_geometry_helper import get_flat_list

class ProteinPointCloudData:
    """
    Handles parsing of PDB files and extraction of atomic and point cloud data.

    Attributes:
        pdb_id (str): Path to the PDB file.
        chain_id (str): Protein chain identifier.
        atom_type (str): Atom type (optional).
        atomic_df (PandasPdb): Parsed atomic dataframe.
    """

    def __init__(self, pdb_id, chain_id):
        """
        Initialize ProteinPointCloudData.

        Args:
            pdb_id (str): Path to the PDB file.
            chain_id (str): Protein chain identifier.
        """
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.atom_type = None
        self.atomic_df = None
        self._parse_file()

    def _parse_file(self):
        """
        Parse the PDB file using BioPandas.

        Handles exceptions and sets self.atomic_df.
        """
        bio_pandas_parser = PandasPdb()
        try:
            print(f"Attempting to parse the PDB file {self.pdb_id}.")
            self.atomic_df = bio_pandas_parser.read_pdb(self.pdb_id)
            print("Parsing successful!")
        except FileNotFoundError:
            print("File not found.", file=sys.stderr)
            self.atomic_df = None
        except Exception as e:
            print("Unknown exception occurred. Check validity of input file.", file=sys.stderr)
            self.atomic_df = None

    def get_calpha_info(self):
        """
        Extract C-alpha carbon atom information from the atomic dataframe.

        Returns:
            pd.DataFrame: DataFrame containing C-alpha atoms for the specified chain.
        """
        if self.atomic_df is None:
            print("\tBioPandas Object not found!", file=sys.stderr)
            return None
        try:
            ca_atomic_df = self.atomic_df.get_model(1)
            ca_atomic_df = ca_atomic_df.df['ATOM']
            ca_atomic_df = ca_atomic_df[ca_atomic_df['chain_id'] == self.chain_id]
            ca_atomic_df = ca_atomic_df[ca_atomic_df['alt_loc'] != "B"]
            ca_atomic_df = ca_atomic_df[ca_atomic_df['atom_name'] == "CA"]
            return ca_atomic_df
        except Exception as e:
            print("\tUnknown exception occurred. Check validity of input file.", file=sys.stderr)
            return None

    def get_all_atom_info(self):
        """
        Extract all atom information from the atomic dataframe.

        Returns:
            pd.DataFrame: DataFrame containing all atoms for the specified chain.
        """
        if self.atomic_df is None:
            print("\tBioPandas DataFrame Object not found!", file=sys.stderr)
            return None
        try:
            all_atomic_df = self.atomic_df.get_model(1)
            all_atomic_df = all_atomic_df.df['ATOM']
            all_atomic_df = all_atomic_df[all_atomic_df['chain_id'] == self.chain_id]
            all_atomic_df = all_atomic_df[all_atomic_df['alt_loc'] != "B"]
            return all_atomic_df
        except Exception as e:
            print("\tUnknown exception occurred.", file=sys.stderr)
            return None

    def get_point_cloud(self, atomic_df):
        """
        Extract point cloud data (coordinates, sequence, temperature factors) from atomic dataframe.

        Args:
            atomic_df (pd.DataFrame): Atomic dataframe to extract from.

        Returns:
            tuple: (atomic_df, coords, aa_seq, temp_factors)
                atomic_df (pd.DataFrame): Input atomic dataframe.
                coords (np.ndarray): Nx3 array of coordinates.
                aa_seq (list): List of amino acid residue names.
                temp_factors (np.ndarray): Nx1 array of temperature factors.
        """
        if atomic_df is None:
            print("\tBioPandas DataFrame Object not found!", file=sys.stderr)
            return None
        try:
            print(f"\tExtracting point cloud for {self.pdb_id} chain {self.chain_id}...")
            aa_seq = []
            coords = atomic_df[['x_coord', 'y_coord', 'z_coord']].to_numpy()
            temp_factors = atomic_df[['b_factor']].to_numpy()
            for res_id in atomic_df['residue_number'].unique():
                res_name = atomic_df[atomic_df['residue_number'] == res_id]['residue_name'].values
                aa_seq.append(list(res_name))
            aa_seq = get_flat_list(aa_seq)
            return atomic_df, coords, aa_seq, temp_factors
        except Exception as e:
            print("\tUnknown exception occurred.", file=sys.stderr)
            return None