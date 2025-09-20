"""
PDB file parsing utilities for the computational geometry dashboard.
Includes functions for converting amino acid codes and extracting features from PDB files.
"""

import os
from biopandas.pdb import PandasPdb

def aa3_to_aa1(aa3):
    """
    Convert a three-letter amino acid code to a one-letter code.

    Args:
        aa3 (str): Three-letter amino acid code.

    Returns:
        str: One-letter amino acid code, or '?' if unknown.
    """
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z',
        'XAA': 'X', 'TER': '*'
    }
    aa3 = aa3.upper()
    return aa_dict.get(aa3, '?')

def biopandas_pdb_parser(pdb_id, pdb_file_loc):
    """
    Parse a Protein Data Bank (PDB) file and extract features using Biopandas.

    Args:
        pdb_id (str): PDB file name or identifier.
        pdb_file_loc (str): Directory containing the PDB file.

    Returns:
        tuple: (
            header,
            atomic_df_ca,
            atomic_df_aa,
            bfactor_ca,
            bfactor_aa,
            coords_ca,
            coords_aa,
            aa_seq_ca,
            aa_seq_aa,
            aa1_ca,
            aa1_aa
        )
    """
    pdb_path = os.path.join(pdb_file_loc, pdb_id)
    bio_pandas_parser = PandasPdb()
    atomic_df = bio_pandas_parser.read_pdb(pdb_path)
    header = atomic_df.header

    atomic_df = atomic_df.get_model(1)
    atomic_df = atomic_df.df['ATOM']
    atomic_df = atomic_df[atomic_df['chain_id'] == "A"]
    atomic_df = atomic_df[atomic_df['alt_loc'] != "B"]
    atomic_df_aa = atomic_df.copy()
    atomic_df_ca = atomic_df_aa[atomic_df_aa['atom_name'] == "CA"].copy()

    atomic_df_ca.loc[:, 'label'] = (
        atomic_df_ca['residue_number'].astype(str) + '-' +
        atomic_df_ca['residue_name'].astype(str) + '-' +
        atomic_df_ca['atom_name'].astype(str) + '-' +
        atomic_df_ca['atom_number'].astype(str)
    )

    atomic_df_aa.loc[:, 'label'] = (
        atomic_df_aa['residue_number'].astype(str) + '-' +
        atomic_df_aa['residue_name'].astype(str) + '-' +
        atomic_df_aa['atom_name'].astype(str) + '-' +
        atomic_df_aa['atom_number'].astype(str)
    )

    coords_ca = atomic_df_ca[['x_coord', 'y_coord', 'z_coord']].to_numpy()
    coords_aa = atomic_df_aa[['x_coord', 'y_coord', 'z_coord']].to_numpy()

    bfactor_ca = atomic_df_ca[['b_factor']].to_numpy()
    bfactor_aa = atomic_df_aa[['b_factor']].to_numpy()

    aa_seq_ca = list(atomic_df_ca['residue_name'])
    aa_seq_aa = list(atomic_df_aa['residue_name'])

    aa1_ca = "".join([aa3_to_aa1(i) for i in aa_seq_ca])
    aa1_aa = "".join([aa3_to_aa1(i) for i in aa_seq_aa])

    return (
        header, atomic_df_ca, atomic_df_aa, bfactor_ca, bfactor_aa,
        coords_ca, coords_aa, aa_seq_ca, aa_seq_aa, aa1_ca, aa1_aa
    )