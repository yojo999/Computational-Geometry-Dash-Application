"""
Visualization helper functions for the computational geometry dashboard.
Includes color maps and utilities for extracting edge coordinates from
tetrahedra and triangles, and for selecting arrays by residue indices.
"""

import itertools
import numpy as np

atom_cmap = {
    'C': 'limegreen', 'CA': 'lime', 'CB': 'antiquewhite',
    'CD': 'floralwhite', 'CD1': 'floralwhite', 'CD2': 'floralwhite',
    'CE1': 'ivory', 'CE2': 'ivory', 'CG': 'navajowhite', 'CG1': 'navajowhite',
    'CG2': 'navajowhite', 'CZ': 'navajowhite', 'N': 'limegreen', 'ND2': 'lime',
    'NE': 'lime', 'NH1': 'lime', 'NH2': 'lime', 'O': 'gold', 'OD1': 'lightyellow',
    'OD2': 'lightyellow', 'OE1': 'lightyellow', 'OE2': 'lightyellow',
    'OG': 'lightyellow', 'OG1': 'lightyellow', 'OH': 'cornsilk',
    'OXT': 'forestgreen', 'SG': 'green'
}

def arrays_of_selected(
    res_index_arr: np.ndarray,
    to_select_lst: list
) -> np.ndarray:
    """
    Select rows from an array where any element matches the selection list.

    Args:
        res_index_arr (np.ndarray): Array of residue indices (NxM).
        to_select_lst (list): List of indices to select.

    Returns:
        np.ndarray: Selected rows from res_index_arr.
    """
    bool_mask = np.isin(res_index_arr, to_select_lst).any(axis=1)
    req_tetras = res_index_arr[bool_mask]
    return req_tetras

def get_edges_tetra(
    tetra_coordinates: np.ndarray
) -> tuple[list, list, list]:
    """
    Extract edge coordinates from tetrahedra.

    Args:
        tetra_coordinates (np.ndarray): Array of tetrahedra coordinates (Nx4x3).

    Returns:
        tuple: Lists of x, y, z coordinates for all edges.
    """
    xline, yline, zline = [], [], []
    for simplex in tetra_coordinates:
        points = [simplex[0], simplex[1], simplex[2], simplex[3]]
        pairs = np.array(list(itertools.combinations(points, 2)))
        for pair in pairs:
            for coord in pair:
                xline.append(coord[0])
                yline.append(coord[1])
                zline.append(coord[2])
    return xline, yline, zline

def get_edges_tri(
    tri_coordinates: np.ndarray
) -> tuple[list, list, list]:
    """
    Extract edge coordinates from triangles.

    Args:
        tri_coordinates (np.ndarray): Array of triangle coordinates (Nx3x3).

    Returns:
        tuple: Lists of x, y, z coordinates for all edges.
    """
    xline, yline, zline = [], [], []
    for simplex in tri_coordinates:
        points = [simplex[0], simplex[1], simplex[2]]
        pairs = np.array(list(itertools.combinations(points, 2)))
        for pair in pairs:
            for coord in pair:
                xline.append(coord[0])
                yline.append(coord[1])
                zline.append(coord[2])
    return xline, yline, zline