import os
import dash
import plotly
import dash_bio
import numpy as np
import pandas as pd
import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go
from dash import callback_context
from dash_iconify import DashIconify
from dash_bio.utils import PdbParser
import dash_mantine_components as dmc
from dash_bio.utils import protein_reader
from dash.exceptions import PreventUpdate
import dash_bio.utils.ngl_parser as ngl_parser
from dash import (dcc, html, Dash, _dash_renderer, dash_table,
    Input, Output, State, callback, clientside_callback)

from descriptions import (desc_1A7I, desc_1AAZ, desc_1A32,
    desc_1AG6, desc_1AOJ, desc_1CRN, desc_1AL1, desc_about)
from pdb_parser import biopandas_pdb_parser, aa3_to_aa1
from biopandas.pdb import PandasPdb
from sequence_chart import custom_sequence_chart
from computational_geometry_chart import cg_chart
from mol3d_styles import create_mol3d_style
from custom_rama_chart import custom_ramachandran_chart

_dash_renderer._set_react_version("18.2.0")
app = dash.Dash(__name__)
app.title = "Computational Geometry of three dimensional Protein Structures"

CARD_BG = "#2e2e2e"
APP_BG = "#242424"
TEXT_COLOR = "white"

modal_text_mapping = {
    "1A7I.pdb": desc_1A7I,
    "1A32.pdb": desc_1A32,
    "1AAZ.pdb": desc_1AAZ,
    "1AG6.pdb": desc_1AG6,
    "1AL1.pdb": desc_1AL1,
    "1AOJ.pdb": desc_1AOJ,
    "1CRN.pdb": desc_1CRN,
}

pdb_options = [
    {"label": "1A7I : A", "value": "1A7I.pdb"},
    {"label": "1A32 : A", "value": "1A32.pdb"},
    {"label": "1AAZ : A", "value": "1AAZ.pdb"},
    {"label": "1AG6 : A", "value": "1AG6.pdb"},
    {"label": "1AL1 : A", "value": "1AL1.pdb"},
    {"label": "1AOJ : A", "value": "1AOJ.pdb"},
    {"label": "1CRN : A", "value": "1CRN.pdb"},
]

csv_options = [
    "PDB file", 
    "C-Alpha Convex Hull Dataset", 
    "All Atom Convex Hull Dataset", 
    "C-Alpha Delaunay Tessellation Dataset", 
    "All Atom Delaunay Tessellation Dataset"
]

app.layout = dmc.MantineProvider(
    theme={"colorScheme": "dark"},
    children=html.Div([
        # Header
        html.Div([
            dmc.Drawer(
                id="info-drawer",
                padding="xs",
                size="100vh",
                opened=False,
                children=[
                    dmc.ScrollArea(
                        dcc.Markdown(desc_about, dangerously_allow_html=True),
                        style={"height": "90vh"}
                    )]
                ),
            ]),
        dmc.Flex(
            justify="space-between", align="center", wrap="wrap", gap="md",
            direction="row", p="xs",
            children=[
                dmc.Title("Computational Geometry of three dimensional Protein Structures", order=2, c="dodgerblue", textWrap='balance'),
                dmc.Group(
                    children=[
                        dmc.Button("Reset view", size="compact-xs", id="reset-btn", color="dodgerblue"),
                        html.A(
                            dmc.Button("GitHub", color="dodgerblue", size="compact-xs", id="github-btn"),
                            href="https://github.com/yojo999/Computational-Geometry-Dash-Application",
                            target="_blank",
                            style={"textDecoration": "none"}
                        ),
                        dmc.Button("About", size="compact-xs", id="about-btn", color="dodgerblue"),
                    ]
                )
            ]
        ),
        # Grids
        dmc.Grid(
            gutter="xs", p="xs",
            children=[
                dmc.GridCol(span=2,
                    style={"display": "flex", "justifyContent": "center"},
                    children=[
                        dmc.Group([
                            dmc.Card(
                                shadow="xs",
                                padding="xs",
                                radius="xs",
                                withBorder=False,
                                style={"backgroundColor": CARD_BG},
                                children=[
                                    dmc.Select(
                                        id="pdb-selector",
                                        value="1A7I.pdb",
                                        data=pdb_options,
                                        w=175,
                                        mb=10),
                            ]),
                        dmc.ActionIcon(
                            id="modal-btn",
                            variant="light",
                            color="dodgerblue",
                            children=DashIconify(icon="mdi:information-outline", width=20),
                            style={"cursor": "pointer"},
                        ),
                        dmc.Modal(
                            id="pdb-modal",
                            children=html.Div(
                                dcc.Markdown(
                                    id="pdb-modal-text",
                                    dangerously_allow_html=True)),
                            opened=False),
                    ])]
                    ),
                dcc.Store(id='seq-store'),
                dmc.GridCol(
                    span=10, 
                    children=dmc.Card(
                        shadow="xs",
                        padding="xs",
                        radius="xs",
                        withBorder=False,
                        style={"backgroundColor": CARD_BG},
                        children=dcc.Graph(
                                id="seq-viewer",
                                config={
                                    'scrollZoom': False,
                                    'displayModeBar': False,
                                    'staticPlot': False
                                },
                                style={'height': '100px', 'padding': '10'})
                    )),
                ]),
        dmc.Grid(
            gutter="xs", p="xs",
            children=[
                dcc.Store(
                    id='mol-store',
                    data = {
                        'modelData': {'atoms': [], 'bonds': []},
                        'styles': [],
                        'selectedAtomIds': []}),
                dcc.Store(id='della-store'),
                dcc.Store(id='selection-store', data=[]),
                dcc.Store(id='rama-store'),
                dmc.GridCol(span=3, 
                    children=dmc.Card(
                        shadow="xs",
                        padding="xs",
                        radius="xs",
                        withBorder=False,
                        style={"backgroundColor": CARD_BG},
                        children=dmc.Group(
                            children=[
                                dmc.Text("Select Residue(s)", c='grey'),
                                dash_bio.Molecule3dViewer(
                                    id='mol-viewer',
                                    modelData={'atoms': [], 'bonds': []},
                                    styles=[],
                                    selectedAtomIds=[])
                                ])
                    ),),
                dmc.GridCol(span=5, 
                    children=dmc.Card(
                        shadow="xs",
                        padding="xs",
                        radius="xs",
                        withBorder=False,
                        style={"backgroundColor": CARD_BG},
                        children=dcc.Graph(id="della-chart")
                        )),
                dmc.GridCol(span=4, 
                    children=dmc.Card(
                        shadow="xs",
                        padding="xs",
                        radius="xs",
                        withBorder=False,
                        style={"backgroundColor": CARD_BG},
                        children=dcc.Graph(id="rama-chart"))
                        ),
                ]),
        dmc.Grid(
            gutter="xs", p="xs",
            children=[
                dmc.GridCol(
                    span=12,
                    children=dmc.Card(
                        shadow="sm", padding="xs", radius="xs", withBorder=False,
                        style={"backgroundColor": CARD_BG},
                        children=dmc.Flex(
                            justify="space-between", align="center",
                            direction="row", wrap="wrap", gap="md",
                            children=[
                                dmc.Text("Selected Residues:", id="show-selection", c="white"),
                                dmc.Select(
                                    id="csv-selector",
                                    data=csv_options,
                                    placeholder="Download Options",
                                    comboboxProps={"dropdownPadding": 0},
                                    clearable=True,
                                    value=None),
                            dcc.Download(id="csv-download"),
                            dcc.Store(id='table-store'),
                            ]
                        )
                    )
                ),
                dmc.GridCol(
                    span=12,
                    children=dmc.Card(
                        shadow="sm", padding="xs", radius="xs", withBorder=False,
                        style={"backgroundColor": CARD_BG},
                        children=dmc.Flex(
                            justify="space-between", align="center",
                            direction="row", wrap="wrap", gap="md",
                            children=[dash_table.DataTable(
                            id='della-table',
                            columns=[],
                            data=[],
                            editable=False,
                            filter_action='native',
                            sort_action='native',
                            sort_mode='multi',
                            column_selectable=False,
                            row_selectable=False,
                            row_deletable=False,
                            page_action='native',
                            page_current=0,
                            page_size=10,
                            style_header={
                                'textAlign': 'center',
                                'backgroundColor': '#242424',
                                'color': 'white'},
                            style_data={
                                'textAlign': 'center',
                                'backgroundColor': '#242424',
                                'color': 'white'},
                            style_filter={
                                "backgroundColor": "#242424", 
                                "color": "white", 
                                },
                            style_data_conditional=[
                                {"if": {"state": "active"},
                                "color": "white",
                                "backgroundColor": "#242424"
                                }
                                ],
                            ),
                            html.Div(id="test-out",
                                style={"display": "none"})
                            ]
                        )
                    )
                )
            ]
        )
    ], style={"backgroundColor": APP_BG, "minHeight": "100vh", "color": TEXT_COLOR})
)

@callback(
    Output("info-drawer", "opened"),
    Output("info-drawer", "position"),
    Input("about-btn", "n_clicks"),
    State("info-drawer", "opened"),
    prevent_initial_call=True,
)
def toggle_drawer(n_clicks, position):
    """
    Opens the info drawer when the About button is clicked.

    Args:
        n_clicks (int): Number of clicks on the About button.
        position (str): Current drawer position.

    Returns:
        tuple: (opened state, position)
    """
    return True, position

@app.callback(
    Output("pdb-modal", "opened"),
    Input("modal-btn", "n_clicks"),
    prevent_initial_call=True
)
def open_modal(n_clicks):
    """
    Opens the PDB modal when the info icon is clicked.

    Args:
        n_clicks (int): Number of clicks on the info icon.

    Returns:
        bool: True to open the modal.
    """
    return True

@app.callback(
    Output("pdb-modal-text", "children"),
    Input("pdb-selector", "value"),
)
def update_modal_text(pdb_id):
    """
    Updates the modal text based on the selected PDB file.

    Args:
        pdb_id (str): Selected PDB file name.

    Returns:
        str: Description text for the selected PDB.
    """
    print("from update_modal: ", pdb_id)
    return modal_text_mapping.get(pdb_id)

@app.callback(
    Output('mol-store', 'data'),
    Output('seq-store', 'data'),
    Output('table-store', 'data'),
    Output('della-store', 'data'),
    Output('rama-store', 'data'),
    Input('pdb-selector', 'value'),
    Input('reset-btn', 'n_clicks')
)
def update_store(pdb_id, reset_clicks):
    """
    Loads and parses PDB data, updates stores for molecule viewer, sequence chart,
    Delaunay chart, Ramachandran chart, and table.

    Args:
        pdb_id (str): Selected PDB file name.
        reset_clicks (int): Number of clicks on the reset button.

    Returns:
        tuple: Data for molecule viewer, sequence chart, table, Delaunay chart, Ramachandran chart.
    """
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "pdb_files")
    pdb_path = os.path.join(data_dir, pdb_id)
    print("pdb_path", pdb_path)

    triggered_id = callback_context.triggered[0]['prop_id'].split('.')[0]
    if not pdb_id:
        raise dash.exceptions.PreventUpdate

    print(f"[INFO] Reloading: {pdb_id} (triggered by {triggered_id})")

    hed, df_ca, df_aa, bf_ca, bf_aa, xyz_ca, xyz_aa, aa_seq_ca, aa_seq_aa, aa1_ca, aa1_aa = biopandas_pdb_parser(pdb_id, data_dir)
    seq_chart = custom_sequence_chart(aa1_ca, hed, highlight_indices=None)

    parser = PdbParser(pdb_path)
    data = parser.mol3d_data()
    styles = create_mol3d_style(
        data['atoms'],
        visualization_type='cartoon',
        color_element='chain',
    )
    mol_store_data = {
        'modelData': data,
        'styles': styles,
        'selectedAtomIds': []
    }

    della_chart = cg_chart(pdb_id, ".pdb", "A", data_dir, [])
    rama_chart = custom_ramachandran_chart(pdb_id, pdb_path, [])

    della_cols = ["amino_acid_residues","simplex_indices","simplex_type","L1","L2","L3",
        "L4","L5","L6","SumL","AvgL","Tvol","T1_area","T2_area","T3_area","T4_area",
        "Tsurface_area","O1","O2","O3","O4"]
    csv_name = f"{pdb_id[0:4]}_A_della_calpha.csv" 
    csv_path = os.path.join(data_dir, csv_name)
    della_df = pd.read_csv(csv_path)
    table_data = della_df.to_dict("records")
    return mol_store_data, seq_chart, table_data, della_chart, rama_chart

@app.callback(
    Output('mol-viewer', 'modelData'),
    Output('mol-viewer', 'styles'),
    Output('seq-viewer', 'figure'),
    Output('della-table', 'columns'),
    Output('della-table', 'data'),
    Output('della-chart', 'figure'),
    Output('show-selection', 'children'),
    Output('rama-chart', 'figure'),
    Input('mol-store', 'data'),
    Input('seq-store', 'data'),
    Input('table-store', 'data'),
    Input('della-store', 'data'),
    Input('mol-viewer', 'selectedAtomIds'),  
    Input('reset-btn', 'n_clicks'),
    Input('pdb-selector', 'value'),
    prevent_initial_call=True
)
def update_view(
    mol_data, seq_data, table_data, della_chart, 
    selected_atom_ids, reset_clicks, pdb_id):
    """
    Updates the main view components based on user interaction and store data.

    Args:
        mol_data (dict): Molecule viewer data.
        seq_data (dict): Sequence chart data.
        table_data (list): Table data.
        della_chart (plotly.graph_objs.Figure): Delaunay chart.
        selected_atom_ids (list): Selected atom IDs in molecule viewer.
        reset_clicks (int): Number of clicks on reset button.
        pdb_id (str): Selected PDB file name.

    Returns:
        tuple: Updated modelData, styles, sequence chart, table columns, table data,
               Delaunay chart, selection text, Ramachandran chart.
    """
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "pdb_files")
    pdb_path = os.path.join(data_dir, pdb_id)
    print("pdb_path", pdb_path)
    ctx = dash.callback_context
    if not ctx.triggered:
        raise dash.exceptions.PreventUpdate
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

    model_data = mol_data.get('modelData', {'atoms': [], 'bonds': []}) if mol_data else {}
    atom_styles = mol_data.get('styles', []) if mol_data else []
    sel_atms = mol_data.get('selectedAtomIds', []) if mol_data else []

    df_loaded = pd.DataFrame(table_data) if table_data else pd.DataFrame()
    table_columns = [{"name": i, "id": i, 'selectable': True} for i in df_loaded.columns[1: -1]] if not df_loaded.empty else []

    if not pdb_id:
        pdb_id = "default"

    if triggered_id == 'mol-viewer' and selected_atom_ids:
        res_index = [model_data['atoms'][aid]['residue_index'] for aid in selected_atom_ids]
        print(res_index)
        hed, df_ca, df_aa, bf_ca, bf_aa, xyz_ca, xyz_aa, aa_seq_ca, aa_seq_aa, aa1_ca, aa1_aa = biopandas_pdb_parser(pdb_id, data_dir)
        seq_chart = custom_sequence_chart(aa1_ca, hed, res_index)
        df = pd.DataFrame(table_data)
        filtered_df = df[df["simplex_indices"].apply(
            lambda x: any(str(val) in x.split("-") for val in res_index)
            )]
        table_data = filtered_df.to_dict("records")
        della_chart = cg_chart(pdb_id, ".pdb", "A", data_dir, res_index)
        rama_chart = custom_ramachandran_chart(pdb_id, pdb_path, res_index)
        residues = ", ".join(str(i) for i in res_index)
        return model_data, atom_styles, seq_chart, table_columns, table_data, della_chart, f"Selected Residues : {residues}",  rama_chart

    else:
        model_data = mol_data.get('modelData', {'atoms': [], 'bonds': []}) if mol_data else {}
        atom_styles = mol_data.get('styles', []) if mol_data else []
        hed, df_ca, df_aa, bf_ca, bf_aa, xyz_ca, xyz_aa, aa_seq_ca, aa_seq_aa, aa1_ca, aa1_aa = biopandas_pdb_parser(pdb_id, data_dir)
        seq_chart = custom_sequence_chart(aa1_ca, hed, highlight_indices=None)
        rama_chart = custom_ramachandran_chart(pdb_id, pdb_path, [])
    return model_data, atom_styles, seq_chart, table_columns, table_data, della_chart, [], rama_chart

@app.callback(
    Output('mol-viewer', 'selectedAtomIds'),
    Input('reset-btn', 'n_clicks'),
    prevent_initial_call=True
)
def clear_selection(n_clicks):
    """
    Clears the selected atoms in the molecule viewer.

    Args:
        n_clicks (int): Number of clicks on the reset button.

    Returns:
        list: Empty list to clear selection.
    """
    return []

@app.callback(
    Output("csv-download", "data"),
    Output("csv-selector", "value"),
    Input("csv-selector", "value"),
    Input('pdb-selector', 'value'),
    prevent_initial_call=True
)
def download_csv_and_reset(selected_file, pdb_id):
    """
    Handles CSV and PDB file download based on user selection and resets dropdown.

    Args:
        selected_file (str): Selected file type for download.
        pdb_id (str): Selected PDB file name.

    Returns:
        tuple: Download data, None to reset dropdown.
    """
    if selected_file and pdb_id:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(base_dir, "pdb_files")
        pdb_dl_file = f"{pdb_id}"
        cx_calpha_csv = f"{pdb_id[0:4]}_A_cx_hull_calpha.csv"
        cx_all_atom_csv = f"{pdb_id[0:4]}_A_cx_hull_all_atom.csv"
        della_calpha_csv = f"{pdb_id[0:4]}_A_della_calpha.csv"
        della_all_atom_csv = f"{pdb_id[0:4]}_A_della_all_atom.csv"
        match selected_file:
            case "PDB file":
                pdb_fpath = os.path.join(data_dir, pdb_id)
                print("Selected PDB file for download", pdb_fpath)
                return dcc.send_file(pdb_fpath), None
            case "C-Alpha Convex Hull Dataset":
                fpath = os.path.join(data_dir, cx_calpha_csv)
                csv_df = pd.read_csv(fpath)
                return dcc.send_data_frame(csv_df.to_csv, f"{cx_calpha_csv}", index=False), None
            case "All Atom Convex Hull Dataset":
                fpath = os.path.join(data_dir, cx_all_atom_csv)
                csv_df = pd.read_csv(fpath)
                return dcc.send_data_frame(csv_df.to_csv, f"{cx_all_atom_csv}", index=False), None
            case "C-Alpha Delaunay Tessellation Dataset":
                fpath = os.path.join(data_dir, della_calpha_csv)
                csv_df = pd.read_csv(fpath)
                return dcc.send_data_frame(csv_df.to_csv, f"{della_calpha_csv}", index=False), None
            case "All Atom Delaunay Tessellation Dataset":
                fpath = os.path.join(data_dir, della_all_atom_csv)
                csv_df = pd.read_csv(fpath)
                return dcc.send_data_frame(csv_df.to_csv, f"{della_all_atom_csv}", index=False), None
            case _:
                print("Selected file not found!")
                return None
    return dash.no_update, dash.no_update

if __name__ == "__main__":
    app.run(debug=True)