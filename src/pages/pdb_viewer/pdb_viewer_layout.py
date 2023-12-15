# Package imports
import dash
import dash_bio as dashbio
import os

from dash import html, dcc
from dash_bootstrap_templates import load_figure_template

# Local imports
from src.utils.pdb_viewer import get_most_recent_directory
from . import pdb_viewer_callbacks


dash.register_page(
    __name__,
    path='/molview',
    title='molview',
)


load_figure_template('CYBORG')

parent_directory = '/Users/rwalker/Documents/TT_lab/python/dash_prot_char_v1_1/silent_test_data/extracted_pdbs'

directory = get_most_recent_directory(parent_directory)

file_names = [file for file in os.listdir(directory) if file.endswith('.pdb')]

# Create dropdown options
pdb_options = [{'label': name, 'value': os.path.join(directory, name)} for name in file_names]


layout = html.Div([
    html.H1('3D Molecular Structure Viewer'),

    dcc.Dropdown(
        id='pdb-file-dropdown',
        options=pdb_options,
        value=None,
        placeholder='Select a PDB to view'
    ),

    dashbio.Molecule3dViewer(
        id='molecule-3d-viewer',
        modelData={
            'atoms': [],
            'bonds': [],
        },
        styles={},
    ),
    html.H1('NGL Molecular Structure Viewer'),

    dcc.Dropdown(
        id='ngl-file-dropdown',
        options=pdb_options,
        value=None,
        placeholder='Select a NGL model to view'
    ),

    dashbio.NglMoleculeViewer(
        id='ngl-molecule-viewer',
        height='1500px',
        width='1500px',
    ),

])

