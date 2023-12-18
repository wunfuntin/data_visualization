# Package imports
import dash
import dash_bio as dashbio
import os

from dash import html, dcc
from dash_bootstrap_templates import load_figure_template

# Local imports
from src.utils.pdb_viewer_functions import get_most_recent_directory
from . import pdb_viewer_callbacks


dash.register_page(
    __name__,
    path="/molview",
    title="molview",
)


load_figure_template("CYBORG")

layout = html.Div(
    [
        html.H1("3D Molecular Structure Viewer"),
        dcc.Input(
            id="parent-pdb-directory",
            placeholder="Enter file path to PDB directories",
        ),
        dcc.Dropdown(
            id="directory-dropdown",
            placeholder="Select the PDB directory to view",
        ),
        dcc.Dropdown(
            id="pdb-file-dropdown",
            value=None,
            placeholder="Select a PDB to view",
        ),
        dashbio.Molecule3dViewer(
            id="molecule-3d-viewer",
            modelData={
                "atoms": [],
                "bonds": [],
            },
            styles={},
        ),
        html.H1("NGL Molecular Structure Viewer"),
        dcc.Dropdown(
            id="pdb-file-dropdown",
            value=None,
            placeholder="Select a NGL model to view",
        ),
        dashbio.NglMoleculeViewer(
            id="ngl-molecule-viewer",
            height="1500px",
            width="1500px",
        ),
    ]
)
