# Package imports
from dash import Input, Output, callback
from dash_bio.utils import PdbParser

# Local imports
from src.utils.pdb_viewer import modify_style


@callback(
    Output('molecule-3d-viewer', 'modelData'),
    Input('pdb-file-dropdown', 'value'),
)
def update_molecule_viewer(selected_file):
    if selected_file:
        try:
            parser = PdbParser(selected_file)
            model_data = parser.mol3d_data()

            return model_data

        except Exception as e:
            print(f"Error processing PDB file: {e}")

            return {'atoms': [], 'bonds': []}

    return {'atoms': [], 'bonds': []}

@callback(
    Output('molecule-3d-viewer', 'styles'),
    [Input('pdb-file-dropdown', 'value'),
     ],
)
def style_molecule(selected_file):
    if selected_file:
        apply_style = modify_style(selected_file)
        return apply_style if apply_style is not None else {}
    return {}
