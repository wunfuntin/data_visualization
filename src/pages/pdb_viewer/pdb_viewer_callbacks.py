# Package imports
import dash_bio.utils.ngl_parser as ngl_parser
from dash import Input, Output, callback
from dash_bio.utils import PdbParser
from dash.exceptions import PreventUpdate

# Local imports
from src.utils.pdb_viewer import modify_style, get_most_recent_directory


@callback(
    Output('molecule-3d-viewer', 'modelData'),
    Input('pdb-file-dropdown', 'value'),
)
def update_molecule_viewer(selected_file):
    print(selected_file)
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


@callback(
    Output('ngl-molecule-viewer', 'data'),
    Output('ngl-molecule-viewer', 'molStyles'),
    Input('ngl-file-dropdown', 'value'),
)
def update_ngl_viewer(selected_file):
    if selected_file is None:

        raise PreventUpdate

    parent_directory = ('/Users/rwalker/Documents/TT_lab/python/dash_prot_char_v1_1/'
                        'silent_test_data/extracted_pdbs')

    directory = get_most_recent_directory(parent_directory)
    if not directory.endswith('/'):
        directory += '/'

    selected_file_with_extension = selected_file.split('/')[-1]
    file_name = selected_file_with_extension.split('.')[0]

    mol_styles_dict = {
        'representations': ['cartoon'],
        'chosenAtomsColor': 'white',
        'chosenAtomsRadius': 1,
        'molSpacingXaxis': 100,
    }

    data_list = [ngl_parser.get_data(data_path=directory,
                                     pdb_id=file_name,
                                     color='red',
                                     reset_view=True,
                                     local=True,
                                     )
                 ]

    return data_list, mol_styles_dict

