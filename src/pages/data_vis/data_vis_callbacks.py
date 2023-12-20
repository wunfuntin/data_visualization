# Package imports
import base64
import json

import dash
from dash import dcc, html, callback, DiskcacheManager
from dash.dependencies import Input, Output, State
import diskcache
import io
import pandas as pd
import plotly.graph_objs as go
import re
import sys
from dash.exceptions import PreventUpdate

from utils.common_functions import (generate_hover_text, create_correlation_plot,
                                        filter_dataframe, create_3d_correlation_plot,
                                        filter_dataframe_3d_scatter)

from utils.pdb_management import get_top_pdbs_from_silent

cache = diskcache.Cache('./cache')
background_callback_manager = DiskcacheManager(cache)

# Callback function to process CSV file upload
@callback(
    Output('dataframe-json', 'children'),
    [Input('upload-csv', 'contents')]
)
def handle_csv_upload(contents):
    if contents is not None:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        try:
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))

            # Convert to JSON for intermediate storage
            return df.to_json(date_format='iso', orient='split')

        except Exception as e:
            return html.Div(['There was an error processing the CSV file.'])
        
# Callback for Inputing and storing file paths for Silent Files and output directories
@callback(
    [Output('storedSil', 'data'),
     Output('storedOut', 'data')],
    [Input('submitbutton_id', 'n_clicks')],
    [State('storedSil', 'data'),
     State('silent', 'value'),
     State('outDir', 'value'),
     State('storedOut', 'data')]
)
def update_dirs(n_clicks, storedSil, silent, outDir, storedOut ):
    if n_clicks is None:
        raise PreventUpdate
    if silent not in storedSil:
        storedSil.append(silent)
    if outDir not in storedOut:
        storedOut.append(outDir)
    return storedSil, storedOut
    
# @callback(
#     Output('xaxis-column-1', 'options'),
#     [Input('dataframe-json', 'children')]
# )
# def update_xaxis_column_options(csv_json):
#     if csv_json:
#         df = pd.read_json(io.StringIO(csv_json), orient='split')
#         return [{'label': col, 'value': col} for col in df.columns if col != 'description']
#     else:
#         return []  # Return empty options if no data is available
#
#
# @callback(
#     Output('yaxis-column-1', 'options'),
#     [Input('dataframe-json', 'children')]
# )
# def update_yaxis_column_options(csv_json):
#     if csv_json:
#         df = pd.read_json(io.StringIO(csv_json), orient='split')
#         return [{'label': col, 'value': col} for col in df.columns if col != 'description']
#     else:
#         return []  # Return empty options if no data is available
#
#
# @callback(
#     Output('xaxis-column-2', 'options'),
#     [Input('dataframe-json', 'children')]
# )
# def update_xaxis_column_2_options(csv_json):
#     if csv_json:
#         df = pd.read_json(io.StringIO(csv_json), orient='split')
#         return [{'label': col, 'value': col} for col in df.columns if col != 'description']
#     else:
#         return []  # Return empty options if no data is available
#
#
# @callback(
#     Output('yaxis-column-2', 'options'),
#     [Input('dataframe-json', 'children')]
# )
# def update_yaxis_column_2_options(csv_json):
#     if csv_json:
#         df = pd.read_json(io.StringIO(csv_json), orient='split')
#         return [{'label': col, 'value': col} for col in df.columns if col != 'description']
#     else:
#         return []  # Return empty options if no data is available
#
#
# Callbacks to update scatter and correlation plots based on dropdown selection
# Update dropdown menu
@callback(
    Output('dropdown', 'options'),
    [Input('dataframe-json', 'children')]
)
def update_dropdown_options(csv_json):
    if csv_json is not None:
        df = pd.read_json(io.StringIO(csv_json), orient='split')
        return [{'label': col, 'value': col} for col in df.columns if col != 'description']
    else:
        return []  # Return empty or default options if merged_json is None


@callback(
    [Output('xaxis-column-1', 'options'),
     Output('yaxis-column-1', 'options'),
     Output('xaxis-column-2', 'options'),
     Output('yaxis-column-2', 'options'),
     Output('xaxis-column-3d', 'options'),
     Output('yaxis-column-3d', 'options'),
     Output('zaxis-column-3d', 'options'),],
    [Input('dataframe-json', 'children',)]
)
def update_axis_column_options(csv_json):
    if csv_json:
        df = pd.read_json(io.StringIO(csv_json), orient='split')
        options = [{'label': col, 'value': col} for col in df.columns if col != 'description']
        return [options, options, options, options, options, options, options]
    else:
        return [[] for _ in range(7)]  # Return empty options for all dropdowns if no data


# Callbacks for updating plots
@callback(
    Output('scatter-plot', 'figure'),
    [Input('dataframe-json', 'children'),  # Use merged data as input
     Input('dropdown', 'value')]
)
def update_scatter_plot(csv_json, selected_column):
    if csv_json is not None and selected_column is not None:
        df = pd.read_json(io.StringIO(csv_json), orient='split')
        fig = go.Figure(data=[
            go.Scatter(
                x=df['description'],
                y=df[selected_column],
                mode='markers',
                marker=dict(color=df[selected_column], colorscale='Viridis'),
                text=df.apply(lambda row: generate_hover_text(row, df.columns), axis=1),
                hoverinfo='text'
            )
        ])
        fig.update_layout(
            title=f"Scatter Plot of {selected_column} vs Description",
            xaxis_title="Description",
            yaxis_title=selected_column
        )
        return fig
    else:
        return go.Figure()

@callback(
    Output(component_id='scatter-plot-click', component_property='children'),
    [Input(component_id='scatter-plot', component_property='clickData')]
)
def detect_click(click_data):
    if click_data is not None:
        point_data = click_data['points'][0]
        message = f'you clicked point {point_data}'
        print(message)

        # Extract the text field
        text_field = click_data['points'][0]['text']

        # Split the text field into individual data entries
        data_entries = text_field.split('<br>')

        # Find the entry for 'design_name' and extract its value
        design_name = None
        for entry in data_entries:
            if entry.startswith('design_name:'):
                design_name = entry.split(': ')[1].strip()  # Split by ': ' and strip whitespace
                break

        print(design_name)  # This will print the value of design_name


    else:
        print('no click detected.')


@callback(
    Output('radar-plot', 'figure'),
    [Input('dataframe-json', 'children'),
     Input('scatter-plot', 'hoverData')]
)
def update_radar_plot(csv_json, hoverData):
    if csv_json is not None and hoverData is not None:
        df = pd.read_json(io.StringIO(csv_json), orient='split')
        if hoverData and 'points' in hoverData and hoverData['points']:
            hover_index = hoverData['points'][0]['pointIndex']
            row = df.iloc[hover_index]
            radar_attributes = [
                'RMSD', 'int_area_to_len_ratio', 'hydrogen_bonds',
                'salt_bridges'
            ]
            radar_data = [row[attr] for attr in radar_attributes if attr in row]
            fig = go.Figure(data=go.Scatterpolar(
                r=radar_data,
                theta=radar_attributes,
            ))
            fig.update_layout(
                polar=dict(
                    radialaxis=dict(visible=True)
                ),
                showlegend=False
            )
            return fig
    return go.Figure()


@callback(
    Output('correlation-plot-1', 'figure'),
    [Input('dataframe-json', 'children'),
     Input('xaxis-column-1', 'value'),
     Input('yaxis-column-1', 'value'),
     Input('x-min-input-1', 'value'),
     Input('x-max-input-1', 'value'),
     Input('y-min-input-1', 'value'),
     Input('y-max-input-1', 'value')]
)
def update_correlation_plot_1(csv_json, x_col, y_col, x_min, x_max, y_min, y_max):
    if csv_json is not None and x_col is not None and y_col is not None:
        df = pd.read_json(io.StringIO(csv_json), orient='split')
        filtered_df = filter_dataframe(df, x_col, y_col, x_min, x_max, y_min, y_max)
        return create_correlation_plot(filtered_df, x_col, y_col)
    else:
        return go.Figure()  # Return an empty figure if any input is None


@callback(
    Output('correlation-plot-2', 'figure'),
    [Input('dataframe-json', 'children'),
     Input('xaxis-column-2', 'value'),
     Input('yaxis-column-2', 'value'),
     Input('x-min-input-2', 'value'),
     Input('x-max-input-2', 'value'),
     Input('y-min-input-2', 'value'),
     Input('y-max-input-2', 'value')]
)
def update_correlation_plot_2(csv_json, x_col, y_col, x_min, x_max, y_min, y_max):
    if csv_json is not None and x_col is not None and y_col is not None:
        df = pd.read_json(io.StringIO(csv_json), orient='split')
        filtered_df = filter_dataframe(df, x_col, y_col, x_min, x_max, y_min, y_max)
        return create_correlation_plot(filtered_df, x_col, y_col)
    else:
        return go.Figure()  # Return an empty figure if any input is None

@callback(
    Output('3d-scatter-plot', 'figure'),
    [Input('dataframe-json', 'children'),
     Input('xaxis-column-3d', 'value'),
     Input('yaxis-column-3d', 'value'),
     Input('zaxis-column-3d', 'value'),
     Input('x-min-input-3d', 'value'),
     Input('x-max-input-3d', 'value'),
     Input('y-min-input-3d', 'value'),
     Input('y-max-input-3d', 'value'),
     Input('z-min-input-3d', 'value'),
     Input('z-max-input-3d', 'value'),]
)
def update_scatter3d_plot(csv_json, x_col, y_col, z_col, x_min, x_max, y_min, y_max, z_min, z_max):
    if csv_json is not None and x_col is not None and y_col is not None and z_col is not None:
        df = pd.read_json(io.StringIO(csv_json), orient='split')
        filtered_df = filter_dataframe_3d_scatter(df, x_col, y_col, z_col, x_min, x_max, y_min, y_max, z_min, z_max)
        return create_3d_correlation_plot(filtered_df, x_col, y_col, z_col)
    else:
        return go.Figure()

@callback(
    [Output('common-data-table', 'data'),
     Output('filtered-data', 'children')],
    [Input('dataframe-json', 'children'),
     Input('correlation-plot-1', 'selectedData'),
     Input('correlation-plot-2', 'selectedData'),
     Input('x-min-input-1', 'value'),
     Input('x-max-input-1', 'value'),
     Input('y-min-input-1', 'value'),
     Input('y-max-input-1', 'value'),
     Input('x-min-input-2', 'value'),
     Input('x-max-input-2', 'value'),
     Input('y-min-input-2', 'value'),
     Input('y-max-input-2', 'value'),
     Input('z-min-input-3d', 'value'),
     Input('z-max-input-3d', 'value'),
     Input('xaxis-column-1', 'value'),
     Input('yaxis-column-1', 'value'),
     Input('xaxis-column-2', 'value'),
     Input('yaxis-column-2', 'value'),
     Input('zaxis-column-3d', 'value')]
)
def update_common_data_table(csv_json, selected_correlation_1, selected_correlation_2,
                             x_min_1, x_max_1, y_min_1, y_max_1,
                             x_min_2, x_max_2, y_min_2, y_max_2,
                             z_min_1, z_max_1, x_col_1, y_col_1,
                             x_col_2, y_col_2, z_col_1):
    try:
        if csv_json:
            df = pd.read_json(io.StringIO(csv_json), orient='split')
            # Apply filters based on input box values and dropdown selections
            filtered_df_1 = filter_dataframe(df, x_col_1, y_col_1, x_min_1, x_max_1, y_min_1, y_max_1)
            filtered_df_2 = filter_dataframe(df, x_col_2, y_col_2, x_min_2, x_max_2, y_min_2, y_max_2)
            filtered_df_3 = filter_dataframe_3d_scatter(df, x_col_1, y_col_1, z_col_1, x_min_1, x_max_1, y_min_1,
                                                        y_max_1, z_min_1, z_max_1)

            # Merging with suffixes

            common_data = pd.merge(pd.merge(filtered_df_1, filtered_df_2, on='description', how='inner',
                                            suffixes=('', '_drop')), filtered_df_3, on='description',
                                            how='inner', suffixes=('', '_drop'),
                                   )
            common_data = common_data.loc[:, ~common_data.columns.str.contains('_drop')]

            # Check if common_data is valid and not None
            if common_data is not None and not common_data.empty:
                data_for_table = common_data.to_dict('records')
                json_data = common_data.to_json(date_format='iso', orient='split')
                return data_for_table, json_data
            else:
                return [], None  # Return empty list and None for no data
        else:
            return [], None
    except Exception as e:
        # Handle any exceptions
        print(f"Error in update_common_data_table: {e}")
        return [], None  # Return empty list and None in case of an error


# Callback to download CSV
@callback(
    Output("download-dataframe-csv", "data"),
    [Input("btn_csv", "n_clicks"),
     Input('common-data-table', 'data')],  # Use filtered data
    prevent_initial_call=True
)
def download_csv(n_clicks, filtered_json):
    if n_clicks is not None and filtered_json is not None:
        # Check if filtered_json is a string (JSON)
        if isinstance(filtered_json, str):
            filtered_df = pd.read_json(io.StringIO(filtered_json), orient='split')
        # If filtered_json is a list or dict, convert directly to DataFrame
        else:
            filtered_df = pd.DataFrame(filtered_json)

        return dcc.send_data_frame(filtered_df.to_csv, filename="common_data_points.csv")
    else:
        raise dash.exceptions.PreventUpdate

# Silent File dropdown
@callback(
    Output('my-dropdownSil', 'options'),
    [Input('storedSil', 'data')],
    prevent_initial_call=True
)
def update_dropdown_options(data):
    options = [{'label': value, 'value': value} for value in data]
    return options

# Callback to display the selected item from the dropdown
@callback(
    Output('selectedSil', 'children'),
    [Input('my-dropdownSil', 'value')],
    prevent_initial_call=True
)
def display_selected_item(selected_value):
    return selected_value

# output directory dropdown
@callback(
    Output('my-dropdownDir', 'options'),
    [Input('storedOut', 'data')],
    prevent_initial_call=True
)
def update_dropdown_options(data):
    options = [{'label': value, 'value': value} for value in data]
    return options

# Callback to display the selected item from the dropdown
@callback(
    Output('selectedDir', 'children'),
    [Input('my-dropdownDir', 'value')],
    prevent_initial_call=True
)
def display_selected_item(selected_value):
    return selected_value

# Callback to get the PDBS from the silent files based on the common_data_points
@callback(
    Output('get-pdb-list', 'children'),
    [Input('get-pdb-button', 'n_clicks')],
    [State('common-data-table', 'data')],
)
def get_pdb_list(n_clicks, column_data):
    if n_clicks > 0:
        column_contents = [row['design_name'].rstrip('.pdb') for row in column_data]
        return json.dumps(column_contents)
    return None

# Callback to eventually retrieve the associated pdb files with their design_name
@callback([
    Output('pdb-file-list', 'children'),
    Output('pdb-path-list', 'children')],
    [Input('get-pdb-list', 'children'),
     Input('selectedSil','children'),
     Input('selectedDir','children')
     ],
    background=True, manager=background_callback_manager
)
def get_pdb_files(pdb_list,selectedSil, selectedDir):
    if pdb_list is not None:
        pdb_file_list = json.loads(pdb_list)
        pdb_path_list = get_top_pdbs_from_silent(pdb_file_list, selectedSil, selectedDir)
        return pdb_file_list, pdb_path_list
    return 'Click the button and wait for data to be processed.', None

