import dash
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template
from . import data_vis_callbacks

dash.register_page(
    __name__,
    path='/plots',
    title='plots',
    )

load_figure_template("CYBORG")

layout = html.Div([
    html.H1("Data Visualization Dashboard"),

    html.H2("Upload Files"),
    dcc.Upload(
        id='upload-csv',
        children=html.Div(['Drag and Drop or ', html.A('Select csv File')]),
        style={
            'width': '25%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px',
        },
        # Allow only a single file to be uploaded
        multiple=False,
        accept='.csv'
    ),

    # Scatter Plot and Radar Plot
    html.Div([
        html.Div([
            html.H2('Scatter Plot'),
            dcc.Graph(id='scatter-plot'),
            dcc.Dropdown(
                id='dropdown',
                style={'width': '60%',
                       },
                options=[],  # Set initial options to empty
                value=None
            ),
        ], style={'width': '48%',
                  'display': 'inline-block',
                  }
        ),

        html.Div([
            html.H2('Radar Plot'),
            dcc.Graph(id='radar-plot')
        ], style={'width': '48%',
                  'display': 'inline-block',
                  },
        ),
    ]),

    # Correlation Plot 1
    html.Div([
        html.H2("Correlation Plot 1"),
        dcc.Dropdown(
            id='xaxis-column-1',
            options=[],  # Set initial options to empty
            value=None,
            style={'width': '35%',
                   },
        ),
        dcc.Dropdown(
            id='yaxis-column-1',
            options=[],  # Set initial options to empty
            value=None,
            style={'width': '35%',
                   },
        ),
        dcc.Graph(id='correlation-plot-1'),
        html.Div([
            html.Label("X-axis range: "),
            dcc.Input(
                id='x-min-input-1',
                type='number',
                placeholder='Min X',
                style={'width': '10%',
                       },
                ),
            dcc.Input(id='x-max-input-1',
                      type='number',
                      placeholder='Max X',
                      style={'width': '10%',
                             },
                      ),
            html.Label("Y-axis range: "),
            dcc.Input(id='y-min-input-1',
                      type='number',
                      placeholder='Min Y',
                      style={'width': '10%',
                             },
                      ),
            dcc.Input(id='y-max-input-1',
                      type='number',
                      placeholder='Max Y',
                      style={'width': '10%',
                             },
                      ),
        ]),
    ]),
    # Correlation Plot 2
    html.Div([
        html.H2("Correlation Plot 2"),
        dcc.Dropdown(
            id='xaxis-column-2',
            options=[],  # Set initial options to empty
            value=None,
            style={'width': '35%',
                   },
        ),
        dcc.Dropdown(
            id='yaxis-column-2',
            options=[],  # Set initial options to empty
            value=None,
            style={'width': '35%',
                   },
        ),
        dcc.Graph(id='correlation-plot-2'),
        html.Div([
            html.Label("X-axis range: "),
            dcc.Input(id='x-min-input-2',
                        type='number',
                        placeholder='Min X',
                        style={'width': '10%',
                               },
                        ),
            dcc.Input(id='x-max-input-2',
                        type='number',
                        placeholder='Max X',
                        style={'width': '10%',
                               },
                        ),
            html.Label("Y-axis range: "),
            dcc.Input(id='y-min-input-2',
                        type='number',
                        placeholder='Min Y',
                        style={'width': '10%',
                               },
                        ),
            dcc.Input(id='y-max-input-2',
                        type='number',
                        placeholder='Max Y',
                        style={'width': '10%',
                               },
                        ),
        ]),
    ]),

    html.Div([
        html.H2("3d Scatter Plot"),
        dcc.Dropdown(
            id='xaxis-column-3d',
            options=[],  # Set initial options to empty
            value=None,
            style={'width': '35%',
                   },
        ),
        dcc.Dropdown(
            id='yaxis-column-3d',
            options=[],  # Set initial options to empty
            value=None,
            style={'width': '35%',
                   },
        ),
        dcc.Dropdown(
            id='zaxis-column-3d',
            options=[],  # Set initial options to empty
            value=None,
            style={'width': '35%',
                   },
        ),
        dcc.Graph(id='3d-scatter-plot',
                  style={
                      'width': '90vw',
                      'height': '90vh',
                  }),
        html.Div([
                    html.Label("X-axis range: "),
                    dcc.Input(id='x-min-input-3d',
                                type='number',
                                placeholder='Min X',
                                style={'width': '10%',
                                       },
                                ),
                    dcc.Input(id='x-max-input-3d',
                                type='number',
                                placeholder='Max X',
                                style={'width': '10%',
                                       },
                                ),
                    html.Label("Y-axis range: "),
                    dcc.Input(id='y-min-input-3d',
                                type='number',
                                placeholder='Min Y',
                                style={'width': '10%',
                                       },
                                ),
                    dcc.Input(id='y-max-input-3d',
                                type='number',
                                placeholder='Max Y',
                                style={'width': '10%',
                                       },
                                ),
                    html.Label("z-axis range: "),
                    dcc.Input(id='z-min-input-3d',
                                type='number',
                                placeholder='Min Z',
                                style={'width': '10%',
                                       },
                                ),
                    dcc.Input(id='z-max-input-3d',
                                type='number',
                                placeholder='Max Z',
                                style={'width': '10%',
                                       },
                                ),
    ]),

    # Common Data Points Table
    html.H2("Common Data Points"),
    dash_table.DataTable(
        id='common-data-table',
        style_data=
        {
            'color': 'black',
            'backgroundColor': 'white'
        },
        style_data_conditional=[
            {
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(220, 220, 220)',
            }
        ],
        style_header={
            'backgroundColor': 'rgb(210, 210, 210)',
            'color': 'black',
            'fontWeight': 'bold'
        }
    ),

    html.Div([
        dbc.Button('Download CSV', id='btn_csv', color='secondary'),
        dcc.Download(id='download-dataframe-csv')],
        style={'padding-bottom': '20px'}
    ),

    html.Div([
       dbc.Button('Get pdb files', id='get-pdb-button', color='secondary', n_clicks=0)
    ]),

    html.Div([
        html.P(id='pdb-file-list')
    ]),

    html.Div([
        html.P(id='pdb-path-list')
    ]),

    # Hidden Divs
    html.Div(id='dataframe-json', style={'display': 'none'}),
    html.Div(id='filtered-data', style={'display': 'none'}),

    html.Div(id='scatter-plot-click', style={'display': 'none'}),

    html.Div(id='get-pdb-list', style={'display': 'none'}),
    ]),
])
