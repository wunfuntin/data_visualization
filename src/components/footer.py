# notes
"""
This file is for creating a simple footer element.
This component will sit at the bottom of each page of the application.
"""

# package imports
from dash import html
import dash_bootstrap_components as dbc

footer = html.Footer(
    dbc.Container(
        [
            html.Hr(),
            'Copyright© 2023 Ryan G. Walker',
            html.Br(),
            'All rights reserved.'
        ]
    )
)