# Package imports
import dash
from dash import html, dcc
import dash_bootstrap_components as dbc
from flask import Flask
import multiprocessing as mp

# Local imports
from components import navbar, footer


server = Flask(__name__)

app = dash.Dash(
    __name__,
    server=server,
    use_pages=True,
    external_stylesheets=[
        dbc.themes.CYBORG,
        dbc.icons.FONT_AWESOME,
    ],
    meta_tags=[
        {
            'name': 'viewport',
            'content': 'width=device-width, initial-scale=1'
        }
    ],
    suppress_callback_exceptions=True,
    title='RFdiffusion Data Visualization',
)

server = app.server

app.layout = html.Div(
    [
        dcc.Location(id='url', refresh=False),
        navbar,
        dbc.Container(dash.page_container, class_name='my-2'),
        footer,
    ]
)


if __name__ == '__main__':
    mp.set_start_method('spawn')
    app.run_server(
        debug=True
    )
