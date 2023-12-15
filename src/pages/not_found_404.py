import dash
from dash import html, dcc

dash.register_page(
    __name__,
    path='/404')

layout = html.Div(
    [
        html.H1('404'),
        html.H2('Page not found'),
        html.H2('Oh, something went wrong!'),
        dcc.Link('Go back to home', href='/'),
    ]
)
