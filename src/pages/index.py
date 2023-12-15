import dash
from dash import html


dash.register_page(
    __name__,
    path='/',
    redirect_from=['/index'],
    title='RFdiffusion Analytics',
)


layout = html.Div(
    [
        html.H1('Protein Analytics'),
        html.Div(
            html.A('Go to Dashboard', href='/plots')
        ),
        html.Div(id='content'),
    ]
)
