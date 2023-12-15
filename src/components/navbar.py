import dash_bootstrap_components as dbc

navbar = dbc.NavbarSimple(
    children=[
        # dbc.NavItem(dbc.NavLink("Page 1", href="#")),
        dbc.DropdownMenu(
            children=[
                # dbc.DropdownMenuItem("Data Processing", header=True),
                dbc.DropdownMenuItem("Data Processing", href="/plots"),
                dbc.DropdownMenuItem("PDB Inspection", href="/molview"),
            ],
            nav=True,
            in_navbar=True,
            label="Tools",
        ),
    ],
    brand="RFdiffusion Data Visualization",
    brand_href="/",
    color="primary",
    dark=True,
)