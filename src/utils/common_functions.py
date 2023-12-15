import plotly.graph_objs as go


def generate_hover_text(row, columns):
    return '<br>'.join([f'{col}: {row[col]}' for col in columns])


def create_correlation_plot(df, x_column, y_column):
    sorted_df = df.sort_values(by=x_column)

    # Check if y_column is numeric and handle NaNs
    if sorted_df[y_column].dtype in ['float64', 'int64']:
        sorted_df = sorted_df.dropna(subset=[y_column])
        color = sorted_df[y_column]
    else:
        color = 'blue'  # Fallback to fixed color if non-numeric

    fig = go.Figure(data=[
        go.Scatter(
            x=sorted_df[x_column],
            y=sorted_df[y_column],
            mode='markers',
            marker=dict(color=color, colorscale='Viridis', showscale=True),
            text=sorted_df.apply(lambda row: generate_hover_text(row, df.columns), axis=1),
            hoverinfo='text'
        )
    ])
    fig.update_layout(
        title=f"Correlation between {x_column} and {y_column}",
        xaxis_title=x_column,
        yaxis_title=y_column
    )
    return fig


def create_3d_correlation_plot(df, x_column, y_column, z_column):
    sorted_df = df.sort_values(by=x_column)

    # Check if y_column is numeric and handle NaNs
    if sorted_df[y_column].dtype in ['float64', 'int64']:
        sorted_df = sorted_df.dropna(subset=[y_column])
        color = sorted_df[y_column]
    else:
        color = 'blue'

    fig = go.Figure(data=[
        go.Scatter3d(
            x=sorted_df[x_column],
            y=sorted_df[y_column],
            z=sorted_df[z_column],
            mode='markers',
            marker=dict(color=color,
                        colorscale='Viridis',
                        showscale=True,
                        size=1),
            text=sorted_df.apply(lambda row: generate_hover_text(row, df.columns), axis=1),
            hoverinfo='text',
        )
    ])
    fig.update_layout(
        title=f'Correlation between {x_column}, {y_column}, and {z_column}',
        scene=dict(
            xaxis=dict(title=x_column),
            yaxis=dict(title=y_column),
            zaxis=dict(title=z_column),
        )
    )
    return fig

def filter_dataframe(df, x_col, y_col, x_min, x_max, y_min, y_max):
    filtered_df = df.copy()
    if x_min is not None:
        filtered_df = filtered_df[filtered_df[x_col] >= x_min]
    if x_max is not None:
        filtered_df = filtered_df[filtered_df[x_col] <= x_max]
    if y_min is not None:
        filtered_df = filtered_df[filtered_df[y_col] >= y_min]
    if y_max is not None:
        filtered_df = filtered_df[filtered_df[y_col] <= y_max]
    return filtered_df

def filter_dataframe_3d_scatter(df, x_col, y_col, z_col, x_min, x_max, y_min, y_max, z_min, z_max):
    filtered_df = df.copy()
    if x_min is not None:
        filtered_df = filtered_df[filtered_df[x_col] >= x_min]
    if x_max is not None:
        filtered_df = filtered_df[filtered_df[x_col] <= x_max]
    if y_min is not None:
        filtered_df = filtered_df[filtered_df[y_col] >= y_min]
    if y_max is not None:
        filtered_df = filtered_df[filtered_df[y_col] <= y_max]
    if z_min is not None:
        filtered_df = filtered_df[filtered_df[z_col] >= z_min]
    if z_max is not None:
        filtered_df = filtered_df[filtered_df[z_col] <= z_max]
    return filtered_df
