import plotly.offline as opy
import plotly.express as px
from plotly.offline import plot
import plotly.graph_objs as go


import pandas as pd

STATIC_URL = "static/"


def make_stacked_bar(csv_name):
    df = pd.read_csv(STATIC_URL + csv_name)
    fig = go.Figure()
    for label in df.columns.values[1:]:
        fig.add_trace(go.Bar(x=df['date'], y=df[label], name=label))

    # Change the bar mode
    fig.update_layout(
        barmode='stack',
        xaxis=go.layout.XAxis(
            title=go.layout.xaxis.Title(
                text="Date",
                font=dict(
                    size=14,
                    color="#7f7f7f"))),
        yaxis=go.layout.YAxis(
            title=go.layout.yaxis.Title(
                text="Number of Reactions",
                font=dict(
                    size=14,
                    color="#7f7f7f")),
            type="log"))
    plot_div = plot(fig, output_type='div', include_plotlyjs=False)

    return plot_div


def make_valid_reaction_bar_char(csv_name):
    df = pd.read_csv(STATIC_URL + csv_name)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df['date'], y=df['Valid'],
                             mode='lines',
                             name='Valid',
                             fill='tozeroy'))
    fig.add_trace(go.Scatter(x=df['date'], y=df['Invalid'],
                             mode='lines',
                             name='Invalid',
                             fill='tozeroy'))

    fig.update_layout(
        xaxis=go.layout.XAxis(
            title=go.layout.xaxis.Title(
                text="Date",
                font=dict(
                    size=14,
                    color="#7f7f7f"))),
        yaxis=go.layout.YAxis(
            title=go.layout.yaxis.Title(
                text="Number of Reactions",
                font=dict(
                    size=14,
                    color="#7f7f7f"))))
    plot_div = plot(fig, output_type='div', include_plotlyjs=False)
    return plot_div


def make_area_chart(csv_name):
    print(csv_name)
    df = pd.read_csv(STATIC_URL + csv_name)
    fig = go.Figure()
    for i, label in enumerate(df.columns.values[1:]):
        if i == 0:
            fig.add_trace(go.Scatter(
                x=df['date'], y=df[label], name=label, fill='tozeroy'))
        else:
            fig.add_trace(go.Scatter(
                x=df['date'], y=df[label], name=label, fill='tozeroy'))

    plot_div = plot(fig, output_type='div', include_plotlyjs=False)
    return plot_div
