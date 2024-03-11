import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as pg


def generate_truthtables(indexes_to_exclude=[]):
    correct_matrix = np.eye(39)
    correct_entries = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8),
                   (2, 9), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (5, 6),
                   (5, 7), (5, 8), (5, 9), (6, 7),
                   (6, 8), (6, 9), (7, 8), (7, 9),
                   (8, 9), (10, 11), (12, 13), (14, 15), (16, 17), (18, 19), (20, 21), (22, 23), (22, 24), (23, 24), (25, 26), (27,28), (29, 30), (31, 32), (33, 34), (35, 36), (37, 38)]
    for point in correct_entries:
        correct_matrix[point] = 1 

    correct_matrix += correct_matrix.T
    correct_matrix = correct_matrix.astype(bool)

    mask = np.ones(correct_matrix.shape)
    for i in range(39):
        if i in indexes_to_exclude:
            mask[i,:] = 0
            mask[:,i] = 0
    mask = mask.astype(bool)

    positive_matrix = correct_matrix.copy()*mask
    negative_matrix = np.invert(correct_matrix)*mask
    return positive_matrix, negative_matrix


def plot_matrix(mat, midpoint=0.5):
    fig = pg.Figure()
    fig.add_trace(
        pg.Heatmap(
            z=mat.astype(float),
            colorscale=[(0, "#d7191c"), (0.4, "#fdae61"), (0.5, "#ffffbf"), (0.6, "#a6d96a"),(1, "#1a9641")],
            zmid=midpoint,
        )
    )

    fig.update_yaxes(autorange="reversed")
    fig.update_xaxes(side="top")

    fig.update_xaxes(
        tickfont_size=28/3, 
    )
    fig.update_yaxes(
        tickfont_size=28/3,  
    )
    fig.update_layout(
        width=300,
        height=300,
        template="simple_white",
        margin=dict(t=30, b=0, l=0, r=10),
        showlegend=False,
        title_font_size=10,
        title_font_family="Inter", 
    )

    fig.update_layout(
        template="simple_white",
        width=300,
        height=300,
        showlegend=False,
        margin=dict(t=10, b=10, l=10, r=10),
        font_family="Inter",
        legend_font_size=28/3,
    )
    fig.update_traces(
        showscale=False, 
    )

    return fig



def plot_hists(mat, pos_matrix, neg_matrix, title="Unknown"):

    posdata = np.where(pos_matrix, mat, np.nan)
    negdata = np.where(neg_matrix, mat, np.nan)

    threshold_width = np.nanmin(posdata) - np.nanmax(negdata)
    threshold_middle = (np.nanmin(posdata) + np.nanmax(negdata))/2
    print(threshold_middle)

    df = pd.DataFrame.from_dict({'score': mat.flatten(), 'pos': pos_matrix.flatten(), 'neg': neg_matrix.flatten(), 'diag': np.eye(39, dtype=bool).flatten()})
    df['cat'] = "None"
    df.loc[df.pos == True, 'cat'] = "Pos"
    df.loc[df.neg == True, 'cat'] = "Neg"
    df.loc[df.diag == True, 'cat'] = "None"

    df = df.drop(df.loc[df.cat == "None"].index)

    fig = px.histogram(
        df, 
        x="score", 
        color="cat",
        histnorm="percent", 
        barmode='overlay',
        opacity=0.7,
        color_discrete_map={"Pos": "#f3ec17", "Neg": "#0e71b9", "None": "#bdbdbd"},
        title=f"{title}, d={threshold_width:.3f}",
        nbins=20,
    )
    fig.update_traces(xbins=dict(start=0.0, end=1.05, size=0.05))
    fig.update_xaxes(
        range=[0, 1.05],
        title_text='Similarity score', 
        title_font_family="Inter", 
        title_font_size=28/3, 
        tickfont_size=28/3, 
    )
    fig.update_yaxes(
        title_text='Percent', 
        title_font_family="Inter", 
        title_font_size=28/3, 
        tickfont_size=28/3, 
        minor_ticks="outside", 
    )
    fig.update_layout(
        width=300,
        height=300,
        template="simple_white",
        margin=dict(t=30, b=45, l=0, r=10),
        showlegend=False,
        title_font_size=10,
        title_font_family="Inter", 
        font_family="Inter",
        legend_font_size=28/3,
    )

    return fig



def reference_matrix(pos_matrix, neg_matrix):
    reference_data = np.empty_like(pos_matrix, dtype=float)
    reference_data[:] = np.NaN
    reference_data[pos_matrix] = 1
    reference_data[neg_matrix] = 0
    return plot_matrix(reference_data)