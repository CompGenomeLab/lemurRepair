import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd


def percentage(row, reads_sum):
    return row['count']/reads_sum*100


def plot_length_dist(infile_human_cpd, infile_lemur_cpd, infile_human_64, infile_lemur_64, outfile):
    fig = make_subplots(rows=2, cols=2)

    with open(infile_human_cpd) as in_handler:
        df = pd.read_csv(in_handler, sep=',')
        reads_sum = df['count'].sum()
        df['read_percentage'] = df.apply(lambda row: percentage(row, reads_sum), axis=1)
        fig.add_trace(
            go.Bar(name='length', x=df['index'], y=df['read_percentage'], marker={'color':'#1f76b4'}),
            row=1, col=1
            )

    with open(infile_lemur_cpd) as in_handler:
        df = pd.read_csv(in_handler, sep=',')
        reads_sum = df['count'].sum()
        df['read_percentage'] = df.apply(lambda row: percentage(row, reads_sum), axis=1)
        fig.add_trace(
            go.Bar(name='length', x=df['index'], y=df['read_percentage'], marker={'color':'#1f76b4'}),
            row=1, col=2
            )

    with open(infile_human_64) as in_handler:
        df = pd.read_csv(in_handler, sep=',')
        reads_sum = df['count'].sum()
        df['read_percentage'] = df.apply(lambda row: percentage(row, reads_sum), axis=1)
        fig.add_trace(
            go.Bar(name='length', x=df['index'], y=df['read_percentage'], marker={'color':'#1f76b4'}),
            row=2, col=1
            )

    with open(infile_lemur_64) as in_handler:
        df = pd.read_csv(in_handler, sep=',')
        reads_sum = df['count'].sum()
        df['read_percentage'] = df.apply(lambda row: percentage(row, reads_sum), axis=1)
        fig.add_trace(
            go.Bar(name='length', x=df['index'], y=df['read_percentage'], marker={'color':'#1f76b4'}),
            row=2, col=2
            )

    fig.update_xaxes(title_text="Read length (nt)", row=1, col=1)
    fig.update_xaxes(title_text="Read length (nt)", row=1, col=2)
    fig.update_xaxes(title_text="Read length (nt)", row=2, col=1)
    fig.update_xaxes(title_text="Read length (nt)", row=2, col=2)

    fig.update_yaxes(title_text="Read count (%)", range=[0, 25], row=1, col=1)
    fig.update_yaxes(title_text="Read count (%)", range=[0, 25], row=1, col=2)
    fig.update_yaxes(title_text="Read count (%)", range=[0, 25], row=2, col=1)
    fig.update_yaxes(title_text="Read count (%)", range=[0, 25], row=2, col=2)

    fig.update_layout(showlegend=False, title_text="Read length distribution", height=1200, width=1600)
    fig.write_image(outfile)


if __name__ == "__main__":
    plot_length_dist("human_cpd.csv", "human_64.csv", "lemur_cpd.csv", "lemur_64.csv", "fig-1.A.pdf")
