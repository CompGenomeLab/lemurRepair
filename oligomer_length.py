import plotly.graph_objects as go
import pandas as pd


def plot_length(infile, outfile):
    with open(infile) as in_handler:
        df = pd.read_csv(in_handler, sep=',')
        fig = go.Figure(
            data=[
                go.Bar(name='length', x=df['index'], y=df['count'], marker={'color':'#1f76b4'})
                ],
            layout=go.Layout(
                title='Length distribution',
                yaxis={'title': 'Count'},
                xaxis={'title': 'Length'}
            )
        )
        fig.write_image(outfile)


if __name__ == "__main__":
    plot_length("lemur_nuc_len_dist.csv", "lemur_nuc_len_dist.pdf")
    plot_length("human_nuc_len_dist.csv", "human_nuc_len_dist.pdf")