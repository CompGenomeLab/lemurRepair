import plotly.graph_objects as go
import pandas as pd


def plot_frequency(infile, outfile):
    with open(infile) as in_handler:
        df = pd.read_json(in_handler)
        fig = go.Figure(
            data=[
                go.Bar(name='T', x=df.index, y=df['T'], marker={'color':'#1f76b4'}),
                go.Bar(name='C', x=df.index, y=df['C'], marker={'color':'#ff7e0e'}),
                go.Bar(name='G', x=df.index, y=df['G'], marker={'color':'#2da02c'}),
                go.Bar(name='A', x=df.index, y=df['A'], marker={'color':'#d62728'})
                ],
            layout=go.Layout(
                title='Nucleotide distribution of 26nt oligomer',
                barmode='stack',
                yaxis={'title': 'Frequency'},
                xaxis={'title': 'Position'}
            )
        )
        fig.write_image(outfile)

if __name__ == "__main__":
    plot_frequency("lemur_frequency.json", "lemur_frequency.pdf")
    plot_frequency("human_frequency.json", "human_frequency.pdf")