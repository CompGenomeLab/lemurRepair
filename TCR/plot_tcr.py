import plotly.graph_objects as go
from plotly.subplots import make_subplots


def calc_tcr(infile, tsnts):
    with open(infile, 'r') as bed_handler:
        tcr_array_ts = [0 for _ in range(201)]
        tcr_array_nts = [0 for _ in range(201)]
        for line in bed_handler:
            line = line.strip().split()
            if f'{line[6]}' == tsnts:
                xrmid = int(line[1]) + (int(line[2])-int(line[1])) // 2
                wbin = (xrmid - int(line[8]))//100 + 1
                if wbin > 201:
                    wbin = 201
                tcr = 'ts' if f'{line[-1]}' != f'{line[5]}' else 'nts'
                if tcr == 'ts':
                    tcr_array_ts[wbin-1] += 1
                else:
                    tcr_array_nts[wbin-1] += 1
                wbin = wbin-100
                if f'{line[-1]}' == '-':
                    wbin = 0 - wbin
                liststring = '\t'.join(line)
        tcr_array_nts = tcr_array_nts[::-1]
        return tcr_array_ts, tcr_array_nts


def plot_tcr(outfile):

    human_tss_cpd_xr_ts, human_tss_cpd_xr_nts = calc_tcr("human_tss_cpd_tcr.bed", "tss")
    human_tss_64_xr_ts, human_tss_64_xr_nts = calc_tcr("human_tss_64_tcr.bed", "tss")

    human_tes_cpd_xr_ts, human_tes_cpd_xr_nts = calc_tcr("human_tes_cpd_tcr.bed", "tes")
    human_tes_64_xr_ts, human_tes_64_xr_nts = calc_tcr("human_tes_64_tcr.bed", "tes")


    lemur_tss_cpd_xr_ts, lemur_tss_cpd_xr_nts = calc_tcr("lemur_tss_cpd_tcr.bed", "tss")
    lemur_tss_64_xr_ts, lemur_tss_64_xr_nts = calc_tcr("lemur_tss_64_tcr.bed", "tss")

    lemur_tes_cpd_xr_ts, lemur_tes_cpd_xr_nts = calc_tcr("lemur_tes_cpd_tcr.bed", "tes")
    lemur_tes_64_xr_ts, lemur_tes_64_xr_nts = calc_tcr("lemur_tes_64_tcr.bed", "tes")

    ##---------------------------------------------------------------

    shuffled_human_tss_cpd_xr_ts, shuffled_human_tss_cpd_xr_nts = calc_tcr("human_tss_cpd_shuffled_tcr.bed", "tss")
    shuffled_human_tss_64_xr_ts, shuffled_human_tss_64_xr_nts = calc_tcr("human_tss_64_shuffled_tcr.bed", "tss")

    shuffled_human_tes_cpd_xr_ts, shuffled_human_tes_cpd_xr_nts = calc_tcr("human_tes_cpd_shuffled_tcr.bed", "tes")
    shuffled_human_tes_64_xr_ts, shuffled_human_tes_64_xr_nts = calc_tcr("human_tes_64_shuffled_tcr.bed", "tes")


    shuffled_lemur_tss_cpd_xr_ts, shuffled_lemur_tss_cpd_xr_nts = calc_tcr("lemur_tss_cpd_shuffled_tcr.bed", "tss")
    shuffled_lemur_tss_64_xr_ts, shuffled_lemur_tss_64_xr_nts = calc_tcr("lemur_tss_64_shuffled_tcr.bed", "tss")

    shuffled_lemur_tes_cpd_xr_ts, shuffled_lemur_tes_cpd_xr_nts = calc_tcr("lemur_tes_cpd_shuffled_tcr.bed", "tes")
    shuffled_lemur_tes_64_xr_ts, shuffled_lemur_tes_64_xr_nts = calc_tcr("lemur_tes_64_shuffled_tcr.bed", "tes")

    ##---------------------------------------------------------------

    human_tss_cpd_xr_ts = [i / j for i, j in zip(human_tss_cpd_xr_ts, shuffled_human_tss_cpd_xr_ts)]
    human_tss_cpd_xr_nts = [i / j for i, j in zip(human_tss_cpd_xr_nts, shuffled_human_tss_cpd_xr_nts)]
    human_tss_64_xr_ts = [i / j for i, j in zip(human_tss_64_xr_ts, shuffled_human_tss_64_xr_ts)]
    human_tss_64_xr_nts = [i / j for i, j in zip(human_tss_64_xr_nts, shuffled_human_tss_64_xr_nts)]
    human_tes_cpd_xr_ts = [i / j for i, j in zip(human_tes_cpd_xr_ts, shuffled_human_tes_cpd_xr_ts)]
    human_tes_cpd_xr_nts = [i / j for i, j in zip(human_tes_cpd_xr_nts, shuffled_human_tes_cpd_xr_nts)]

    human_tes_64_xr_ts = [i / j for i, j in zip(human_tes_64_xr_ts, shuffled_human_tes_64_xr_ts)]
    human_tes_64_xr_nts = [i / j for i, j in zip(human_tes_64_xr_nts, shuffled_human_tes_64_xr_nts)]

    lemur_tss_cpd_xr_ts = [i / j for i, j in zip(lemur_tss_cpd_xr_ts, shuffled_lemur_tss_cpd_xr_ts)]
    lemur_tss_cpd_xr_nts = [i / j for i, j in zip(lemur_tss_cpd_xr_nts, shuffled_lemur_tss_cpd_xr_nts)]
    lemur_tss_64_xr_ts = [i / j for i, j in zip(lemur_tss_64_xr_ts, shuffled_lemur_tss_64_xr_ts)]
    lemur_tss_64_xr_nts = [i / j for i, j in zip(lemur_tss_64_xr_nts, shuffled_lemur_tss_64_xr_nts)]
    lemur_tes_cpd_xr_ts = [i / j for i, j in zip(lemur_tes_cpd_xr_ts, shuffled_lemur_tes_cpd_xr_ts)]
    lemur_tes_cpd_xr_nts = [i / j for i, j in zip(lemur_tes_cpd_xr_nts, shuffled_lemur_tes_cpd_xr_nts)]

    lemur_tes_64_xr_ts = [i / j for i, j in zip(lemur_tes_64_xr_ts, shuffled_lemur_tes_64_xr_ts)]
    lemur_tes_64_xr_nts = [i / j for i, j in zip(lemur_tes_64_xr_nts, shuffled_lemur_tes_64_xr_nts)]
    ##---------------------------------------------------------------

    x_axis = [i/10 for i in range(-100,101)]

    fig = make_subplots(rows=2, cols=4, shared_yaxes=True, horizontal_spacing=0.03, vertical_spacing=0.2)

    fig.add_traces(
        [go.Scatter(name='TS', x=x_axis, y=human_tss_cpd_xr_ts, line={'color':'#d62728'}),
        go.Scatter(name='NTS', x=x_axis, y=human_tss_cpd_xr_nts, line={'color':'#1f76b4'})],
        rows=[1,1], cols=[1,1]
    )

    fig.add_traces(
        [go.Scatter(name='TS', x=x_axis, y=human_tes_cpd_xr_ts, line={'color':'#d62728'}, showlegend=False),
        go.Scatter(name='NTS', x=x_axis, y=human_tes_cpd_xr_nts, line={'color':'#1f76b4'}, showlegend=False)],
        rows=[1,1], cols=[2,2]
    )

    fig.add_traces(
        [go.Scatter(name='TS', x=x_axis, y=lemur_tss_cpd_xr_ts, line={'color':'#d62728'}, showlegend=False),
        go.Scatter(name='NTS', x=x_axis, y=lemur_tss_cpd_xr_nts, line={'color':'#1f76b4'}, showlegend=False)],
        rows=[1,1], cols=[3,3]
    )

    fig.add_traces(
        [go.Scatter(name='TS', x=x_axis, y=lemur_tes_cpd_xr_ts, line={'color':'#d62728'}, showlegend=False),
        go.Scatter(name='NTS', x=x_axis, y=lemur_tes_cpd_xr_nts, line={'color':'#1f76b4'}, showlegend=False)],
        rows=[1,1], cols=[4,4]
    )

    #---------------------------------------------------------------------------

    fig.add_traces(
        [go.Scatter(name='TS', x=x_axis, y=human_tss_64_xr_ts, line={'color':'#d62728'}, showlegend=False),
        go.Scatter(name='NTS', x=x_axis, y=human_tss_64_xr_nts, line={'color':'#1f76b4'}, showlegend=False)],
        rows=[2,2], cols=[1,1]
    )

    fig.add_traces(
        [go.Scatter(name='TS', x=x_axis, y=human_tes_64_xr_ts, line={'color':'#d62728'}, showlegend=False),
        go.Scatter(name='NTS', x=x_axis, y=human_tes_64_xr_nts, line={'color':'#1f76b4'}, showlegend=False)],
        rows=[2,2], cols=[2,2]
    )

    fig.add_traces(
        [go.Scatter(name='TS', x=x_axis, y=lemur_tss_64_xr_ts, line={'color':'#d62728'}, showlegend=False),
        go.Scatter(name='NTS', x=x_axis, y=lemur_tss_64_xr_nts, line={'color':'#1f76b4'}, showlegend=False)],
        rows=[2,2], cols=[3,3]
    )

    fig.add_traces(
        [go.Scatter(name='TS', x=x_axis, y=lemur_tes_64_xr_ts, line={'color':'#d62728'}, showlegend=False),
        go.Scatter(name='NTS', x=x_axis, y=lemur_tes_64_xr_nts, line={'color':'#1f76b4'}, showlegend=False)],
        rows=[2,2], cols=[4,4]
    )


    fig.update_xaxes(title_text="Position Relative to TSS (kb)", row=1, col=1)
    fig.update_xaxes(title_text="Position Relative to TSS (kb)", row=1, col=3)
    fig.update_xaxes(title_text="Position Relative to TSS (kb)", row=2, col=1)
    fig.update_xaxes(title_text="Position Relative to TSS (kb)", row=2, col=3)

    fig.update_xaxes(title_text="Position Relative to TES (kb)", row=1, col=2)
    fig.update_xaxes(title_text="Position Relative to TES (kb)", row=2, col=2)
    fig.update_xaxes(title_text="Position Relative to TES (kb)", row=1, col=4)
    fig.update_xaxes(title_text="Position Relative to TES (kb)", row=2, col=4)

    fig.update_yaxes(title_text="Relative Repair - CPD", row=1, col=1)
    fig.update_yaxes(title_text="Relative Repair - (6-4)PP", row=2, col=1)


    fig.update_layout(font={'family': 'Arial', 'size': 20}, template='plotly_white', height=1000, width=2250, margin=dict(
        l=90,
        r=0,
        b=70,
        t=10,
        pad=0
    ))
    fig.write_image(outfile)

if __name__ == "__main__":
    plot_tcr("tcr.pdf")
