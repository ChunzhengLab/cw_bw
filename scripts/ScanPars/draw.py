#!/usr/bin/env python3
"""
Script to draw Δδ (DelDelta) and Δγ (DelGamma) vs a chosen x-axis variable using ROOT.

Usage:
    python draw.py \
        --input-csv <csv_file> \
        --xaxis <xaxis>

Options:
    -i, --input-csv   Path to CSV file (output of extract_profiles.py)
    -x, --xaxis       X-axis variable: betaT, betaScale, rho2Scale, rho2_L, rho2_p, v2_p, v2_L

Notes:
    - If xaxis is betaT or betaScale, only rows with rho2Scale == 1 are used.
    - If xaxis is rho2Scale, rho2_L, rho2_p, v2_p, or v2_L, only rows with betaScale == 1 are used.
    - Each centrality (15,25,35,45,55) is drawn as a separate line in the legend.
    - Left pad: Δδ vs xaxis; right pad: Δγ vs xaxis.
"""
import argparse
import csv
from array import array
import ROOT

def parse_args():
    parser = argparse.ArgumentParser(
        description='Draw DelDelta and DelGamma vs chosen x-axis using ROOT'
    )
    parser.add_argument(
        '-i', '--input-csv',
        required=True,
        help='CSV file with processed results'
    )
    parser.add_argument(
        '-x', '--xaxis',
        required=True,
        choices=['betaT', 'betaScale', 'rho2Scale', 'rho2_L', 'rho2_p', 'v2_p', 'v2_L', 'v2', 'fLBC'],
        help='Variable to plot on the x-axis'
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Force scientific notation on y-axis labels globally
    ROOT.TGaxis.SetMaxDigits(4)

    # Read CSV
    data = []
    with open(args.input_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            # convert numeric fields
            cleaned = {
                'centrality': int(row['centrality']),
                'betaScale': float(row['betaScale']),
                'betaT': float(row['betaT']),
                'rho2Scale': float(row['rho2Scale']),
                'rho2_L': float(row['rho2_L']),
                'rho2_p': float(row['rho2_p']),
                'v2_p': float(row['v2_p']),
                'v2_p_Err': float(row['v2_p_Err']),
                'v2_L': float(row['v2_L']),
                'v2_L_Err': float(row['v2_L_Err']),
                'DelDelta': float(row['DelDelta']),
                'DelDelta_Err': float(row['DelDelta_Err']),
                'DelGamma': float(row['DelGamma']),
                'DelGamma_Err': float(row['DelGamma_Err']),
                'fLBC': float(row['fLBC']),
                'v2': float(row['v2']),
                'v2_Err': float(row['v2_Err']),
            }
            data.append(cleaned)

    # Map each centrality to its fixed fLBC value
    fracLBC_map = {
        15: 0.082,
        25: 0.079,
        35: 0.079,
        45: 0.074,
        55: 0.090,
    }

    # Apply filters based on x-axis selection
    if args.xaxis in ('betaT', 'betaScale'):
        # require rho2Scale == 1 and fLBC matches centrality
        data = [
            d for d in data
            if abs(d['rho2Scale'] - 1.0) < 1e-6
               and abs(d['fLBC'] - fracLBC_map[d['centrality']]) < 1e-6
        ]
    elif args.xaxis in ('rho2Scale', 'rho2_L', 'rho2_p', 'v2_p', 'v2_L', 'v2'):
        # require betaScale == 1 and fLBC matches centrality
        data = [
            d for d in data
            if abs(d['betaScale'] - 1.0) < 1e-6
               and abs(d['fLBC'] - fracLBC_map[d['centrality']]) < 1e-6
        ]
    elif args.xaxis == 'fLBC':
        # require betaScale == 1 and rho2Scale == 1
        data = [
            d for d in data
            if abs(d['betaScale'] - 1.0) < 1e-6
               and abs(d['rho2Scale'] - 1.0) < 1e-6
        ]
    else:
        data = data

    # Unique centrality values sorted
    cents = sorted({d['centrality'] for d in data})
    # Map centrality numeric values to descriptive labels
    cent_labels = {
        15: '10-20%',
        25: '20-30%',
        35: '30-40%',
        45: '40-50%',
        55: '50-60%',
    }
    # Determine dynamic x-axis limits based on data
    x_vals = [d[args.xaxis] for d in data]
    xmin, xmax = min(x_vals) * 0.8, max(x_vals) * 1.2

    # Setup custom colors using hex codes
    c15 = ROOT.TColor.GetColor("#183F5E")  # deep blue
    c25 = ROOT.TColor.GetColor("#21639A")  # light blue
    c35 = ROOT.TColor.GetColor("#F6D565")  # yellow
    c45 = ROOT.TColor.GetColor("#3EAEA4")  # green
    c55 = ROOT.TColor.GetColor("#EC553C")  # red
    colors = {
        15: c15,
        25: c25,
        35: c35,
        45: c45,
        55: c55,
    }

    # Create canvas with two pads
    c = ROOT.TCanvas('c', 'DelDelta & DelGamma', 1200, 600)
    c.Divide(2, 1)

    # Left pad: DelDelta
    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.12)
    graphs_delta = []
    legend1 = ROOT.TLegend(0.15, 0.65, 0.48, 0.85)
    legend1.SetNColumns(2)
    legend1.SetBorderSize(0)
    legend1.SetFillStyle(0)
    legend1.SetHeader("#bf{Centrality Interval}","C")
    for cent in cents:
        # select rows for this cent
        sub = sorted([d for d in data if d['centrality'] == cent], key=lambda x: x[args.xaxis])
        n = len(sub)
        x = array('d', [d[args.xaxis] for d in sub])
        y = array('d', [d['DelDelta'] for d in sub])
        yerr = array('d', [d['DelDelta_Err'] for d in sub])
        zero = array('d', [0.0]*n)

        g = ROOT.TGraphErrors(n, x, y, zero, yerr)
        g.SetLineColor(colors.get(cent, ROOT.kBlack))
        g.SetMarkerColor(colors.get(cent, ROOT.kBlack))
        g.SetMarkerStyle(20)
        g.SetLineWidth(2)
        if cent == cents[0]:
            # Determine axis ranges and labels for DelDelta
            if args.xaxis == 'fLBC':
                x_min, x_max, y_min, y_max = 0.0, 0.18, -0.001, 0.021
                x_label = '#it{f}_{LBC}'
            elif args.xaxis in ('betaT', 'betaScale'):
                x_min, x_max, y_min, y_max = 0.2, 1.00, -0.001, 0.017
                x_label = '#beta_{T}'
            elif args.xaxis == 'rho2Scale':
                x_min, x_max, y_min, y_max = 0.0, 0.28, -0.001, 0.017
                x_label = '#rho_{2}'
            elif args.xaxis == 'rho2_p':
                x_min, x_max, y_min, y_max = 0.0, 0.28, -0.001, 0.017
                x_label = '#rho_{2,p}'
            elif args.xaxis == 'rho2_L':
                x_min, x_max, y_min, y_max = 0.0, 0.28, -0.001, 0.017
                x_label = '#rho_{2,#Lambda}'
            elif args.xaxis == 'v2':
                x_min, x_max, y_min, y_max = 0.06, 0.22, -0.001, 0.017
                x_label = '#it{v}_{2}'
            elif args.xaxis == 'v2_p':
                x_min, x_max, y_min, y_max = 0.06, 0.22, -0.001, 0.017
                x_label = '#it{v}_{2,p}'
            elif args.xaxis == 'v2_L':
                x_min, x_max, y_min, y_max = 0.06, 0.22, -0.001, 0.017
                x_label = '#it{v}_{2,#Lambda}'
            else:
                x_min, x_max, y_min, y_max = xmin, xmax, y_min, y_max
                x_label = args.xaxis
            # Draw axis frame
            frame = ROOT.gPad.DrawFrame(x_min, y_min, x_max, y_max, f';{x_label};#Delta#delta')
            # Ensure y-axis uses exponent notation
            frame.GetYaxis().SetNoExponent(False)
            # Draw graph
            g.Draw('PL')
        else:
            g.Draw('PL SAME')
        legend1.AddEntry(g, cent_labels[cent], 'lp')
        graphs_delta.append(g)
    legend1.Draw()

    # Right pad: DelGamma
    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    graphs_gamma = []
    legend2 = ROOT.TLegend(0.15, 0.65, 0.48, 0.85)
    legend2.SetNColumns(2)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0)
    legend2.SetHeader("#bf{Centrality Interval}","C")
    for cent in cents:
        sub = sorted([d for d in data if d['centrality'] == cent], key=lambda x: x[args.xaxis])
        n = len(sub)
        x = array('d', [d[args.xaxis] for d in sub])
        y = array('d', [d['DelGamma'] for d in sub])
        yerr = array('d', [d['DelGamma_Err'] for d in sub])
        zero = array('d', [0.0]*n)

        g = ROOT.TGraphErrors(n, x, y, zero, yerr)
        g.SetLineColor(colors.get(cent, ROOT.kBlack))
        g.SetMarkerColor(colors.get(cent, ROOT.kBlack))
        g.SetMarkerStyle(21)
        g.SetLineWidth(2)
        if cent == cents[0]:
            # Determine axis ranges and labels for DelGamma
            if args.xaxis == 'fLBC':
                x_min, x_max, y_min, y_max = 0.0, 0.18, -0.0003, 0.0045
                x_label = '#it{f}_{LBC}'
            elif args.xaxis in ('betaT', 'betaScale'):
                x_min, x_max, y_min, y_max = 0.2, 1.00, -0.0002, 0.0040
                x_label = '#beta_{T}'
            elif args.xaxis == 'rho2Scale':
                x_min, x_max, y_min, y_max = 0.0, 0.28, -0.0003, 0.0040
                x_label = '#rho_{2}'
            elif args.xaxis == 'rho2_p':
                x_min, x_max, y_min, y_max = 0.0, 0.28, -0.0003, 0.0040
                x_label = '#rho_{2,p}'
            elif args.xaxis == 'rho2_L':
                x_min, x_max, y_min, y_max = 0.0, 0.28, -0.0003, 0.0040
                x_label = '#rho_{2,#Lambda}'
            elif args.xaxis == 'v2':
                x_min, x_max, y_min, y_max = 0.06, 0.22, -0.0002, 0.0040
                x_label = '#it{v}_{2}'
            elif args.xaxis == 'v2_p':
                x_min, x_max, y_min, y_max = 0.06, 0.22, -0.0002, 0.0040
                x_label = '#it{v}_{2,p}'
            elif args.xaxis == 'v2_L':
                x_min, x_max, y_min, y_max = 0.06, 0.22, -0.0002, 0.0040
                x_label = '#it{v}_{2,#Lambda}'
            else:
                x_min, x_max, y_min, y_max = xmin, xmax, y_min, y_max
                x_label = args.xaxis
            # Draw axis frame
            frame = ROOT.gPad.DrawFrame(x_min, y_min, x_max, y_max, f';{x_label};#Delta#gamma')
            # Ensure y-axis uses exponent notation
            frame.GetYaxis().SetNoExponent(False)
            # Draw graph
            g.Draw('PL')
        else:
            g.Draw('PL SAME')
        legend2.AddEntry(g, cent_labels[cent], 'lp')
        graphs_gamma.append(g)
    legend2.Draw()

    # Update and save
    c.Update()
    outname = f'draw_{args.xaxis}.pdf'
    c.SaveAs(outname)
    print(f'Saved plot to {outname}')

    # Keep the canvas open if not in batch mode
    if not ROOT.gROOT.IsBatch():
        c.Draw()
        ROOT.gApplication.Run()

if __name__ == '__main__':
    main()
