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
        choices=['betaT', 'betaScale', 'rho2Scale', 'rho2_L', 'rho2_p', 'v2_p', 'v2_L'],
        help='Variable to plot on the x-axis'
    )
    return parser.parse_args()


def main():
    args = parse_args()

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
            }
            data.append(cleaned)

    # Apply filters
    if args.xaxis in ('betaT', 'betaScale'):
        data = [d for d in data if abs(d['rho2Scale'] - 1.0) < 1e-6]
    else:
        data = [d for d in data if abs(d['betaScale'] - 1.0) < 1e-6]

    # Unique centrality values sorted
    cents = sorted({d['centrality'] for d in data})
    # Determine dynamic x-axis limits based on data
    x_vals = [d[args.xaxis] for d in data]
    xmin, xmax = min(x_vals) * 0.8, max(x_vals) * 1.2

    # Setup colors for each centrality
    colors = {
        15: ROOT.kRed,
        25: ROOT.kBlue,
        35: ROOT.kGreen+2,
        45: ROOT.kMagenta,
        55: ROOT.kCyan
    }

    # Create canvas with two pads
    c = ROOT.TCanvas('c', 'DelDelta & DelGamma', 1200, 600)
    c.Divide(2, 1)

    # Left pad: DelDelta
    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.12)
    graphs_delta = []
    legend1 = ROOT.TLegend(0.65, 0.7, 0.9, 0.9)
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
            g.SetTitle(f'#Delta#delta vs {args.xaxis}')
            g.GetXaxis().SetTitle(args.xaxis)
            g.GetYaxis().SetTitle('#Delta#delta')
            g.GetYaxis().SetRangeUser(0, 1)
            # Set dynamic x-axis limits
            g.GetXaxis().SetLimits(xmin, xmax)
            g.Draw('APL')
        else:
            g.Draw('PL SAME')
        legend1.AddEntry(g, f'cent={cent}%', 'lp')
        graphs_delta.append(g)
    legend1.Draw()

    # Right pad: DelGamma
    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    graphs_gamma = []
    legend2 = ROOT.TLegend(0.65, 0.7, 0.9, 0.9)
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
            g.SetTitle(f'#Delta#gamma vs {args.xaxis}')
            g.GetXaxis().SetTitle(args.xaxis)
            g.GetYaxis().SetTitle('#Delta#gamma')
            g.GetYaxis().SetRangeUser(0, 0.1)
            # Set dynamic x-axis limits
            g.GetXaxis().SetLimits(xmin, xmax)
            g.Draw('APL')
        else:
            g.Draw('PL SAME')
        legend2.AddEntry(g, f'cent={cent}%', 'lp')
        graphs_gamma.append(g)
    legend2.Draw()

    # Update and save
    c.Update()
    outname = f'draw_{args.xaxis}.png'
    c.SaveAs(outname)
    print(f'Saved plot to {outname}')

    # Keep the canvas open if not in batch mode
    if not ROOT.gROOT.IsBatch():
        c.Draw()
        ROOT.gApplication.Run()

if __name__ == '__main__':
    main()