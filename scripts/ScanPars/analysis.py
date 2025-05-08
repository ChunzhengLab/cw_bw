#!/usr/bin/env python3
"""
Script to process ROOT TProfile objects from result files and output a CSV summary.
Usage:
    python extract_profiles.py \
        --root-dir <root_dir> \
        --config-yaml <config_yaml> \
        --output-csv <output_csv>

Options:
    -r, --root-dir      Directory containing ROOT files matching
                        results_cent<cent>_fLBC<...>[_rho2scale<rho>][_betascale<beta>].root
    -c, --config-yaml   Path to default.yaml with betaT, rho2_p, rho2_L arrays
    -o, --output-csv    Filename for the generated CSV report
"""
import os
import re
import yaml
import ROOT
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract TProfile observables from ROOT files and save to CSV'
    )
    parser.add_argument(
        '-r', '--root-dir',
        required=True,
        help='Directory containing ROOT result files'
    )
    parser.add_argument(
        '-c', '--config-yaml',
        required=True,
        help='Path to default.yaml with betaT, rho2_p, rho2_L arrays'
    )
    parser.add_argument(
        '-o', '--output-csv',
        required=True,
        help='Output CSV filename'
    )
    return parser.parse_args()


def main():
    args = parse_args()
    root_dir    = args.root_dir
    config_yaml = args.config_yaml
    output_csv  = args.output_csv

    # Load YAML configuration
    with open(config_yaml) as yf:
        config = yaml.safe_load(yf)
    betaT_array   = config.get('betaT', [])
    rho2_p_array  = config.get('rho2_p', [])
    rho2_L_array  = config.get('rho2_L', [])

    # Map centrality suffix to index in arrays
    cent_to_index = {'15': 0, '25': 1, '35': 2, '45': 3, '55': 4}

    # Regex to parse filenames with required fLBC and optional scales
    pattern = re.compile(
        r'results_cent(?P<cent>\d+)_fLBC(?P<fLBC>[\d.]+)'
        r'(?:_rho2scale(?P<rho2scale>[\d.]+))?'
        r'(?:_betascale(?P<betascale>[\d.]+))?\.root$'
    )

    rows = []

    for fname in os.listdir(root_dir):
        m = pattern.match(fname)
        if not m:
            continue
        cent_str = m.group('cent')
        if cent_str not in cent_to_index:
            print(f"Warning: unknown centrality '{cent_str}' in file {fname}, skipping.")
            continue

        idx        = cent_to_index[cent_str]
        centrality = int(cent_str)
        fLBC       = float(m.group('fLBC'))
        betaScale  = float(m.group('betascale') or 1.0)
        rho2Scale  = float(m.group('rho2scale') or 1.0)

        # Raw parameters from YAML
        raw_betaT  = betaT_array[idx]
        raw_rho2_p = rho2_p_array[idx]
        raw_rho2_L = rho2_L_array[idx]

        # Scaled parameters
        betaT   = raw_betaT  * betaScale
        rho2_p  = raw_rho2_p * rho2Scale
        rho2_L  = raw_rho2_L * rho2Scale

        # Open ROOT file
        path = os.path.join(root_dir, fname)
        f    = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            print(f"Error opening ROOT file: {path}")
            continue

        # Retrieve TProfiles
        profiles = {}
        observables = ['delta', 'gamma']
        pairs = [
            'antilambda_proton',
            'antiproton_lambda',
            'antilambda_antiproton',
            'proton_lambda'
        ]
        for obs in observables:
            for pair in pairs:
                key  = f"p_{obs}_{pair}"
                prof = f.Get(key)
                if not prof:
                    print(f"Warning: {key} not found in {fname}")
                else:
                    profiles[key] = prof

        # Compute SS, OS, and difference histograms
        results = {}
        for obs in observables:
            p_ss = profiles[f"p_{obs}_antilambda_antiproton"].Clone(f"p_{obs}_SS")
            p_ss.Add(profiles[f"p_{obs}_proton_lambda"])
            p_os = profiles[f"p_{obs}_antilambda_proton"].Clone(f"p_{obs}_OS")
            p_os.Add(profiles[f"p_{obs}_antiproton_lambda"])
            h_ss   = p_ss.ProjectionX(f"h_{obs}_SS")
            h_os   = p_os.ProjectionX(f"h_{obs}_OS")
            h_diff = h_os.Clone(f"h_{obs}_diff")
            h_diff.Add(h_ss, -1)

            bin_ss   = h_ss.FindBin(centrality)
            bin_os   = h_os.FindBin(centrality)
            bin_diff = h_diff.FindBin(centrality)

            content_ss, err_ss   = h_ss.GetBinContent(bin_ss),   h_ss.GetBinError(bin_ss)
            content_os, err_os   = h_os.GetBinContent(bin_os),   h_os.GetBinError(bin_os)
            content_d,  err_d    = h_diff.GetBinContent(bin_diff), h_diff.GetBinError(bin_diff)

            results[f"{obs.upper()}SS"]      = (content_ss, err_ss)
            results[f"{obs.upper()}OS"]      = (content_os, err_os)
            results[f"Del{obs.capitalize()}"] = (content_d,   err_d)

        # Retrieve v2 profiles and compute means
        prof_v2_p = f.Get("p_v2_pt_proton")
        prof_v2_L = f.Get("p_v2_pt_lambda")
        v2_p_mean, v2_p_err = (0, 0)
        if prof_v2_p:
            v2_p_mean = prof_v2_p.GetMean(2)
            v2_p_err = prof_v2_p.GetMeanError(2)
        else:
            print(f"Warning: p_v2_pt_proton not found in {fname}")
        v2_L_mean, v2_L_err = (0, 0)
        if prof_v2_L:
            v2_L_mean = prof_v2_L.GetMean(2)
            v2_L_err = prof_v2_L.GetMeanError(2)
        else:
            print(f"Warning: p_v2_pt_lambda not found in {fname}")
        results["v2_p"] = (v2_p_mean, v2_p_err)
        results["v2_L"] = (v2_L_mean, v2_L_err)

        # Compute combined v2
        v2_comb, v2_comb_err = (0, 0)
        if prof_v2_p and prof_v2_L:
            prof_v2_sum = prof_v2_p.Clone("p_v2_sum")
            prof_v2_sum.Add(prof_v2_L)
            v2_comb = prof_v2_sum.GetMean(2)
            v2_comb_err = prof_v2_sum.GetMeanError(2)
        else:
            print(f"Warning: Cannot compute combined v2 for {fname}")
        results["v2"]      = (v2_comb, v2_comb_err)

        f.Close()

        # Assemble row
        row = {
            'centrality':   centrality,
            'fLBC':         fLBC,
            'betaScale':    betaScale,
            'betaT':        betaT,
            'rho2Scale':    rho2Scale,
            'rho2_L':       rho2_L,
            'rho2_p':       rho2_p,
            'DeltaSS':      results['DELTASS'][0],
            'DeltaSS_Err':  results['DELTASS'][1],
            'DeltaOS':      results['DELTAOS'][0],
            'DeltaOS_Err':  results['DELTAOS'][1],
            'GammaSS':      results['GAMMASS'][0],
            'GammaSS_Err':  results['GAMMASS'][1],
            'GammaOS':      results['GAMMAOS'][0],
            'GammaOS_Err':  results['GAMMAOS'][1],
            'DelDelta':     results['DelDelta'][0],
            'DelDelta_Err': results['DelDelta'][1],
            'DelGamma':     results['DelGamma'][0],
            'DelGamma_Err': results['DelGamma'][1],
            'v2_p':         results['v2_p'][0],
            'v2_p_Err':     results['v2_p'][1],
            'v2_L':         results['v2_L'][0],
            'v2_L_Err':     results['v2_L'][1],
            'v2':           results['v2'][0],
            'v2_Err':       results['v2'][1],
        }
        rows.append(row)

    # Write CSV
    import csv
    headers = [
        'centrality','fLBC','betaScale','betaT','rho2Scale','rho2_L','rho2_p',
        'DeltaSS','DeltaSS_Err','DeltaOS','DeltaOS_Err',
        'GammaSS','GammaSS_Err','GammaOS','GammaOS_Err',
        'DelDelta','DelDelta_Err','DelGamma','DelGamma_Err',
        'v2_p','v2_p_Err','v2_L','v2_L_Err','v2','v2_Err'
    ]
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"Wrote summary for {len(rows)} files to {output_csv}")

if __name__ == '__main__':
    main()
