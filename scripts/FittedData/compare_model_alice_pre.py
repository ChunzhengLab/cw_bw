#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
compare.py

Builds Δδ and Δγ vs. centrality, comparing ALICE‐preliminary results from
cve_alice_preliminary.csv with model outputs contained in
results_cent*_fLBC*.root files.

The script:
1. Reads Δ_OS, Δ_SS, Γ_OS, Γ_SS (and their stat & syst errors) from the CSV,
   computing total uncertainties and the OS – SS differences.
2. Scans the ../build directory for ROOT files matching
   “results_cent<cent>_fLBC<val>.root”; for every file it:
      • combines the relevant TProfiles into OS and SS profiles,
      • projects to TH1D, forms OS – SS histograms, and extracts their means.
3. Collects these numbers for each fLBC setting across all centralities.
4. Plots two side‑by‑side pads:
      • left: Δδ vs. centrality;
      • right: Δγ vs. centrality,
   overlaying ALICE points (black) with curves for every fLBC (coloured).
5. Saves the canvas as deltaDelta_deltaGamma_comparison.pdf.
"""

import ROOT # type: ignore
import csv
import math
import glob
import os
import re
import argparse

ROOT.gROOT.SetBatch(True)
DEBUG = True  # set to False to silence debug prints


# ----------------------------------------------------------------------
# 1) ALICE preliminary
# ----------------------------------------------------------------------
def read_alice_preliminary(csv_path):
    cent, dDelta, dDeltaErr, dGamma, dGammaErr = [], [], [], [], []

    with open(csv_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            c = float(row["Centrality (%)"])

            # Δ values
            dss      = float(row["Delta_SS"])
            dss_stat = float(row["Delta_SS_StatErr"])
            dss_sys  = float(row["Delta_SS_SystErr"])
            dos      = float(row["Delta_OS"])
            dos_stat = float(row["Delta_OS_StatErr"])
            dos_sys  = float(row["Delta_OS_SystErr"])

            # Γ values
            gss      = float(row["Gamma_SS"])
            gss_stat = float(row["Gamma_SS_StatErr"])
            gss_sys  = float(row["Gamma_SS_SystErr"])
            gos      = float(row["Gamma_OS"])
            gos_stat = float(row["Gamma_OS_StatErr"])
            gos_sys  = float(row["Gamma_OS_SystErr"])

            # total uncertainties
            err_dss = math.hypot(dss_stat, dss_sys)
            err_dos = math.hypot(dos_stat, dos_sys)
            err_gss = math.hypot(gss_stat, gss_sys)
            err_gos = math.hypot(gos_stat, gos_sys)

            # OS – SS
            delta_delta = dos - dss
            delta_gamma = gos - gss
            err_delta   = math.hypot(err_dos, err_dss)
            err_gamma   = math.hypot(err_gos, err_gss)

            cent.append(c)
            dDelta.append(delta_delta)
            dDeltaErr.append(err_delta)
            dGamma.append(delta_gamma)
            dGammaErr.append(err_gamma)

    return cent, dDelta, dDeltaErr, dGamma, dGammaErr


# ----------------------------------------------------------------------
# 2) Model ROOT files
# ----------------------------------------------------------------------
def process_root_files(root_dir: str):
    pattern = os.path.join(root_dir, "results_cent*_fLBC*.root")
    files  = glob.glob(pattern)
    if DEBUG:
        print(f"[DEBUG] Found {len(files)} ROOT files matching pattern '{pattern}'")
        for f in files:
            print(f"[DEBUG]  └─ {f}")
    regex  = re.compile(r"results_cent(\d+)_fLBC([0-9.]+)\.root")
    result = {}

    for path in files:
        if DEBUG:
            print(f"\n[DEBUG] Processing file: {path}")
        m = regex.search(os.path.basename(path))
        if not m:
            continue

        cent  = int(m.group(1))
        flbc  = float(m.group(2))
        if DEBUG:
            print(f"[DEBUG]  • Parsed centrality={cent}, fLBC={flbc}")
        froot = ROOT.TFile.Open(path)
        if not froot or froot.IsZombie():
            print("⚠️  Could not open", path)
            continue

        # Grab all eight profiles
        names = [
            "p_delta_antiproton_lambda",
            "p_delta_antilambda_proton",
            "p_delta_antilambda_antiproton",
            "p_delta_proton_lambda",
            "p_gamma_antiproton_lambda",
            "p_gamma_antilambda_proton",
            "p_gamma_antilambda_antiproton",
            "p_gamma_proton_lambda",
        ]
        prof = {n: froot.Get(n) for n in names}
        if DEBUG:
            available = [k for k, v in prof.items() if v]
            print("[DEBUG]  • Retrieved TProfiles:", available)
        if any(p is None for p in prof.values()):
            print("⚠️  Missing TProfiles in", path)
            if DEBUG:
                missing = [k for k, v in prof.items() if v is None]
                print(f"[DEBUG]  × Missing profiles: {missing}")
            froot.Close()
            continue

        # Build OS and SS profiles ------------------------------
        p_d_SS = prof["p_delta_antilambda_antiproton"].Clone(f"p_d_SS_{cent}_{flbc}")
        p_d_SS.Add(prof["p_delta_proton_lambda"])
        p_d_OS = prof["p_delta_antiproton_lambda"].Clone(f"p_d_OS_{cent}_{flbc}")
        p_d_OS.Add(prof["p_delta_antilambda_proton"])

        p_g_SS = prof["p_gamma_antilambda_antiproton"].Clone(f"p_g_SS_{cent}_{flbc}")
        p_g_SS.Add(prof["p_gamma_proton_lambda"])
        p_g_OS = prof["p_gamma_antiproton_lambda"].Clone(f"p_g_OS_{cent}_{flbc}")
        p_g_OS.Add(prof["p_gamma_antilambda_proton"])

        # Project to TH1D ---------------------------------------
        h_d_SS = p_d_SS.ProjectionX(f"h_d_SS_{cent}_{flbc}")
        h_d_OS = p_d_OS.ProjectionX(f"h_d_OS_{cent}_{flbc}")
        h_g_SS = p_g_SS.ProjectionX(f"h_g_SS_{cent}_{flbc}")
        h_g_OS = p_g_OS.ProjectionX(f"h_g_OS_{cent}_{flbc}")

        h_dDelta = h_d_OS.Clone(f"h_dDelta_{cent}_{flbc}")
        h_dDelta.Add(h_d_SS, -1)

        h_dGamma = h_g_OS.Clone(f"h_dGamma_{cent}_{flbc}")
        h_dGamma.Add(h_g_SS, -1)

        # Use the bin corresponding to the current centrality value
        xaxis = h_dDelta.GetXaxis()
        bin_idx = xaxis.FindBin(cent)
        val_dDelta = h_dDelta.GetBinContent(bin_idx)
        err_dDelta = h_dDelta.GetBinError(bin_idx)

        xaxis_g = h_dGamma.GetXaxis()
        bin_idx_g = xaxis_g.FindBin(cent)
        val_dGamma = h_dGamma.GetBinContent(bin_idx_g)
        err_dGamma = h_dGamma.GetBinError(bin_idx_g)

        if DEBUG:
            print(f"[DEBUG]  • Using bin {bin_idx} (cent={cent}) for Δδ and bin {bin_idx_g} for Δγ")

        if flbc not in result:
            result[flbc] = {"cent": [], "dDelta": [], "dDeltaErr": [], "dGamma": [], "dGammaErr": []}

        result[flbc]["cent"].append(cent)
        result[flbc]["dDelta"].append(val_dDelta)
        result[flbc]["dDeltaErr"].append(err_dDelta)
        result[flbc]["dGamma"].append(val_dGamma)
        result[flbc]["dGammaErr"].append(err_dGamma)

        if DEBUG:
            print(f"[DEBUG]  • Δδ(mean)={val_dDelta:.4e} ± {err_dDelta:.4e}, "
                  f"Δγ(mean)={val_dGamma:.4e} ± {err_dGamma:.4e}")
        froot.Close()

    # sort by centrality
    for flbc, d in result.items():
        z = sorted(zip(d["cent"], d["dDelta"], d["dDeltaErr"], d["dGamma"], d["dGammaErr"]))
        d["cent"], d["dDelta"], d["dDeltaErr"], d["dGamma"], d["dGammaErr"] = map(list, zip(*z))

    if DEBUG:
        print("\n[DEBUG] Summary of collected points:")
        for f, d in result.items():
            print(f"  fLBC={f} ⇒ centralities={d['cent']}")
    return result


# ----------------------------------------------------------------------
# 3) Helpers
# ----------------------------------------------------------------------
def make_graph(x, y, yerr, colour, marker=20, name="gr"):
    g = ROOT.TGraphErrors(len(x))
    g.SetName(name)
    for i in range(len(x)):
        g.SetPoint(i, x[i], y[i])
        g.SetPointError(i, 0, yerr[i])
    g.SetMarkerStyle(marker)
    g.SetMarkerColor(colour)
    g.SetLineColor(colour)
    return g


# ----------------------------------------------------------------------
# 4) Plot
# ----------------------------------------------------------------------
def plot_comparison(al_cent, al_dD, al_dDE, al_dG, al_dGE, model):
    canvas = ROOT.TCanvas("c", "Δδ / Δγ comparison", 1200, 600)
    palette = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2,
               ROOT.kMagenta + 2, ROOT.kOrange + 7, ROOT.kCyan + 1,
               ROOT.kAzure + 7]
    canvas.Divide(2,1);

    # -- Δδ ----------------------------------
    canvas.cd(1)
    mg_dD = ROOT.TMultiGraph()
    leg1  = ROOT.TLegend(0.12, 0.72, 0.42, 0.88)

    g_al_dD = make_graph(al_cent, al_dD, al_dDE, ROOT.kBlack, 21, "g_alice_dD")
    mg_dD.Add(g_al_dD, "P")
    leg1.AddEntry(g_al_dD, "ALICE preliminary", "p")

    for idx, (flbc, d) in enumerate(sorted(model.items())):
        col = palette[idx % len(palette)]
        g = make_graph(d["cent"], d["dDelta"], d["dDeltaErr"], col,
                       20, f"g_dD_{flbc}")
        mg_dD.Add(g, "P")
        leg1.AddEntry(g, f"f_{{LBC}} = {flbc:.3f}", "p")

    mg_dD.Draw("A P")
    mg_dD.GetHistogram().GetXaxis().SetLimits(0, 60)
    mg_dD.GetHistogram().GetYaxis().SetRangeUser(-0.002, 0.015)
    mg_dD.GetXaxis().SetTitle("Centrality (%)")
    mg_dD.GetYaxis().SetTitle("#delta#delta")
    leg1.Draw()

    canvas.cd(2)
    mg_dG = ROOT.TMultiGraph()
    leg2  = ROOT.TLegend(0.12, 0.72, 0.42, 0.88)

    g_al_dG = make_graph(al_cent, al_dG, al_dGE, ROOT.kBlack, 21, "g_alice_dG")
    mg_dG.Add(g_al_dG, "P")

    for idx, (flbc, d) in enumerate(sorted(model.items())):
        col = palette[idx % len(palette)]
        g = make_graph(d["cent"], d["dGamma"], d["dGammaErr"], col,
                       20, f"g_dG_{flbc}")
        mg_dG.Add(g, "P")
        leg2.AddEntry(g, f"f_{{LBC}} = {flbc:.3f}", "p")

    mg_dG.Draw("A P")
    mg_dG.GetHistogram().GetXaxis().SetLimits(0, 60)
    mg_dG.GetHistogram().GetYaxis().SetRangeUser(-0.0003, 0.003)
    mg_dG.GetXaxis().SetTitle("Centrality (%)")
    mg_dG.GetYaxis().SetTitle("#delta#gamma")
    leg2.Draw()

    canvas.SaveAs("deltaDelta_deltaGamma_comparison.pdf")
    print("✅  Saved: deltaDelta_deltaGamma_comparison.pdf")

# ----------------------------------------------------------------------
# 6) Plot dDelta/dGamma vs fLBC for each centrality
# ----------------------------------------------------------------------
def plot_dDeltaGamma_vs_fLBC(model, al_cent, al_dD, al_dDE, al_dG, al_dGE, filename="dDelta_dGamma_vs_fLBC.pdf"):
    import collections
    from collections import defaultdict

    # Reorganize by centrality
    cent_dict_dD = defaultdict(list)
    cent_dict_dG = defaultdict(list)

    for flbc, data in model.items():
        for c, val_dD, err_dD, val_dG, err_dG in zip(data["cent"], data["dDelta"], data["dDeltaErr"],
                                                     data["dGamma"], data["dGammaErr"]):
            cent_dict_dD[c].append((flbc, val_dD, err_dD))
            cent_dict_dG[c].append((flbc, val_dG, err_dG))

    canvas = ROOT.TCanvas("c_fLBC", "Δδ and Δγ vs fLBC", 1200, 600)
    canvas.Divide(2, 1)

    palette = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2,
               ROOT.kMagenta + 2, ROOT.kOrange + 7, ROOT.kCyan + 1,
               ROOT.kAzure + 7]

    # Δδ pad
    canvas.cd(1)
    ROOT.gPad.SetGrid()
    mg_dD = ROOT.TMultiGraph()
    leg1 = ROOT.TLegend(0.12, 0.68, 0.45, 0.88)

    # --- Initialize results dictionaries ---
    results_dD = {}
    results_dG = {}

    for idx, (cent, lst) in enumerate(sorted(cent_dict_dD.items())):
        flbc_vals, vals, errs = zip(*sorted(lst))
        col = palette[idx % len(palette)]
        gr = make_graph(flbc_vals, vals, errs, col, 20, f"gr_dD_{cent}")
        mg_dD.Add(gr, "LP")
        leg1.AddEntry(gr, f"Centrality = {cent}%", "p")

        # --- Polynomial (order 2) fit and marker at ALICE value ---
        f1 = ROOT.TF1(f"fit_dD_{cent}", "pol2", 0.0, 1.5)
        f1.SetLineColor(col)
        f1.SetLineWidth(2)
        gr.Fit(f1, "Q")
        alice_val_dD = al_dD[al_cent.index(cent)] if cent in al_cent else None
        if alice_val_dD is not None:
            # Find root(s) of quadratic: f1(x) = alice_val_dD
            a = f1.GetParameter(2)
            b = f1.GetParameter(1)
            c = f1.GetParameter(0) - alice_val_dD
            x0 = float('nan')
            if abs(a) > 1e-12:
                disc = b*b - 4*a*c
                if disc >= 0:
                    sqrt_disc = math.sqrt(disc)
                    sol1 = (-b + sqrt_disc) / (2*a)
                    sol2 = (-b - sqrt_disc) / (2*a)
                    # Choose solution in [0, 1.5]
                    candidates = [sol for sol in (sol1, sol2) if 0.0 <= sol <= 1.5]
                    x0 = candidates[0] if candidates else float('nan')
            elif abs(b) > 1e-12:
                x0 = -c / b
            # Store in results_dD and draw marker
            results_dD[cent] = x0
            tm = ROOT.TMarker(x0, alice_val_dD, 29)
            tm.SetMarkerColor(col)
            tm.SetMarkerSize(1.2)
            tm.Draw()


    mg_dD.Draw("A LP")
    mg_dD.GetHistogram().GetYaxis().SetRangeUser(-0.002, 0.015)
    mg_dD.GetHistogram().GetXaxis().SetTitle("f_{LBC}")
    mg_dD.GetHistogram().GetYaxis().SetTitle("#Delta#delta")
    leg1.Draw()

    # Δγ pad
    canvas.cd(2)
    ROOT.gPad.SetGrid()
    mg_dG = ROOT.TMultiGraph()
    leg2 = ROOT.TLegend(0.12, 0.68, 0.45, 0.88)

    for idx, (cent, lst) in enumerate(sorted(cent_dict_dG.items())):
        flbc_vals, vals, errs = zip(*sorted(lst))
        col = palette[idx % len(palette)]
        gr = make_graph(flbc_vals, vals, errs, col, 24, f"gr_dG_{cent}")
        mg_dG.Add(gr, "LP")
        leg2.AddEntry(gr, f"Centrality = {cent}%", "p")

        # --- Polynomial (order 2) fit and marker at ALICE value ---
        f2 = ROOT.TF1(f"fit_dG_{cent}", "pol2", 0.0, 1.5)
        f2.SetLineColor(col)
        f2.SetLineWidth(2)
        gr.Fit(f2, "Q")
        alice_val_dG = al_dG[al_cent.index(cent)] if cent in al_cent else None
        if alice_val_dG is not None:
            # Find root(s) of quadratic: f2(x) = alice_val_dG
            a = f2.GetParameter(2)
            b = f2.GetParameter(1)
            c = f2.GetParameter(0) - alice_val_dG
            x1 = float('nan')
            if abs(a) > 1e-12:
                disc = b*b - 4*a*c
                if disc >= 0:
                    sqrt_disc = math.sqrt(disc)
                    sol1 = (-b + sqrt_disc) / (2*a)
                    sol2 = (-b - sqrt_disc) / (2*a)
                    candidates = [sol for sol in (sol1, sol2) if 0.0 <= sol <= 1.5]
                    x1 = candidates[0] if candidates else float('nan')
            elif abs(b) > 1e-12:
                x1 = -c / b
            # Store in results_dG and draw marker
            results_dG[cent] = x1
            tm2 = ROOT.TMarker(x1, alice_val_dG, 29)
            tm2.SetMarkerColor(col)
            tm2.SetMarkerSize(1.2)
            tm2.Draw()


    mg_dG.Draw("SAME LP")
    mg_dG.GetHistogram().GetYaxis().SetRangeUser(-0.001, 0.003)
    mg_dG.GetHistogram().GetXaxis().SetTitle("f_{LBC}")
    mg_dG.GetHistogram().GetYaxis().SetTitle("#Delta#gamma")
    leg2.Draw()

    # --- Print combined results for all centralities ---
    print("fLBC fit results by centrality:")
    for cent in sorted(results_dD.keys()):
        d0 = results_dD.get(cent, float('nan'))
        g0 = results_dG.get(cent, float('nan'))
        print(f" cent={cent}% ⇒ fLBC_Δδ={d0:.4f}, fLBC_Δγ={g0:.4f}")

    canvas.SaveAs(filename)
    print(f"✅  Saved: {filename}")


# ----------------------------------------------------------------------
# 5) Argument parsing
# ----------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Compare ALICE preliminary data with model outputs.")
    parser.add_argument("--root-dir", default=".", help="Directory containing ROOT files.")
    parser.add_argument("--csv-alice", default="../../refdata/cve_alice_preliminary.csv", help="CSV file path for ALICE preliminary data.")
    parser.add_argument("--output-dir", default=".", help="Directory to save output plots.")
    return parser.parse_args()


# ----------------------------------------------------------------------
# 6) Main
# ----------------------------------------------------------------------
def main():
    args = parse_args()
    alice = read_alice_preliminary(args.csv_alice)
    model = process_root_files(args.root_dir)

    os.makedirs(args.output_dir, exist_ok=True)
    comparison_plot = os.path.join(args.output_dir, "deltaDelta_deltaGamma_comparison.pdf")
    flbc_plot = os.path.join(args.output_dir, "dDelta_dGamma_vs_fLBC.pdf")

    plot_comparison(*alice, model)
    os.rename("deltaDelta_deltaGamma_comparison.pdf", comparison_plot)

    plot_dDeltaGamma_vs_fLBC(model, *alice, flbc_plot)


if __name__ == "__main__":
    main()
