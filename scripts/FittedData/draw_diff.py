import ROOT
import sys
import csv
import math

def convert_profile(p):
    nbins = p.GetNbinsX()
    xmin = p.GetXaxis().GetXmin()
    xmax = p.GetXaxis().GetXmax()
    h = ROOT.TH1D(f"{p.GetName()}_clone", "", nbins, xmin, xmax)
    for b in range(1, nbins + 1):
        h.SetBinContent(b, p.GetBinContent(b))
        h.SetBinError(b, p.GetBinError(b))
    return h

def read_alice_preliminary(csv_file):
    cent = []
    dDelta = []
    dDeltaErr = []
    dGamma = []
    dGammaErr = []
    with open(csv_file, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Parse centrality
            c = float(row["Centrality (%)"])
            # Parse Delta SS and OS and their stat/systematic errors
            dss = float(row["Delta_SS"])
            dss_stat = float(row["Delta_SS_StatErr"])
            dss_sys  = float(row["Delta_SS_SystErr"])
            dos = float(row["Delta_OS"])
            dos_stat = float(row["Delta_OS_StatErr"])
            dos_sys  = float(row["Delta_OS_SystErr"])
            # Parse Gamma SS and OS and their errors
            gss = float(row["Gamma_SS"])
            gss_stat = float(row["Gamma_SS_StatErr"])
            gss_sys  = float(row["Gamma_SS_SystErr"])
            gos = float(row["Gamma_OS"])
            gos_stat = float(row["Gamma_OS_StatErr"])
            gos_sys  = float(row["Gamma_OS_SystErr"])
            # Total uncertainties
            err_dss = math.hypot(dss_stat, dss_sys)
            err_dos = math.hypot(dos_stat, dos_sys)
            err_gss = math.hypot(gss_stat, gss_sys)
            err_gos = math.hypot(gos_stat, gos_sys)
            # Compute OS - SS and combined errors
            delta_delta = dos - dss
            delta_gamma = gos - gss
            err_delta = math.hypot(err_dos, err_dss)
            err_gamma = math.hypot(err_gos, err_gss)
            # Append to lists
            cent.append(c)
            dDelta.append(delta_delta)
            dDeltaErr.append(err_delta)
            dGamma.append(delta_gamma)
            dGammaErr.append(err_gamma)
    return cent, dDelta, dDeltaErr, dGamma, dGammaErr

def draw_diff(filename="results_combine.root"):
    import argparse
    parser = argparse.ArgumentParser(description="Plot diff or intg mode")
    parser.add_argument('--mode', choices=['diff','intg'], default='diff')
    parser.add_argument('--alice-csv', help='Path to ALICE preliminary CSV')
    args = parser.parse_args()
    mode = args.mode
    alice_csv = args.alice_csv

    f = ROOT.TFile.Open(filename)
    if not f or f.IsZombie():
        print(f"Cannot open file {filename}")
        return

    # Define centrality points for graphs
    centroids = [15, 25, 35, 45, 55]

    if mode == 'diff':
        types = ["delta", "gamma"]
        pairsSS = ["antilambda_antiproton", "proton_lambda"]
        pairsOS = ["antiproton_lambda", "antilambda_proton"]

        nSumPt = 3
        nEtaGap = 4

        # Use ROOT built-in color numbers
        colorsEtaGap = [ROOT.kBlue + 2, ROOT.kRed + 1, ROOT.kOrange + 1, ROOT.kCyan + 1]
        colorsSumPt = colorsEtaGap[:3]

        sumPtLabels = ["1<#Sigma p_{T}<3 GeV/c", "3<#Sigma p_{T}<5 GeV/c", "5<#Sigma p_{T}<8 GeV/c"]
        etaGapLabels = ["0.2<#Delta#eta<0.4", "0.4<#Delta#eta<0.6", "0.6<#Delta#eta<0.8", "#Delta#eta>0.8"]

        h_diff_sumPt = {}
        h_diff_etaGap = {}
        h_diff_int = {}

        # Keep legend objects alive to prevent garbage collection
        legends = []
        graphs_sumPt = []
        graphs_etaGap = []

        for type_ in types:
            for i in range(nSumPt):
                p_ss1 = f.Get(f"p_{type_}_{pairsSS[0]}_sumPt_{i}")
                p_ss2 = f.Get(f"p_{type_}_{pairsSS[1]}_sumPt_{i}")
                if not p_ss1 or not p_ss2: continue
                p_ss1.Add(p_ss2)
                h_ss = convert_profile(p_ss1)

                p_os1 = f.Get(f"p_{type_}_{pairsOS[0]}_sumPt_{i}")
                p_os2 = f.Get(f"p_{type_}_{pairsOS[1]}_sumPt_{i}")
                if not p_os1 or not p_os2: continue
                p_os1.Add(p_os2)
                h_os = convert_profile(p_os1)

                h_diff = h_os.Clone(f"h_{type_}_sumPt_{i}")
                h_diff.Add(h_ss, -1.0)
                h_diff_sumPt[f"{type_}_{i}"] = h_diff

            for j in range(nEtaGap):
                p_ss1 = f.Get(f"p_{type_}_{pairsSS[0]}_etaGap_{j}")
                p_ss2 = f.Get(f"p_{type_}_{pairsSS[1]}_etaGap_{j}")
                if not p_ss1 or not p_ss2: continue
                p_ss1.Add(p_ss2)
                h_ss = convert_profile(p_ss1)

                p_os1 = f.Get(f"p_{type_}_{pairsOS[0]}_etaGap_{j}")
                p_os2 = f.Get(f"p_{type_}_{pairsOS[1]}_etaGap_{j}")
                if not p_os1 or not p_os2: continue
                p_os1.Add(p_os2)
                h_os = convert_profile(p_os1)

                h_diff = h_os.Clone(f"h_{type_}_etaGap_{j}")
                h_diff.Add(h_ss, -1.0)
                h_diff_etaGap[f"{type_}_{j}"] = h_diff

            # integrated
            p_ssi1 = f.Get(f"p_{type_}_{pairsSS[0]}")
            p_ssi2 = f.Get(f"p_{type_}_{pairsSS[1]}")
            if not p_ssi1 or not p_ssi2: continue
            p_ssi1.Add(p_ssi2)
            h_ssi = convert_profile(p_ssi1)

            p_osi1 = f.Get(f"p_{type_}_{pairsOS[0]}")
            p_osi2 = f.Get(f"p_{type_}_{pairsOS[1]}")
            if not p_osi1 or not p_osi2: continue
            p_osi1.Add(p_osi2)
            h_osi = convert_profile(p_osi1)

            h_int = h_osi.Clone(f"h_int_{type_}")
            h_int.Add(h_ssi, -1.0)
            h_int.SetFillColorAlpha(ROOT.kRed + 2, 0.3)
            h_int.SetLineColor(ROOT.kRed + 2)
            h_int.SetLineStyle(1)
            h_int.SetLineWidth(2)
            h_int.SetMarkerColor(ROOT.kRed + 2)
            h_diff_int[type_] = h_int

        # === Plot SumPt dependence ===
        c1 = ROOT.TCanvas("c1", "SumPt dependence", 1200, 500)
        c1.Divide(2,1)

        for idx, type_ in enumerate(types):
            pad = c1.cd(idx+1)
            pad.SetGrid()
            if type_ == "delta":
                frame = pad.DrawFrame(10., -0.002, 60., 0.02, ";centrality (%);#Delta#delta")
            else:
                frame = pad.DrawFrame(10., -0.0010, 60., 0.010, ";centrality (%);#Delta#gamma")
            frame.GetXaxis().SetNdivisions(5, ROOT.kFALSE)
            frame.GetXaxis().SetTickLength(0.03)
            leg = ROOT.TLegend(0.2,0.6,0.4,0.8)
            legends.append(leg)
            leg.SetHeader("#Sigma p_{T} selection", "C")
            if type_ in h_diff_int:
                # Convert integrated histogram to TGraphErrors in black
                h_int = h_diff_int[type_]
                g_int = ROOT.TGraphErrors(len(centroids))
                for idx_c, cent in enumerate(centroids):
                    bin_idx = h_int.FindBin(cent)
                    y = h_int.GetBinContent(bin_idx)
                    err = h_int.GetBinError(bin_idx)
                    g_int.SetPoint(idx_c, cent, y)
                    g_int.SetPointError(idx_c, 0, err)
                g_int.SetLineColor(ROOT.kBlack)
                g_int.SetMarkerColor(ROOT.kBlack)
                g_int.SetLineStyle(1)
                g_int.SetLineWidth(2)
                g_int.SetMarkerStyle(21)
                g_int.Draw("LP SAME")
                legends.append(g_int)  # keep alive
                leg.AddEntry(g_int, "Integrated", "lp")
            for i in range(nSumPt):
                key = f"{type_}_{i}"
                if key not in h_diff_sumPt: continue
                h_hist = h_diff_sumPt[key]
                # Convert histogram to TGraphErrors at fixed centroids
                graph = ROOT.TGraphErrors(len(centroids))
                for idx_c, cent in enumerate(centroids):
                    bin_idx = h_hist.FindBin(cent)
                    y = h_hist.GetBinContent(bin_idx)
                    err = h_hist.GetBinError(bin_idx)
                    graph.SetPoint(idx_c, cent, y)
                    graph.SetPointError(idx_c, 0, err)
                graph.SetLineColor(colorsSumPt[i])
                graph.SetMarkerColor(colorsSumPt[i])
                graph.SetLineStyle(1)
                graph.SetLineWidth(2)
                graph.SetMarkerStyle(20)
                graph.SetMarkerSize(1.0)
                graph.Draw("LP SAME")
                graphs_sumPt.append(graph)
                leg.AddEntry(graph, sumPtLabels[i], "p")
            leg.Draw()
            pad.Modified()
            pad.Update()

        # === Plot EtaGap dependence ===
        c2 = ROOT.TCanvas("c2", "EtaGap dependence", 1200, 500)
        c2.Divide(2,1)

        for idx, type_ in enumerate(types):
            pad = c2.cd(idx+1)
            pad.SetGrid()
            if type_ == "delta":
                frame = pad.DrawFrame(10., -0.002, 60., 0.02, ";centrality (%);#Delta#delta")
            else:
                frame = pad.DrawFrame(10., -0.00035, 60., 0.0035, ";centrality (%);#Delta#gamma")
            frame.GetXaxis().SetNdivisions(5, ROOT.kFALSE)
            frame.GetXaxis().SetTickLength(0.03)
            leg = ROOT.TLegend(0.2,0.6,0.4,0.8)
            legends.append(leg)
            leg.SetHeader("#Delta#eta selection", "C")
            if type_ in h_diff_int:
                # Convert integrated histogram to TGraphErrors in black
                h_int = h_diff_int[type_]
                g_int = ROOT.TGraphErrors(len(centroids))
                for idx_c, cent in enumerate(centroids):
                    bin_idx = h_int.FindBin(cent)
                    y = h_int.GetBinContent(bin_idx)
                    err = h_int.GetBinError(bin_idx)
                    g_int.SetPoint(idx_c, cent, y)
                    g_int.SetPointError(idx_c, 0, err)
                g_int.SetLineColor(ROOT.kBlack)
                g_int.SetMarkerColor(ROOT.kBlack)
                g_int.SetLineStyle(1)
                g_int.SetLineWidth(2)
                g_int.SetMarkerStyle(21)
                g_int.Draw("LP SAME")
                graphs_etaGap.append(g_int)  # keep alive
                leg.AddEntry(g_int, "Integrated", "lp")
            for j in range(nEtaGap):
                key = f"{type_}_{j}"
                if key not in h_diff_etaGap: continue
                h_hist = h_diff_etaGap[key]
                # Convert histogram to TGraphErrors at fixed centroids
                graph = ROOT.TGraphErrors(len(centroids))
                for idx_c, cent in enumerate(centroids):
                    bin_idx = h_hist.FindBin(cent)
                    y = h_hist.GetBinContent(bin_idx)
                    err = h_hist.GetBinError(bin_idx)
                    graph.SetPoint(idx_c, cent, y)
                    graph.SetPointError(idx_c, 0, err)
                graph.SetLineColor(colorsEtaGap[j])
                graph.SetMarkerColor(colorsEtaGap[j])
                graph.SetLineStyle(1)
                graph.SetLineWidth(2)
                graph.SetMarkerStyle(20)
                graph.SetMarkerSize(1.0)
                graph.Draw("LP SAME")
                graphs_etaGap.append(graph)
                leg.AddEntry(graph, etaGapLabels[j], "p")
            leg.Draw()
            pad.Modified()
            pad.Update()

        # Ensure SumPt canvas is fully updated before saving
        c1.cd()
        c1.Modified()
        c1.Update()
        c1.SaveAs("fig_diff_SumPt.pdf")

        # Ensure EtaGap canvas is fully updated before saving
        c2.cd()
        c2.Modified()
        c2.Update()
        c2.SaveAs("fig_diff_EtaGap.pdf")

    elif mode == 'intg':
        # Read ALICE preliminary data
        cent, dDelta, dDeltaErr, dGamma, dGammaErr = read_alice_preliminary(alice_csv)
        # Prepare integrated model histograms
        types = ["delta", "gamma"]
        h_model_int = {}
        for type_ in types:
            p_ss = f.Get(f"p_{type_}_antilambda_antiproton")
            p_ss.Add(f.Get(f"p_{type_}_proton_lambda"))
            h_ss = convert_profile(p_ss)
            p_os = f.Get(f"p_{type_}_antilambda_proton")
            p_os.Add(f.Get(f"p_{type_}_antiproton_lambda"))
            h_os = convert_profile(p_os)
            h_diff = h_os.Clone(f"h_model_int_{type_}")
            h_diff.Add(h_ss, -1.0)
            h_model_int[type_] = h_diff
        # Determine y-axis range across both observables
        y_min, y_max = float('inf'), float('-inf')
        for type_, yerrs in zip(types, [dDeltaErr, dGammaErr]):
            h = h_model_int[type_]
            for cent_val in cent:
                bin_idx = h.FindBin(cent_val)
                val = h.GetBinContent(bin_idx)
                err = h.GetBinError(bin_idx)
                y_min = min(y_min, val - err)
                y_max = max(y_max, val + err)
        # Create single canvas
        c = ROOT.TCanvas("c_intg", "Integrated", 800, 600)
        pad = c.cd()
        pad.SetGrid()
        # Draw frame with fixed axis ranges and labels
        frame = pad.DrawFrame(10.0, -0.0005, 60.0, 0.013, ";Centrality (%);Correlations")
        # Force scientific notation on y-axis
        frame.GetYaxis().SetNoExponent(False)
        frame.GetYaxis().SetMaxDigits(3)
        # frame.GetXaxis().SetNdivisions(5, ROOT.kFALSE)  # Remove or comment out to avoid conflict
        # frame.GetXaxis().SetTickLength(0.03)            # Remove or comment out to avoid conflict
        # Plot model bands and lines with exactly 5 centrality points
        color_map = {"delta": ROOT.kRed+2, "gamma": ROOT.kBlue+2}
        graphs_intg = []
        centroids = [15, 25, 35, 45, 55]
        for type_ in types:
            h = h_model_int[type_]
            g_mod = ROOT.TGraphErrors(len(centroids))
            for i_c, cent_val in enumerate(centroids):
                bin_idx = h.FindBin(cent_val)
                val = h.GetBinContent(bin_idx)
                err = h.GetBinError(bin_idx)
                g_mod.SetPoint(i_c, cent_val, val)
                g_mod.SetPointError(i_c, 0, err)
            col = color_map[type_]
            g_mod.SetFillColorAlpha(col, 0.3)
            g_mod.SetLineColor(col)
            g_mod.SetLineWidth(2)
            g_mod.Draw("3 SAME")
            graphs_intg.append(g_mod)
        # Plot ALICE data points
        data_vals = {"delta": (dDelta, dDeltaErr), "gamma": (dGamma, dGammaErr)}
        for type_ in types:
            yvals, yerrs = data_vals[type_]
            n = len(cent)
            g_data = ROOT.TGraphErrors(n)
            for i_c, cent_val in enumerate(cent):
                g_data.SetPoint(i_c, cent_val, yvals[i_c])
                g_data.SetPointError(i_c, 0, yerrs[i_c])
            col = color_map[type_]
            g_data.SetLineColor(col)
            g_data.SetMarkerColor(col)
            g_data.SetMarkerStyle(20)    # solid circle
            g_data.SetMarkerSize(1.5)    # larger marker
            g_data.SetLineWidth(0)
            g_data.Draw("P SAME")
            graphs_intg.append(g_data)
        # Add legend at top-left for model bands and data points
        legend = ROOT.TLegend(0.2, 0.65, 0.7, 0.85)
        legend.SetNColumns(2)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        # entries: first two are model bands, next two are data points
        legend.AddEntry(graphs_intg[0], "BW + LBC #Delta#delta", "f")
        legend.AddEntry(graphs_intg[2], "ALICE Preliminary #Delta#delta", "p")
        legend.AddEntry(graphs_intg[1], "BW + LBC #Delta#gamma", "f")
        legend.AddEntry(graphs_intg[3], "ALICE Preliminary #Delta#gamma", "p")
        legend.Draw()
        # Save
        c.SaveAs("fig_intg_combined.pdf")

if __name__ == "__main__":
    draw_diff()