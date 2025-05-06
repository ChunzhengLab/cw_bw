import ROOT
import sys

def convert_profile(p):
    nbins = p.GetNbinsX()
    xmin = p.GetXaxis().GetXmin()
    xmax = p.GetXaxis().GetXmax()
    h = ROOT.TH1D(f"{p.GetName()}_clone", "", nbins, xmin, xmax)
    for b in range(1, nbins + 1):
        h.SetBinContent(b, p.GetBinContent(b))
        h.SetBinError(b, p.GetBinError(b))
    return h

def draw_diff(filename="results_combine.root"):
    f = ROOT.TFile.Open(filename)
    if not f or f.IsZombie():
        print(f"Cannot open file {filename}")
        return

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
            pad.DrawFrame(10., -0.002, 60., 0.02, ";centrality (%);#Delta#delta")
        else:
            pad.DrawFrame(10., -0.0010, 60., 0.010, ";centrality (%);#Delta#gamma")
        leg = ROOT.TLegend(0.2,0.6,0.4,0.8)
        legends.append(leg)
        leg.SetHeader("#Sigma p_{T} selection", "C")
        if type_ in h_diff_int:
            h_diff_int[type_].Draw("E2 SAME")
            leg.AddEntry(h_diff_int[type_], "Integrated", "f")
        for i in range(nSumPt):
            key = f"{type_}_{i}"
            if key not in h_diff_sumPt: continue
            h = h_diff_sumPt[key]
            h.SetLineColor(colorsSumPt[i])
            h.SetMarkerColor(colorsSumPt[i])
            h.SetLineStyle(1)
            h.SetLineWidth(2)
            h.SetMarkerStyle(20)
            h.SetMarkerSize(1.0)
            h.Draw("P SAME")
            leg.AddEntry(h, sumPtLabels[i], "p")
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
            pad.DrawFrame(10., -0.002, 60., 0.02, ";centrality (%);#Delta#delta")
        else:
            pad.DrawFrame(10., -0.00035, 60., 0.0035, ";centrality (%);#Delta#gamma")
        leg = ROOT.TLegend(0.2,0.6,0.4,0.8)
        legends.append(leg)
        leg.SetHeader("#Delta#eta selection", "C")
        if type_ in h_diff_int:
            h_diff_int[type_].Draw("E2 SAME")
            leg.AddEntry(h_diff_int[type_], "Integrated", "f")
        for j in range(nEtaGap):
            key = f"{type_}_{j}"
            if key not in h_diff_etaGap: continue
            h = h_diff_etaGap[key]
            h.SetLineColor(colorsEtaGap[j])
            h.SetMarkerColor(colorsEtaGap[j])
            h.SetLineStyle(1)
            h.SetLineWidth(2)
            h.SetMarkerStyle(20)
            h.SetMarkerSize(1.0)
            h.Draw("P SAME")
            leg.AddEntry(h, etaGapLabels[j], "p")
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

if __name__ == "__main__":
    draw_diff()