#include <TFile.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <map>
#include <string>

void Draw_diff(const char* filename = "results_combine.root") {
    // Open the ROOT file
    TFile* f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        return;
    }

    // Types and pair definitions
    const char* types[] = {"delta", "gamma"};
    const char* pairsSS[2] = {"antilambda_antiproton", "proton_lambda"};
    const char* pairsOS[2] = {"antiproton_lambda", "antilambda_proton"};

    // Bin counts
    const int nSumPt = 3;
    const int nEtaGap = 4; // only 0~3 available in data

    // Custom color definitions
    int ci[4];
    TColor* tc[4];
    ci[0] = TColor::GetFreeColorIndex();
    tc[0] = new TColor(ci[0],   0/255.,  24/255., 113/255.);//dark blue
    ci[1] = TColor::GetFreeColorIndex();
    tc[1] = new TColor(ci[1], 255/255.,  88/255.,  93/255.);//red
    ci[2] = TColor::GetFreeColorIndex();
    tc[2] = new TColor(ci[2], 255/255., 181/255.,  73/255.);//yellow
    ci[3] = TColor::GetFreeColorIndex();
    tc[3] = new TColor(ci[3], 65/255.,  182/255., 230/255.);//light blue
    int colorsSumPt[nSumPt] = {ci[0], ci[1], ci[2]};
    int colorsEtaGap[nEtaGap] = {ci[0], ci[1], ci[2], ci[3]};

    // Bin labels for legend
    const char* sumPtLabels[nSumPt] = {"1<#Sigma p_{T}<3 GeV/c", "3<#Sigma p_{T}<5 GeV/c", "5<#Sigma p_{T}<8 GeV/c"};
    const char* etaGapLabels[nEtaGap] = {"0.2<#Delta#eta<0.4", "0.4<#Delta#eta<0.6", "0.6<#Delta#eta<0.8", "#Delta#eta>0.8"};

    // Containers for the final difference histograms
    std::map<std::string, TH1D*> h_diff_sumPt;
    std::map<std::string, TH1D*> h_diff_etaGap;

    // Helper to convert TProfile to TH1D
    auto convertProfile = [](TProfile* p) {
        int nbins = p->GetNbinsX();
        double xmin = p->GetXaxis()->GetXmin();
        double xmax = p->GetXaxis()->GetXmax();
        TH1D* h = new TH1D(Form("%s_clone", p->GetName()), "", nbins, xmin, xmax);
        for (int b = 1; b <= nbins; ++b) {
            h->SetBinContent(b, p->GetBinContent(b));
            h->SetBinError(b, p->GetBinError(b));
        }
        return h;
    };

    // Loop over delta/gamma types
    for (const char* type : types) {
        // === SumPt slices ===
        for (int i = 0; i < nSumPt; ++i) {
            TProfile* p_ss1 = (TProfile*)f->Get(Form("p_%s_%s_sumPt_%d", type, pairsSS[0], i));
            TProfile* p_ss2 = (TProfile*)f->Get(Form("p_%s_%s_sumPt_%d", type, pairsSS[1], i));
            if (!p_ss1 || !p_ss2) continue;
            p_ss1->Add(p_ss2);
            TH1D* h_ss = convertProfile(p_ss1);

            TProfile* p_os1 = (TProfile*)f->Get(Form("p_%s_%s_sumPt_%d", type, pairsOS[0], i));
            TProfile* p_os2 = (TProfile*)f->Get(Form("p_%s_%s_sumPt_%d", type, pairsOS[1], i));
            if (!p_os1 || !p_os2) continue;
            p_os1->Add(p_os2);
            TH1D* h_os = convertProfile(p_os1);

            TH1D* h_diff = (TH1D*)h_os->Clone(Form("h_%s_sumPt_%d", type, i));
            h_diff->Add(h_ss, -1.0);
            h_diff_sumPt[Form("%s_%d", type, i)] = h_diff;
        }

        // === EtaGap slices ===
        for (int j = 0; j < nEtaGap; ++j) {
            TProfile* p_ss1 = (TProfile*)f->Get(Form("p_%s_%s_etaGap_%d", type, pairsSS[0], j));
            TProfile* p_ss2 = (TProfile*)f->Get(Form("p_%s_%s_etaGap_%d", type, pairsSS[1], j));
            if (!p_ss1 || !p_ss2) continue;
            p_ss1->Add(p_ss2);
            TH1D* h_ss = convertProfile(p_ss1);

            TProfile* p_os1 = (TProfile*)f->Get(Form("p_%s_%s_etaGap_%d", type, pairsOS[0], j));
            TProfile* p_os2 = (TProfile*)f->Get(Form("p_%s_%s_etaGap_%d", type, pairsOS[1], j));
            if (!p_os1 || !p_os2) continue;
            p_os1->Add(p_os2);
            TH1D* h_os = convertProfile(p_os1);

            TH1D* h_diff = (TH1D*)h_os->Clone(Form("h_%s_etaGap_%d", type, j));
            h_diff->Add(h_ss, -1.0);
            h_diff_etaGap[Form("%s_%d", type, j)] = h_diff;
        }
    }

    // Set global style
    gStyle->SetOptStat(0);

    // === Compute integrated difference ===
    std::map<std::string, TH1D*> h_diff_int;
    for (const char* type : types) {
        // SS integrated
        TProfile* p_ssi1 = (TProfile*)f->Get(Form("p_%s_%s", type, pairsSS[0]));
        TProfile* p_ssi2 = (TProfile*)f->Get(Form("p_%s_%s", type, pairsSS[1]));
        if (!p_ssi1 || !p_ssi2) continue;
        p_ssi1->Add(p_ssi2);
        TH1D* h_ssi = convertProfile(p_ssi1);
        // OS integrated
        TProfile* p_osi1 = (TProfile*)f->Get(Form("p_%s_%s", type, pairsOS[0]));
        TProfile* p_osi2 = (TProfile*)f->Get(Form("p_%s_%s", type, pairsOS[1]));
        if (!p_osi1 || !p_osi2) continue;
        p_osi1->Add(p_osi2);
        TH1D* h_osi = convertProfile(p_osi1);
        // diff
        TH1D* h_int = (TH1D*)h_osi->Clone(Form("h_int_%s", type));
        h_int->Add(h_ssi, -1.0);
        // style for fill
        h_int->SetFillColorAlpha(kRed+2, 0.3);
        h_int->SetLineColor(kRed+2);
        h_int->SetMarkerColor(kRed+2);
        h_diff_int[Form("%s", type)] = h_int;
    }

    // === Plot SumPt dependence ===
    TCanvas* c1 = new TCanvas("c1", "SumPt dependence", 1200, 500);
    c1->Divide(2,1);

    // Left pad: delta
    c1->cd(1);
    gPad->SetGrid();
    gPad->DrawFrame(10., -0.002, 60., 0.02, ";centrality (%);#Delta#delta");
    TLegend* leg1 = new TLegend(0.2,0.6,0.4,0.8);
    // draw integrated band
    if (h_diff_int.count("delta")) {
        h_diff_int["delta"]->Draw("E2 SAME");
        leg1->AddEntry(h_diff_int["delta"], "Integrated", "f");
    }
    for (int i = 0; i < nSumPt; ++i) {
        TH1D* h = h_diff_sumPt[Form("delta_%d", i)];
        if (!h) continue;
        h->SetLineColor(colorsSumPt[i]);
        h->SetLineStyle(1);
        h->SetLineWidth(2);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1.0);
        h->SetMarkerColor(colorsSumPt[i]);
        h->Draw("LP SAME");
        leg1->AddEntry(h, sumPtLabels[i], "p");
    }
    leg1->Draw();

    // Right pad: gamma
    c1->cd(2);
    gPad->SetGrid();
    gPad->DrawFrame(10., -0.0010, 60., 0.010, ";centrality (%);#Delta#gamma");
    TLegend* leg2 = new TLegend(0.2,0.6,0.4,0.8);
    // draw integrated band
    if (h_diff_int.count("gamma")) {
        h_diff_int["gamma"]->Draw("E2 SAME");
        leg2->AddEntry(h_diff_int["gamma"], "Integrated", "f");
    }
    for (int i = 0; i < nSumPt; ++i) {
        TH1D* h = h_diff_sumPt[Form("gamma_%d", i)];
        if (!h) continue;
        h->SetLineColor(colorsSumPt[i]);
        h->SetLineStyle(1);
        h->SetLineWidth(2);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1.0);
        h->SetMarkerColor(colorsSumPt[i]);
        h->Draw("LP SAME");
        leg2->AddEntry(h, sumPtLabels[i], "p");
    }
    leg2->Draw();

    // === Plot EtaGap dependence ===
    TCanvas* c2 = new TCanvas("c2", "EtaGap dependence", 1200, 500);
    c2->Divide(2,1);

    // Left pad: delta
    c2->cd(1);
    gPad->SetGrid();
    gPad->DrawFrame(10., -0.002, 60., 0.02, ";centrality (%);#Delta#delta");
    TLegend* leg3 = new TLegend(0.2,0.6,0.4,0.8);
    // draw integrated band
    if (h_diff_int.count("delta")) {
        h_diff_int["delta"]->Draw("E2 SAME");
        leg3->AddEntry(h_diff_int["delta"], "Integrated", "f");
    }
    for (int j = 0; j < nEtaGap; ++j) {
        TH1D* h = h_diff_etaGap[Form("delta_%d", j)];
        if (!h) continue;
        h->SetLineColor(colorsEtaGap[j]);
        h->SetLineStyle(1);
        h->SetLineWidth(2);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1.0);
        h->SetMarkerColor(colorsEtaGap[j]);
        h->Draw("LP SAME");
        leg3->AddEntry(h, etaGapLabels[j], "p");
    }
    leg3->Draw();

    // Right pad: gamma
    c2->cd(2);
    gPad->SetGrid();
    gPad->DrawFrame(10., -0.00035, 60., 0.0035, ";centrality (%);#Delta#gamma");
    TLegend* leg4 = new TLegend(0.2,0.6,0.4,0.8);
    // draw integrated band
    if (h_diff_int.count("gamma")) {
        h_diff_int["gamma"]->Draw("E2 SAME");
        leg4->AddEntry(h_diff_int["gamma"], "Integrated", "f");
    }
    for (int j = 0; j < nEtaGap; ++j) {
        TH1D* h = h_diff_etaGap[Form("gamma_%d", j)];
        if (!h) continue;
        h->SetLineColor(colorsEtaGap[j]);
        h->SetLineStyle(1);
        h->SetLineWidth(2);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1.0);
        h->SetMarkerColor(colorsEtaGap[j]);
        h->Draw("LP SAME");
        leg4->AddEntry(h, etaGapLabels[j], "p");
    }
    leg4->Draw();

    // // Cleanup
    // f->Close();
}
