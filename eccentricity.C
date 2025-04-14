#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMath.h>
#include <TH2F.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

void eccentricity() {
    TFile* file = TFile::Open("collisions.root", "READ");
    TTree* tree = (TTree*)file->Get("collisionTree");

    std::vector<double>* xVec = nullptr;
    std::vector<double>* yVec = nullptr;
    double b = 0.0;

    tree->SetBranchAddress("x", &xVec);
    tree->SetBranchAddress("y", &yVec);
    tree->SetBranchAddress("b", &b);

    std::map<double, std::vector<double>> Erp_map;
    std::map<double, std::vector<double>> Epp_map;
    std::map<double, std::vector<int>> Npart_map;

    std::vector<std::vector<int>> v_par;
    std::vector<std::vector<double>> v_srp, v_epp;

    const int nEntries = tree->GetEntries();

    for (int i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        int N = xVec->size();
        if (N < 2) continue;

        double sumx = 0, sumy = 0, sumx2 = 0, sumy2 = 0, sumxy = 0;
        for (int j = 0; j < N; j++) {
            sumx += (*xVec)[j];
            sumy += (*yVec)[j];
            sumx2 += (*xVec)[j] * (*xVec)[j];
            sumy2 += (*yVec)[j] * (*yVec)[j];
            sumxy += (*xVec)[j] * (*yVec)[j];
        }

        double meanx = sumx / N;
        double meany = sumy / N;
        double varx = sumx2 / N - meanx * meanx;
        double vary = sumy2 / N - meany * meany;
        double cov  = sumxy / N - meanx * meany;

        double erp = std::sqrt(std::pow(varx - vary, 2) + 4 * cov * cov) / (varx + vary);
        double epp = (vary - varx) / (vary + varx);
        double srp = TMath::Pi() * std::sqrt(std::abs(varx * vary));

        Erp_map[b].push_back(erp);
        Epp_map[b].push_back(epp);
        Npart_map[b].push_back(N);

        v_par.push_back({N});
        v_srp.push_back({srp});
        v_epp.push_back({epp});
    }

    std::vector<double> bVals, erpMeans, eppMeans, nPartMeans;
    for (const auto& pair : Erp_map) {
        double bval = pair.first;

        const std::vector<double>& erplist = pair.second;
        const std::vector<double>& epplist = Epp_map[bval];
        const std::vector<int>& nplist = Npart_map[bval];

        double media_rp = TMath::Mean(erplist.begin(), erplist.end());
        double media_pp = TMath::Mean(epplist.begin(), epplist.end());
        double media_np = TMath::Mean(nplist.begin(), nplist.end());

        bVals.push_back(bval);
        erpMeans.push_back(media_rp);
        eppMeans.push_back(media_pp);
        nPartMeans.push_back(media_np);
    }

    TCanvas* c1 = new TCanvas("c1", "Excentricidade RP e PP média vs b", 800, 600);
    TGraph* grErp = new TGraph(bVals.size(), &bVals[0], &erpMeans[0]);
    grErp->SetTitle("Excentricidade RP e PP média vs b; b (fm); Média da Excentricidade");
    grErp->SetMarkerStyle(20);
    grErp->SetMarkerColor(kRed);
    grErp->SetLineColor(kRed);

    TGraph* grEpp = new TGraph(bVals.size() + 1);
    grEpp->SetPoint(0, 0, 0);
    for (int i = 0; i < bVals.size(); ++i)
        grEpp->SetPoint(i + 1, bVals[i], eppMeans[i]);
    grEpp->SetMarkerStyle(20);
    grEpp->SetMarkerColor(kBlue);
    grEpp->SetLineColor(kBlue);

    grErp->Draw("AP");
    grEpp->Draw("P SAME");
    grErp->GetXaxis()->SetRangeUser(0, 14);
    grErp->GetYaxis()->SetRangeUser(0, 1);
    c1->SaveAs("eccentricity_vs_b.png");

    TCanvas* c2 = new TCanvas("c2", "Excentricidade vs Participantes", 800, 600);
    TGraph* grErpNp = new TGraph(nPartMeans.size(), &nPartMeans[0], &erpMeans[0]);
    grErpNp->SetTitle("Excentricidade RP e PP média vs N_{part}; N_{part}; Média da Excentricidade");
    grErpNp->SetMarkerStyle(20);
    grErpNp->SetMarkerColor(kRed);
    grErpNp->SetLineColor(kRed);

    TGraph* grEppNp = new TGraph(nPartMeans.size() + 1);
    grEppNp->SetPoint(0, 0, 0);
    for (int i = 0; i < nPartMeans.size(); ++i)
        grEppNp->SetPoint(i + 1, nPartMeans[i], eppMeans[i]);
    grEppNp->SetMarkerStyle(20);
    grEppNp->SetMarkerColor(kBlue);
    grEppNp->SetLineColor(kBlue);

    grErpNp->Draw("AP");
    grEppNp->Draw("P SAME");
    grErpNp->GetXaxis()->SetRangeUser(0, *std::max_element(nPartMeans.begin(), nPartMeans.end()) + 10);
    grErpNp->GetYaxis()->SetRangeUser(0, 1);
    c2->SaveAs("eccentricity_vs_participants.png");

    TH2F* overlap = new TH2F("overlap", "Overlap area and participants histogram; N_{part}; A_{overlap} (fm^{2})", 200, 0, 416, 200, 0, 40);

    for (int s = 0; s < v_par.size(); s++) {
        for (int m = 0; m < v_par[s].size(); m++) {
            double jitter = gRandom->Uniform(-0.5, 0.5);
            overlap->Fill(v_par[s][m] + jitter, v_srp[s][m]);
        }
    }

    TCanvas* ctest = new TCanvas("ctest", "overlap", 800, 600);
    gStyle->SetPalette(kViridis);
    overlap->SetMinimum(0); 
    overlap->SetMaximum(overlap->GetMaximum());
    gPad->Update();
    gPad->SetLogz();
    overlap->Draw("COLZ");
    ctest->SaveAs("overlap.png");

    TH2F* epsilon_heat_map = new TH2F("epsilon_heat_map", "Eccentricity and participants histogram; N_{part}; #varepsilon_{PP}", 200, 0, 416, 200, 0, 1);

    for (int s = 0; s < v_par.size(); s++) {
        for (int m = 0; m < v_par[s].size(); m++) {
            double jitter = gRandom->Uniform(-0.5, 0.5);
            epsilon_heat_map->Fill(v_par[s][m] + jitter, v_epp[s][m]);
        }
    }

    TCanvas* ctest2 = new TCanvas("ctest2", "epsilon_heat_map", 800, 600);
    gStyle->SetPalette(kRainBow);
    epsilon_heat_map->SetMinimum(0);
    epsilon_heat_map->SetMaximum(epsilon_heat_map->GetMaximum());
    gPad->Update();
    gPad->SetLogz();
    epsilon_heat_map->Draw("COLZ");
    ctest2->SaveAs("epsilonheatmap.png");
}

