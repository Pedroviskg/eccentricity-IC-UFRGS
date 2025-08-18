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
    TFile* file = TFile::Open("data.root", "READ");
    TTree* tree = (TTree*)file->Get("Data");

    std::vector<double>* xVec = nullptr;
    std::vector<double>* yVec = nullptr;

    // Estrutura com nCol, nPart, d (igual ao código 1)
    struct Collisions {
        int nCol;
        int nPart;
        double d;
    } data;

    tree->SetBranchAddress("xCol", &xVec);
    tree->SetBranchAddress("yCol", &yVec);
    tree->SetBranchAddress("Collisions", &data);

    const int nEntries = tree->GetEntries();

    // Mapas e vetores globais para acumular resultados
    std::map<double, std::vector<double>> Erp_map;
    std::map<double, std::vector<double>> Epp_map;
    std::map<double, std::vector<int>> Npart_map;

    std::vector<std::vector<int>> v_par;
    std::vector<std::vector<double>> v_srp, v_epp, v_spp;

    // Loop nos eventos
    for (int i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        int N = xVec->size();
        if (N < 2) continue;

        double b = data.d; // parâmetro de impacto

        // Cálculos de momentos
        double sumx = 0, sumy = 0, sumx2 = 0, sumy2 = 0, sumxy = 0;
        for (int j = 0; j < N; j++) {
            sumx  += (*xVec)[j];
            sumy  += (*yVec)[j];
            sumx2 += (*xVec)[j] * (*xVec)[j];
            sumy2 += (*yVec)[j] * (*yVec)[j];
            sumxy += (*xVec)[j] * (*yVec)[j];
        }

        double meanx = sumx / N;
        double meany = sumy / N;
        double varx  = sumx2 / N - meanx * meanx;
        double vary  = sumy2 / N - meany * meany;
        double cov   = sumxy / N - meanx * meany;

        // Excentricidades
        double erp = std::sqrt(std::pow(varx - vary, 2) + 4 * cov * cov) / (varx + vary);
        double epp = (vary - varx) / (vary + varx);
        double srp = TMath::Pi() * std::sqrt(std::abs(varx * vary));
        double spp = TMath::Pi() * std::sqrt(varx * vary - cov*cov);

        // Armazena nos mapas
        Erp_map[b].push_back(erp);
        Epp_map[b].push_back(epp);
        Npart_map[b].push_back(N);

        // Salva para os gráficos 2D
        v_par.push_back({N});
        v_srp.push_back({srp});
        v_epp.push_back({epp});
        v_spp.push_back({spp});
        
    }

    // Médias por valor de b
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

    // --- Gráficos 1D ---
    TCanvas* c1 = new TCanvas("c1", "Excentricidade RP e PP média vs b", 800, 600);
    TGraph* grErp = new TGraph(bVals.size(), &bVals[0], &erpMeans[0]);
    grErp->SetTitle("Excentricidade RP e PP média vs b; b (fm); Média da Excentricidade");
    grErp->SetMarkerStyle(20);
    grErp->SetMarkerColor(kRed);
    grErp->SetLineColor(kRed);

    TGraph* grEpp = new TGraph(bVals.size());
    for (int i = 0; i < bVals.size(); ++i)
        grEpp->SetPoint(i, bVals[i], eppMeans[i]);
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

    TGraph* grEppNp = new TGraph(nPartMeans.size());
    for (int i = 0; i < nPartMeans.size(); ++i)
        grEppNp->SetPoint(i, nPartMeans[i], eppMeans[i]);
    grEppNp->SetMarkerStyle(20);
    grEppNp->SetMarkerColor(kBlue);
    grEppNp->SetLineColor(kBlue);

    grErpNp->Draw("AP");
    grEppNp->Draw("P SAME");
    grErpNp->GetXaxis()->SetRangeUser(0, *std::max_element(nPartMeans.begin(), nPartMeans.end()) + 10);
    grErpNp->GetYaxis()->SetRangeUser(0, 1);
    c2->SaveAs("eccentricity_vs_participants.png");

    // --- Gráficos 2D ---
    TH2F* overlap = new TH2F("overlap", "Overlap area and participants; N_{part}; Overlap area",
                             120, 0, 416, 120, 0, 40);

    for (int s = 0; s < v_par.size(); s++) {
        double jitter = gRandom->Uniform(-0.5, 0.5);
        overlap->Fill(v_par[s][0] + jitter, v_spp[s][0]);
    }
    overlap->SetStats(0);
    TCanvas* ctest = new TCanvas("ctest", "overlap", 800, 600);
    gStyle->SetPalette(kViridis);
    gPad->SetLogz();
    overlap->Draw("COLZ");
    ctest->SaveAs("overlap.png");

    TH2F* epsilon_heat_map = new TH2F("epsilon_heat_map", "Eccentricity vs N_{part}; N_{part}; Eccentricity",
                                      200, 0, 416, 200, 0, 1);

    for (int s = 0; s < v_par.size(); s++) {
        double jitter = gRandom->Uniform(-0.5, 0.5);
        epsilon_heat_map->Fill(v_par[s][0] + jitter, v_epp[s][0]);
    }
    epsilon_heat_map->SetStats(0);
    TCanvas* ctest2 = new TCanvas("ctest2", "epsilon_heat_map", 800, 600);
    gStyle->SetPalette(kRainBow);
    gPad->SetLogz();
    epsilon_heat_map->Draw("COLZ");
    ctest2->SaveAs("epsilonheatmap.png");
}

