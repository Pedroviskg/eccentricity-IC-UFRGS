#include <TF1.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include "globals.h"
#include "eccentricity.h"

using namespace std;

// Setup functions to calculate the parameters
Double_t calcProb(Double_t *x, Double_t *par){
    return (x[0] * x[0] * par[0]) / (1 + exp((x[0] - par[1])/par[2]));
}
Double_t calcD(Double_t *x, Double_t *par){
    return x[0]*2*TMath::Pi();
}

// Function to pretty print the graphs
void config(TGraph graph, Color_t color, const char* name){
    graph.SetMarkerStyle(kFullCircle);
    graph.SetMarkerColor(color);
    graph.SetMarkerSize(6);
    graph.SetName(name);
    graph.DrawClone("P");
}

// Main code
void collisionsDraw (int nucleons = 208, int sim = 20, const string& location = "./sim"){
    gSystem->mkdir(location.c_str());

    // Config RNG
    auto *random = new TRandom();
    random->SetSeed();

    // Simulation variables
    Double_t  xFirst[nucleons] , yFirst[nucleons],  // convert to 2d matrix
            xSecond[nucleons], ySecond[nucleons];
    int nCol = 0;

    // Config the canvas
    TCanvas canvas("B(Z)", "B(Z) linear", 1280, 1000);
    canvas.SetTicks();
    unordered_set<int> nucleonPartTemp;

    vector<Double_t> xFirstPartTemp, yFirstPartTemp;
    vector<Double_t> xSecondPartTemp, ySecondPartTemp;

    for (int p = 0; p < sim; p++) {
		
        nucleonPartTemp.clear();
        xFirstPartTemp.clear(), yFirstPartTemp.clear();
        xSecondPartTemp.clear(), ySecondPartTemp.clear();


        // Config the graph
        canvas.Clear();
        TH1F *f = canvas.DrawFrame(-12, -15, 25, 15);
        f->GetYaxis()->SetTitle("Y (fm)");
        f->GetXaxis()->SetTitle("X (fm)");

        // Set up the functions to be used in the generators
        auto *dist = new TF1("d", calcD, 0, 14, 0);      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        auto *pos = new TF1("pos", calcProb, 0, 14, 3);
        pos->SetParameters(p0, r0, a);

        // Generate collision data
        Double_t d = dist->GetRandom(random);
        for (int i = 0; i < nucleons; i++) {

            Double_t position = pos->GetRandom(random);
            Double_t phi = random->Rndm() * 2 * pi;
            Double_t cTheta = 2 * gRandom->Rndm() - 1 ;
            Double_t sTheta = TMath::Sqrt(1 - cTheta * cTheta);

            xFirst[i] = position * sin(phi) * sTheta;
            yFirst[i] = position * cos(phi) * sTheta;
            position = pos->GetRandom(random);
            phi = random->Rndm() * 2 * pi;
            cTheta = 2 * gRandom->Rndm() - 1 ;
            sTheta = TMath::Sqrt(1 - cTheta * cTheta);
            xSecond[i] = position * sin(phi) * sTheta + d;
            ySecond[i] = position * cos(phi) * sTheta;
        }

        // Verify collision partTemp and number
        for (int i = 0; i < nucleons; i++) {
            bool passThrough = false;
            for (int j = 0; j < nucleons; j++) {
                if (pow(xFirst[i] - xSecond[j], 2) +
                    pow(yFirst[i] - ySecond[j], 2) < (radiusSq)) {
                    nCol++;
                    passThrough = true;
                    if (nucleonPartTemp.insert(j).second) {
                        xSecondPartTemp.emplace_back(xSecond[j]);
                        ySecondPartTemp.emplace_back(ySecond[j]);
                        
                        
                    }
                }
            }
            if (passThrough){
                xFirstPartTemp.emplace_back(xFirst[i]);
		    yFirstPartTemp.emplace_back(yFirst[i]);
             
                
                
                
            }
        }
        
        
        // Convert from list to array for usage in the TGraph module.
        int fSize = static_cast<int>(xFirstPartTemp.size()),
            sSize = static_cast<int>(xSecondPartTemp.size());
        Double_t xfp[fSize], yfp[fSize],
               xsp[sSize], ysp[sSize];
        copy(xFirstPartTemp.begin(), xFirstPartTemp.end(), xfp);
        copy(yFirstPartTemp.begin(), yFirstPartTemp.end(), yfp);
        copy(xSecondPartTemp.begin(), xSecondPartTemp.end(), xsp);
        copy(ySecondPartTemp.begin(), ySecondPartTemp.end(), ysp);
        
        
        // Drawn the nuclei in canvas
        TGraph first(nucleons, xFirst, yFirst);
        TGraph second(nucleons, xSecond, ySecond);
        TGraph firstPart(fSize, xfp, yfp);
        TGraph secondPart(sSize, xsp, ysp);
        config(first, kAzure + 1, "nas");
        config(second, kOrange + 7, "nbs");
        config(firstPart, kAzure + 2, "nap");
        config(secondPart, kRed, "nbp");
	
	  // Create two vectors that will be used to store data at xPos and yPos
	  vector<Double_t> 
		in_x,
		in_y;
	  in_x.insert(in_x.end(), xfp, xfp + fSize);
	  in_x.insert(in_x.end(), xsp, xsp + sSize);
	  in_y.insert(in_y.end(), yfp, yfp + fSize);
	  in_y.insert(in_y.end(), ysp, ysp + sSize);

	  // Insert the vectors in another vector
          xPos.push_back(in_x);
          yPos.push_back(in_y);

        // Create legend
        auto leg = new TLegend(.1, .9, .4, .7);
        leg->SetHeader("Legenda", "c");
        leg->AddEntry("nas", "Espectadores - Nucleo A", "p");
        leg->AddEntry("nbs", "Espectadores - Nucleo B", "p");
        leg->AddEntry("nap", "Participantes - Nucleo A", "p");
        leg->AddEntry("nbp", "Participantes - Nucleo B", "p");
        leg->DrawClone("Same");

        // Export as image the canvas
        string tName = location + "/Sim" + to_string(p) + ".png";
        canvas.SaveAs(tName.c_str());
        
        
    }

 //std::cout << xPos[1];

}
