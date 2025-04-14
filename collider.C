#include <TF1.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

// Constants
const Double_t 
	pi       	= TMath::Pi(),
 	p0       	= 3,
 	r0       	= 6.62,      // fm
 	a        	= 0.542,    // fm
 	sigma    	= 6.5,      // fm^2
 	radiusSq 	= sigma / pi, // fm^2
 	minDisSq 	= 0,          // fm^2
 	nucleons        = 208;
	min             = 0
	max             = 14

Double_t calcProb(Double_t *x, Double_t *par){
    return (x[0] * x[0] * par[0]) / (1 + exp((x[0] - par[1])/par[2]));
}

std::vector<double> values_b()
{
    int divisions = 60;
    double min_b = min;
    double max_b = max + 1;
    double space = (max - min) / divisions;

    std::vector<double> arr;

    for (int i = 0; i <= divisions; i++) {
        arr.push_back(min + i * space);
    }

    return arr;
}
 

void collider(int nucleons = 208, int sim = 1e4) {



    // Creating .root file
    TFile *outputFile = new TFile("collisions.root", "RECREATE");
    TTree *tree = new TTree("collisionTree", "CollideTree");
    std::vector<double> xVec, yVec;
    Double_t impactParameter = 0;
    tree->Branch("x", &xVec);
    tree->Branch("y", &yVec);
    tree->Branch("b", &impactParameter);

    std::vector<double> b = values_b();
    int len = static_cast<int>(b.size());

    auto *pos = new TF1("pos", calcProb, 0, 15, 3);
    pos->SetParameters(p0, r0, a);
    auto *random = new TRandom();
    random->SetSeed();
    for(int iterations = 0; iterations < len; iterations++) { 

    	double d = b[iterations]
    	int localSim = sim;
    	
    	Double_t  xFirst[nucleons] , yFirst[nucleons],
              xSecond[nucleons], ySecond[nucleons];

    	unordered_set<int> nucleonPartTemp;

        std::vector<double> xFirstPartTemp, yFirstPartTemp;
        std::vector<double> xSecondPartTemp, ySecondPartTemp;

        for (int p = 0; p < localSim; p++) 
        {
            bool collided = false;
            while(!collided){
                nucleonPartTemp.clear();
                xFirstPartTemp.clear(); yFirstPartTemp.clear();
                xSecondPartTemp.clear(); ySecondPartTemp.clear();

                

                Double_t d = b[iterations];

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

                for (int i = 0; i < nucleons; i++) {
                    bool passThrough = false;
                    for (int j = 0; j < nucleons; j++) {
                        if (pow(xFirst[i] - xSecond[j], 2) +
                            pow(yFirst[i] - ySecond[j], 2) < radiusSq) {
                            passThrough = true;
                            if (nucleonPartTemp.insert(j).second) {
                                xSecondPartTemp.emplace_back(xSecond[j]);
                                ySecondPartTemp.emplace_back(ySecond[j]);
                            }
                        }
                    }
                    if (passThrough) {
                        xFirstPartTemp.emplace_back(xFirst[i]);
                        yFirstPartTemp.emplace_back(yFirst[i]);
                    }
                }

                if (!xFirstPartTemp.empty()) {
                    collided = true;

                    int fSize = static_cast<int>(xFirstPartTemp.size());
                    int sSize = static_cast<int>(xSecondPartTemp.size());

                    xVec.clear();
                    yVec.clear();

                    xVec.insert(xVec.end(), xFirstPartTemp.begin(), xFirstPartTemp.end());
                    xVec.insert(xVec.end(), xSecondPartTemp.begin(), xSecondPartTemp.end());

                    yVec.insert(yVec.end(), yFirstPartTemp.begin(), yFirstPartTemp.end());
                    yVec.insert(yVec.end(), ySecondPartTemp.begin(), ySecondPartTemp.end());

                    impactParameter = d;
                    tree->Fill();
                }
            }
        }
    }

    // Close the file
    outputFile->cd();
    tree->Write();
    outputFile->Close();
}

