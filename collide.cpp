#include "TF1.h"
#include "TRandom.h"
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <thread>
#include <mutex>

using namespace std;
mutex mtx;

// Structures used to generate and save data.
struct Nucleon {Double_t x, y, z;};

struct Data {
    int nCol;
    int nPart;
    double d;
};

struct Position{
    std::vector<double> xCol;
    std::vector<double> yCol;
};

// Const values used in the simulations
const Double_t
       pi       = TMath::Pi(),
       r0       = 6.62,       // fm
       a        = 0.542,      // fm
       sigma    = 6.5,        // fm^2
       radiusSq = sigma / (10 * pi), // fm^2
       minDisSq = 0;          // fm^2

// Main code
void collide (const string& filename = "./data.root",
              const int nucleons      = 208,
              const int simulations   = 1e6,
              const ::int64_t nThreads =std::thread::hardware_concurrency()){

    // Configure the program to export in a Root file
    Data data{};
    Position position{};
    auto *file = new TFile(filename.c_str(),"recreate");
    auto *save = new TTree("Data", "Simulation");
    save->Branch("Collisions", &data, "NColl/I:NPart/I:Dist/D");
    save->Branch("xCol", &position.xCol);
    save->Branch("yCol", &position.yCol);

    auto *distF = new TF1("distF", "x*2*pi", 0, 2.8*r0);
    auto *posF = new TF1("posF",  "(x*x)/(1+exp((x-[0])/[1]))", 0, 2*r0);
    posF->SetParameters(r0, a);

    // Configure multithreading
    ROOT::EnableImplicitMT();
    ROOT::EnableThreadSafety();
    vector<thread> threads;
    threads.reserve(nThreads);
    for (int k = 0; k < nThreads; k++) {
        threads.emplace_back([&, k] {

            // Run the simulations
            for (int p = k; p < simulations;){
                unordered_set<int> nucleonPartTemp;

                // Create the array of Nucleon to represent the Nuclei.
                Nucleon n1[nucleons];
                Nucleon n2[nucleons];

                Double_t cosTheta, sinTheta,
                         phi,  pos;

                int nCol = 0;
                int nPart = 0;
                std::vector<double> xCol;
                std::vector<double> yCol;
                Double_t d;

                // Generate collision data
                d = distF->GetRandom();
                for (int i = 0; i < nucleons; i++) {

                    // Generate rdm cylindrical coordinates for the n1 nucleus
                    cosTheta = 2 * gRandom->Rndm() - 1;
                    sinTheta = sqrt(1 - cosTheta * cosTheta);
                    phi = 2 * pi * gRandom->Rndm();
                    pos = posF->GetRandom();
                    n1[i].x = pos * sinTheta * cos(phi);
                    n1[i].y = pos * sinTheta * sin(phi);
                    n1[i].z = pos * cosTheta;

                    // Do the same for the n2 nucleus
                    cosTheta = 2 * gRandom->Rndm() - 1;
                    sinTheta = sqrt(1 - cosTheta * cosTheta);
                    phi = 2 * pi * gRandom->Rndm();
                    pos = posF->GetRandom();
                    n2[i].x = pos * sinTheta * cos(phi) + d;
                    n2[i].y = pos * sinTheta * sin(phi);
                    n2[i].z = pos * cosTheta;

                    // Check if no nucleon breaks the pauli exclusion principle
                    for (int j = 0; j < i; j++){
                        Double_t dx = n1[i].x - n1[j].x;
                        Double_t dy = n1[i].y - n1[j].y;
                        Double_t dz = n1[i].z - n1[j].z;
                        if ((dx*dx + dy*dy + dz*dz) < minDisSq){
                            cosTheta = 2 * gRandom->Rndm() - 1;
                            sinTheta = sqrt(1 - cosTheta * cosTheta);
                            phi = 2 * pi * gRandom->Rndm();
                            pos = posF->GetRandom();
                            n1[i].x = pos * sinTheta * cos(phi);
                            n1[i].y = pos * sinTheta * sin(phi);
                            n1[i].z = pos * cosTheta;
                            j = 0;
                        }
                        dx = n2[i].x - n2[j].x;
                        dy = n2[i].y - n2[j].y;
                        dz = n2[i].z - n2[j].z;
                        if ((dx*dx + dy*dy + dz*dz) < minDisSq){
                            cosTheta = 2 * gRandom->Rndm() - 1;
                            sinTheta = sqrt(1 - cosTheta * cosTheta);
                            phi = 2 * pi * gRandom->Rndm();
                            pos = posF->GetRandom();
                            n2[i].x = pos * sinTheta * cos(phi) + d;
                            n2[i].y = pos * sinTheta * sin(phi);
                            n2[i].z = pos * cosTheta;
                            j = 0;
                        }
                    }
                }

                // Verify collision partTemp and number
                Double_t dx, dy;
                bool part1[nucleons];
                bool part2[nucleons];
                for (int i = 0; i < nucleons; i++) {part1[i] = false; part2[i] = false;}
                for (int i = 0; i < nucleons; i++) {
                    for (int j = 0; j < nucleons; j++) {

                        dx = n1[i].x - n2[j].x;
                        dy = n1[i].y - n2[j].y;

                        if ((dx * dx + dy * dy) < radiusSq) {
                            nCol++;
                            part1[i] = true;
                            part2[j] = true;
                            xCol.push_back(n1[i].x);
                            xCol.push_back(n2[j].x);
                            yCol.push_back(n1[i].y);
                            yCol.push_back(n2[j].y);
                        }
                    }
                }
                for (int i = 0; i < nucleons; ++i) {
                    nPart += part1[i] + part2[i];
                }

                // Write the data to the root file
                p+= (int) nThreads;
                if (nPart > 0) {
                    mtx.lock();
                    data.nCol = nCol;
                    data.nPart = nPart;
                    data.d = d;
                    position.xCol = xCol;
                    position.yCol = yCol;
                    save->Fill();
                    mtx.unlock();
                }
            }
        });
    }
    for (auto &t : threads) t.join();
    save->Write();
    file->Close();
}
