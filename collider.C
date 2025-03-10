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

using namespace std;

// Setup functions to calculate the parameters
Double_t calcProb(Double_t *x, Double_t *par){
    return (x[0] * x[0] * par[0]) / (1 + exp((x[0] - par[1])/par[2]));
}
Double_t calcD(Double_t *x, Double_t *par){
    return x[0]*2*TMath::Pi();
}

vec_1d values_b()
{
	int divisions = 50;
	double min = 0.0;
	double max = 14;
	double espacamento = (max - min)/divisions;
	vec_1d arr;
	arr.push_back(min);
	
	for(int i = 0; i < divisions; i++)
	{
		arr.push_back(arr[i] + espacamento);
	}
	return arr;
}

// Main code
void collider (int nucleons = 208, int sim = 1e3){
    
    //values of the b parameter
    vec_1d b = values_b();
    int len = static_cast<int>(b.size());
    for(int iterations = 0; iterations < len; iterations++)
    { 
    	    
	    vec_2d xPos;
	    vec_2d yPos;
	    vec_2d coll;
	    auto *random = new TRandom();
	    random->SetSeed();

	    // Simulation variables
	    Double_t  xFirst[nucleons] , yFirst[nucleons],  // convert to 2d matrix
		      xSecond[nucleons], ySecond[nucleons];
	    unordered_set<int> nucleonPartTemp;

	    vec_1d xFirstPartTemp, yFirstPartTemp;
	    vec_1d xSecondPartTemp, ySecondPartTemp;
	    vec_1d num_coll;

	    for (int p = 0; p < sim; p++) 
	    {
	    bool collided = false;
	    while(collided == false){
		  nucleonPartTemp.clear();
		  xFirstPartTemp.clear(), yFirstPartTemp.clear();
		  xSecondPartTemp.clear(), ySecondPartTemp.clear();

		  // Set up the functions to be used in the generator;
		  auto *pos = new TF1("pos", calcProb, 0, 14, 3);
		  pos->SetParameters(p0, r0, a);

		  // Generate collision data
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
			int nColl = 0;
		  // Verify collision partTemp and number
		  for (int i = 0; i < nucleons; i++) {
		  	
		      bool passThrough = false;
		      for (int j = 0; j < nucleons; j++) {
		          if (pow(xFirst[i] - xSecond[j], 2) +
		              pow(yFirst[i] - ySecond[j], 2) < (radiusSq)) {
		              passThrough = true;
		              nColl += 1;
		              if (nucleonPartTemp.insert(j).second) {
		                  xSecondPartTemp.emplace_back(xSecond[j]);
		                  ySecondPartTemp.emplace_back(ySecond[j]);
		                  
		                  
		                  
		              }
		          }
		      }
		      if (passThrough){
		          xFirstPartTemp.emplace_back(xFirst[i]);
			    yFirstPartTemp.emplace_back(yFirst[i]);
			    nColl += 1;
			   num_coll.push_back(nColl);
			
		          
		      }
		  }
		  
		  
		 if(xFirstPartTemp.empty() == false) {collided = true;}
		  	
		  // Convert from list to array for usage in the TGraph module.
		 if(collided)
		 {
		  int fSize = static_cast<int>(xFirstPartTemp.size()),
		      sSize = static_cast<int>(xSecondPartTemp.size());
		  Double_t xfp[fSize], yfp[fSize],
		         xsp[sSize], ysp[sSize];
		  copy(xFirstPartTemp.begin(), xFirstPartTemp.end(), xfp);
		  copy(yFirstPartTemp.begin(), yFirstPartTemp.end(), yfp);
		  copy(xSecondPartTemp.begin(), xSecondPartTemp.end(), xsp);
		  copy(ySecondPartTemp.begin(), ySecondPartTemp.end(), ysp);
		  
		
		  // Create two vectors that will be used to store data at xPos and yPos
		  
		  vec_1d 
		  	in_x,
		  	in_y;
		  in_x.insert(in_x.end(), xfp, xfp + fSize);
		  in_x.insert(in_x.end(), xsp, xsp + sSize);
		  in_y.insert(in_y.end(), yfp, yfp + fSize);
		  in_y.insert(in_y.end(), ysp, ysp + sSize);
		  
		  // Insert the vectors in another vector
		  xPos.push_back(in_x);
		  yPos.push_back(in_y);
		  coll.push_back(num_coll);
		  }
		  
		 }
	    }
	  	// Adding the results to the 3d vector
	    	xPosSim.push_back(xPos);
	    	yPosSim.push_back(yPos);	    
	    	collisions.push_back(coll);	
    }

}
