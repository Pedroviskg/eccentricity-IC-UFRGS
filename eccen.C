#include "eccentricity.h"
#include "globals.h"

// Calculate standard deviation
std::vector<Double_t> get_sigma(std::vector<vector<Double_t>> pos)
{
	std::vector<Double_t> sigma;

	for(int i = 0; i < pos.size(); i++)
	{
		Double_t 
			sum_mean   = 0,
			sum_sqmean = 0;
			
		int n = pos[i].size();
		for(int j = 0; j < n; j++)
		{
			sum_mean   = sum_mean + pos[i][j];
			sum_sqmean = sum_sqmean + pos[i][j]*pos[i][j];
		}
		
		Double_t mean    = sum_mean/n;
		Double_t mean_sq = mean * mean;
		Double_t sqmean  = sum_sqmean/n;
		Double_t val_sigma = (-mean_sq + sqmean);
		sigma.push_back(val_sigma);
		
	}
	return sigma;
}

// Calculate covarience of x and y
std::vector<Double_t> get_covarience(std::vector<vector<Double_t>> xpos, std::vector<vector<Double_t>> ypos)
{
	std::vector<Double_t> covarience;
	for(int i = 0; i < xpos.size(); i++)
	{
		Double_t 
			sum1 	  = 0,
			sum2 	  = 0,
			prodsum = 0;
		
		// Running the loop to get xmean * ymean
		int n = xpos[i].size();
		
		for(int j = 0; j < n; j++)
		{
			sum1    += xpos[i][j];
			sum2    += ypos[i][j];
			prodsum += xpos[i][j] * ypos[i][j];
		}
	
		Double_t xymean = prodsum/n;
		Double_t x_mean = sum1/n;
		Double_t y_mean = sum2/n;
		Double_t cov    = xymean - x_mean * y_mean;
		covarience.push_back(cov);
	}
	
	return covarience;
}
			
		
void eccen()
{	
	
	int debugging;
	std::cout << "Debugging? (0, 1)"; std::cin >> debugging;
	
	std::vector<Double_t> sqsigmax = get_sigma(xPos);
	std::vector<Double_t> sqsigmay = get_sigma(yPos);
	std::vector<Double_t> cov      = get_covarience(xPos, yPos);
	
	
	int s_size = sqsigmax.size();
	
	
	// Calculating the reaction plane eccentricity
	std::vector<Double_t> e_rp;
	
	
	for(int i = 0; i < s_size; i++)
	{
		Double_t valueserp = abs((-sqsigmay[i] + sqsigmax[i]))/(sqsigmay[i] + sqsigmax[i]);
		e_rp.push_back(valueserp);
	}
	
	
	// Calculating the participant plane eccentricity
	std::vector<Double_t> e_pp;
	
	for(int j = 0; j < s_size; j++)
	{
		Double_t valuesepp = sqrt((sqsigmax[j] - sqsigmay[j]) * (sqsigmax[j] - sqsigmay[j]) + 4 * (cov[j] * cov[j]))/(sqsigmax[j] + sqsigmay[j]);
		e_pp.push_back(valuesepp);
	
	// Priting the eccentricity of erp and epp
	}
	for(int k = 0; k < s_size; k++)
	{
		std::cout << "Erp: " << e_rp[k] << " Epp: "<< e_pp[k] << std::endl;
	}
	
	// Overlap area -> Reaction plane (s_rp) and participant plane(s_pp)

	std::vector<Double_t> s_rp;
	std::vector<Double_t> s_pp;
	
	for(int l = 0; l < s_size; l++)
	{
		Double_t values_rp = pi * sqrt(sqsigmax[l] * sqsigmay[l]);
		s_rp.push_back(values_rp);
		
		Double_t values_pp = pi * sqrt(sqsigmax[l] * sqsigmay[l] - cov[l]);
		s_pp.push_back(values_pp);
	}
		
	
	// Draw a histogram to see the distribuition of the data
	TH1F* histogram = new TH1F("h1", "histograma", 50, 0, 1);
	TCanvas* canvas = new TCanvas();
	for(int i = 0; i < s_size; i++)
	{
		histogram->Fill(e_pp[i]);
	}
	histogram->Draw();
	
	
	
	
	if(debugging == 1)
	{
	// FOR DEBUGING
	
	// DEBUG PURPOSES -> VERIFY IF XMEANSQ AND YMEANSQ HAVE THE SAME SIZE
	std::cerr << "Tamanho de xmeansq: " << sqsigmax.size() << std::endl;
	std::cerr << "Tamanho de ymeansq: " << sqsigmay.size() << std::endl;
	
	// VERIFY IF NO NUCLEUS HAS MORE THAN 416 COLLISIONS (FOR LEAD)
	int erros_x = 0;
	int erros_y = 0;
	for(int i = 0; i < xPos.size(); i++)
	{
		if(xPos[i].size() > 416)
		{
			erros_x += 1;
		}
		if(yPos[i].size() > 416)
		{
			erros_y += 1;
		}
		
		
	}
	
	// 
	for(int m = 0; m < s_size; m++)
	{
		std::cerr << xPos[m].size() << endl;
	}
	
	
	}
	
	
}
