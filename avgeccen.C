#include "globals.h"

// Calculate standard deviation
vec_1d get_sigma(vec_2d pos)
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

// Calculate covariance of x and y
vec_1d get_covarience(vec_2d xpos, vec_2d ypos)
{
	std::vector<Double_t> covariance;
	for(int i = 0; i < xpos.size(); i++)
	{
		Double_t 
			sum1 	  = 0,
			sum2 	  = 0,
			prodsum   = 0;
		
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
		covariance.push_back(cov);
	}
	
	return covariance;
}
	
		
void avgeccen()

// Creating the vectors to store the geometric values
{	
vec_2d
	v_erp,
	v_epp,
	v_srp,
	v_spp,
	v_par;

// Calling the function of collider for the values of the impact parameter (b) from collider.C	
vec_1d b  = values_b();

        // Calculating the square root of the variance
	int len = static_cast<int>(xPosSim.size());
	for(int iterations = 0; iterations < len; iterations++)
	{
		vec_1d sqsigmax = get_sigma(xPosSim[iterations]);
		vec_1d sqsigmay = get_sigma(yPosSim[iterations]);
		vec_1d cov      = get_covarience(xPosSim[iterations], yPosSim[iterations]);
		
		
		int s_size = static_cast<int>(sqsigmax.size());
		
		
		// Calculating the eccentricities of reaction and participants planes
		vec_1d e_rp, e_pp;
		
		for(int j = 0; j < s_size; j++)
		{
			Double_t valuesepp = TMath::Sqrt(abs((sqsigmax[j] - sqsigmay[j]) * (sqsigmax[j] - sqsigmay[j]) + 4 * (cov[j] * cov[j])))/(sqsigmax[j] + sqsigmay[j]);
			e_pp.push_back(valuesepp);
			
			Double_t valueserp = (sqsigmay[j] - sqsigmax[j])/(sqsigmay[j] + sqsigmax[j]);
			e_rp.push_back(valueserp);
	      
		}
		
		v_erp.push_back(e_rp);
		v_epp.push_back(e_pp);
		
		// Creating vectors of the reaction and participant planes
		vec_1d s_rp;
		vec_1d s_pp;
		
		for(int l = 0; l < s_size; l++)
		{
			Double_t values_rp = pi * TMath::Sqrt(sqsigmax[l] * sqsigmay[l]);
			s_rp.push_back(values_rp);
			
			Double_t values_pp = pi * sqrt(sqsigmax[l] * sqsigmay[l] - cov[l]);
			s_pp.push_back(values_pp);
		}
		v_srp.push_back(s_rp);
		v_spp.push_back(s_pp);
		
		vec_1d v_p;
		int teste;
		for(int m = 0; m < xPosSim[iterations].size(); m++)
		{
			 int value = xPosSim[iterations][m].size();
			 v_p.push_back(value);
		}
		
		v_par.push_back(v_p);
	}
	
	// Generate the means to plot in average eccentricity vs average participants
	vec_1d 
		plot_erp,
		plot_epp;
		
	Double_t len_values = static_cast<int>(v_erp.size());
	for(int i = 0; i < len_values; i++)
	{
		Double_t media_rp = TMath::Mean(v_erp[i].size(), v_erp[i].data());
		plot_erp.push_back(media_rp);
		Double_t media_pp = TMath::Mean(v_epp[i].size(), v_epp[i].data());
		plot_epp.push_back(media_pp);
		
	}
	
	// Get the mean of the participants
	vec_1d plot_partNumber;
	for(int k = 0; k < v_par.size(); k++)
	{
		Double_t media_rp = TMath::Mean(v_par[k].size(), v_par[k].data());
		plot_partNumber.push_back(media_rp);
	}
	
	// Creating a 2d histogram for the overlap area
	// Recebe no eixo y -> ecentricidade do RP
	// Recebe no eixo x -> quantidade de participantes
	// Recebe na cor    -> Frequência
	
	
	TH2F* part = new TH2F("part", "overlap; eixox; eixoy",50, 0, 408, 50, 0, 40);
	
	std::cout << v_par.size() << " e " << v_srp.size();
	for(int s = 0; s < v_par.size(); s++)
	{
		for(int m = 0; m < v_par[s].size(); m++)
		{
			std::cout << v_par[s][m];
			part->Fill(v_par[s][m], v_srp[s][m]);
		}
	}
	TCanvas* ctest = new TCanvas("ctest", "overlap", 800, 600);
	gStyle->SetPalette(kRainBow);
	part->Draw("COLZ");
	
	
	
	// Plotting average eccentricities in function of impact parameter
	TCanvas* c1        = new TCanvas();
	TMultiGraph* means = new TMultiGraph();
	TLegend* legend_1    = new TLegend(0.75, 0.79, 0.88, 0.88);
	
	TGraph* graph1     = new TGraph(b.size(), &(b[0]), &(plot_erp[0]));
	graph1->SetMarkerStyle(kFullCircle);
	graph1->SetMarkerColor(2);
	
	TGraph* graph2     = new TGraph(b.size(), &(b[0]), &(plot_epp[0]));
	graph2->SetMarkerStyle(kFullCircle);
	graph2->SetMarkerColor(4);
	
	legend_1->AddEntry(graph1, "<#epsilon_{rp}>", "p");
	legend_1->AddEntry(graph2, "<#epsilon_{pp}>", "p");
	legend_1->SetHeader("Legend", "C");
	
	means->Add(graph1);
	means->Add(graph2);
	means->SetTitle("Average eccentricity as function of impact parameter (Pb+Pb)"); // Note that this may differ depending if the nucleons change
	means->GetXaxis()->SetTitle("b (fm)");
	means->GetXaxis()->CenterTitle(true);
	means->GetXaxis()->SetLimits(0, 14);
	means->GetYaxis()->SetRangeUser(-0.1, 1);
	means->GetYaxis()->SetTitle("<#epsilon>");
	means->GetYaxis()->CenterTitle(true);
	means->Draw("AP");
	legend_1->Draw();
	
	// Plotting average eccentricities in function of the mean of the participants
	TCanvas* c2       = new TCanvas();
	TLegend* legend_2    = new TLegend(0.75, 0.79, 0.88, 0.88);
	
	TMultiGraph* pmg = new TMultiGraph();
	TGraph* graph4    = new TGraph(plot_partNumber.size(), &(plot_partNumber[0]), &(plot_erp[0]));
	graph4->SetMarkerStyle(kFullCircle);
	graph4->SetMarkerColor(2);
	
	TGraph* graph5    = new TGraph(plot_partNumber.size(), &(plot_partNumber[0]), &(plot_epp[0]));
	graph5->SetMarkerStyle(kFullCircle);
	graph5->SetMarkerColor(4);
	
	legend_2->AddEntry(graph1, "<#epsilon_{rp}>", "p");
	legend_2->AddEntry(graph2, "<#epsilon_{pp}>", "p");
	legend_2->SetHeader("Legend", "C");
	
	pmg->Add(graph4);
	pmg->Add(graph5);
	pmg->SetTitle("Average eccentricity as function of average participants (Pb+Pb)");
	pmg->GetXaxis()->SetTitle("<N_{part}>");
	pmg->GetXaxis()->CenterTitle(true);
	pmg->GetXaxis()->SetLimits(0, 408);
	pmg->GetYaxis()->SetRangeUser(-0.1, 1);
	pmg->GetYaxis()->SetTitle("<#epsilon>");
	pmg->GetYaxis()->CenterTitle(true);
	pmg->Draw("AP");
	legend_2->Draw();
	
	c1->SaveAs("eccenVSb.png");
	c2->SaveAs("eccenVSNpart.png");
	
// CC

	/*TF1* NBD = new TF1("NBD", "TMath::Gamma(x+[1])/(TMath::Gamma(x)*TMath::Gamma([1])) *"
                             "(([0]/[1])**(x))/(([0]/[1]+1)**(x+[1]))", 0, 0.1);
                             
	TH1F* hist = new TH1F("hist", "Distribuicao;Valor", 100, 0.0, 12);
	NBD->SetParameters(1.2, 0.77);

	// 
	int b_len = static_cast<int>(b.size());
	std::cout << b_len << endl;
	for(int test = 0; test < b_len; test++)
	{
		for(int test2 = 0; test2 < static_cast<int>(xPosSim[test].size()); test2++)
		{
			int NColl = 0.5 * xPosSim[test][test2].size();
			//std::cout << NColl << endl;
			double sumET = 0;
			for(int k = 0; k <= NColl; k++){
				sumET += NBD->GetRandom();
			}
			std::cout << sumET << endl;
			hist->Fill(sumET * 0.041);
		}
	}

	TCanvas* canvas = new TCanvas("canvas", "Histograma", 800, 600);
	hist->Draw();

	// Salvar histograma em um arquivo
	canvas->SaveAs("histograma.png");

	// Liberar memória
	delete hist;
	delete NBD;
	delete canvas; */
		


	
	
}
	
	
