

void overlap()
{
  // Extracting the information collected from collide.cpp
  TFile *f = new TFile("data.root");
  TTree *t = (TTree*)f->Get("Data");

  std::vector<double> *xCol = nullptr;
  std::vector<double> *yCol = nullptr;
  int NPart;
  t->SetBranchAddress("xCol", &xCol);
  t->SetBranchAddress("yCol", &yCol);
  t->SetBranchAddress("NPart", &NPart);

	std::vector<double> Spp;
	for(int i = 0; i < t->GetEntries(); i++)
	{
		t->GetEntry(i);
		double x2 = 0;
		double y2 = 0;
		double xy = 0;
		for(int k = 0; k < xCol->size(); k++)
		{
			x2 += (*xCol)[k] * (*xCol)[k];
			y2 += (*yCol)[k] * (*yCol)[k];
			xy += (*xCol)[k] * (*yCol)[k];
			
		}
		double x_mean = TMath::Mean(xCol->size(), xCol->data());
		double y_mean = TMath::Mean(yCol->size(), yCol->data());
		double covariance = xy/xCol->size() - x_mean * y_mean;
		Spp.push_back(TMath::Pi() * TMath::Sqrt(abs((x2/xCol->size() - x_mean * x_mean) * (y2/yCol->size() - y_mean * y_mean) - covariance * covariance)));
	}
		
	// Plotting the 2d histogram of the overlap area and number of participants	
	TH2F* OverlapArea = new TH2F("Area", "OverlapArea", 10, 0, 416, 10, 0, 40);
	for(int m = 0; m < t->GetEntries(); m++)
	{
	t->GetEntry(m);
		OverlapArea->Fill(NPart, Spp[m]);
	}
	
	TCanvas* c1 = new TCanvas("c1", "OverlapArea", 800, 600);
	gStyle->SetPalette(kRainBow);
	OverlapArea->Draw("COLZ");
	
}
