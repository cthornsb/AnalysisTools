// TimeAlign.cpp
// C. Thornsberry
// Aug. 26th, 2014
// Scan the TOF for each detector and align the peaks to zero ns
// SYNTAX: ./TimeAlign {root_filename} {root_treename} {time_shift_filename}

#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH2D.h"
#include "TH1D.h"

#include <iostream>
#include <fstream>
#include <sstream>

struct BinContent{
	double content;
	double center;
	
	BinContent(double center_, double content_){
		content = content_; center = center_;
	}
};

// Make a slice along the x-axis of a 2D histogram
// Histograms must have the same number of bins on the x-axis
// _project_x -------------------------------------------------------
bool MakeSliceX(TH2D* hist2D, double low, double high, unsigned int &count, std::vector<BinContent> &contents){
	if(!hist2D){ return false; } // 1D or 2D not initialized
	double range_min = hist2D->GetYaxis()->GetXmin();
	double range_max = hist2D->GetYaxis()->GetXmax();
	int num_bins = hist2D->GetNbinsX();
	if(!(low <= high && low >= range_min && high <= range_max)){ return false; } // Values out of range
	contents.clear();
	
	int low_bin, high_bin, x_bin, z_bin;
	hist2D->GetBinXYZ(hist2D->FindBin(0.0, low), x_bin, low_bin, z_bin);
	hist2D->GetBinXYZ(hist2D->FindBin(0.0, high), x_bin, high_bin, z_bin);
	
	// Make the slice along the x-axis
	count = 0;

	for(int i = low_bin; i <= high_bin; i++){
		for(int j = 1; j <= num_bins; j++){ 
			contents.push_back(BinContent(hist2D->GetXaxis()->GetBinCenter(j) ,hist2D->GetBinContent(j, i)));
			count += hist2D->GetBinContent(j, i);
		}
	}
	
	return true;
}

// For compilation
int main(int argc, char *argv[]){
	if(argc < 4){ 
		std::cout << " Error! Invalid number of arguments. Expected 3, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./TimeAlign {filename} {treename} {shift}\n";
		return 1; 
	}

	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);

	TFile *file = new TFile(argv[1], "READ");
	if(file->IsZombie()){
		std::cout << " Error! Failed to load input file '" << argv[1] << "'\n";
		return 1;
	}

	TTree *tree = (TTree*)file->Get(argv[2]);
	if(!tree){
		std::cout << " Error! Failed to load input tree '" << argv[2] << "'\n";
		file->Close();
		return 1;
	}
	tree->SetMakeClass(1);
	
	std::ofstream output(argv[3]);
	if(!output.good()){ 
		if(argc > 3){ std::cout << " Error! Failed to load output file '" << argv[3] << "'\n"; }
		else{ std::cout << " Error! Failed to load output file 'time_shifts.dat'\n"; }
		file->Close();
		return 1; 
	}

	double max_x, max_y;

	std::vector<BinContent> hist_contents;
	std::vector<double> vandle_tof;
	std::vector<int> vandle_loc;
	TBranch *b_vandle_tof, *b_vandle_loc;
	
	tree->SetBranchAddress("vandle_TOF", &vandle_tof, &b_vandle_tof);
	tree->SetBranchAddress("vandle_loc", &vandle_loc, &b_vandle_loc);

	if(!b_vandle_tof){
		std::cout << " Error! Failed to load input branch for vandle_TOF\n";
		file->Close();
		return 1;
	}
	if(!b_vandle_loc){
		std::cout << " Error! Failed to load input branch for vandle_loc\n";
		file->Close();
		return 1;
	}

	//TCanvas *can1D = new TCanvas("can1D", "1D Canvas");
	TCanvas *can2D = new TCanvas("can2D", "2D Canvas");
	
	std::vector<double>::iterator tof_iter;
	std::vector<int>::iterator loc_iter;
	
	TH2D *hist2 = new TH2D("hist2", "TOF vs. Vandle Bar", 500, -400, 100, 42, 0, 42);
	hist2->SetStats(false);

	// Fill the 2D histogram
	std::cout << " Filling 2D histogram... This may take some time\n";
	unsigned int num_entries = tree->GetEntries();
	unsigned int chunk = (unsigned int)(num_entries*0.1);
	unsigned short percent = 0;
	for(unsigned int i = 0; i < num_entries; i++){
		if(i % chunk == 0){ 
			std::cout << "  Entry number " << i << " (" << percent << "%)\n"; 
			percent += 10;
		}
		tree->GetEntry(i);
		for(tof_iter = vandle_tof.begin(), loc_iter = vandle_loc.begin(); 
		  tof_iter != vandle_tof.end() && loc_iter != vandle_loc.end(); tof_iter++, loc_iter++){
			hist2->Fill(*tof_iter, *loc_iter);
		}
	}
	std::cout << " Filled " << hist2->GetEntries() << " entries into 2D histogram\n";

	// Plot the 2D histogram
	can2D->cd();
	hist2->Draw("COLZ");
	can2D->Update();

	// Find bar maxima
	unsigned int entries;
	for(unsigned int i = 0; i < 42; i++){
		std::cout << " Working on bar number " << i << "... ";
		if(MakeSliceX(hist2, (double)i, (double)i, entries, hist_contents)){
			std::cout << "(" << entries << " entries)\n";
			max_y = -10000.0;
			for(std::vector<BinContent>::iterator iter = hist_contents.begin(); iter != hist_contents.end(); iter++){
				if(iter->content > max_y){ 
					max_y = iter->content;
					max_x = iter->center;
				}
			}
			output << i << "\t" << max_x << std::endl;
		}
		else{ std::cout << "failed!\n"; }
	}
	
	// Cleanup
	can2D->Close();
	file->Close();
	
	delete file;
	
	output.close();
	
	rootapp->Delete();
	
	return 0;
}

// For cint
int TimeAlign(int argc, char* argv[]){ return main(argc, argv); }
