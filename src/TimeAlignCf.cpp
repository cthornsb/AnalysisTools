// TimeAlignCf.cpp
// C. Thornsberry
// Aug. 26th, 2014
// Scan the TOF for each detector and align points input by the user to zero ns
// SYNTAX: ./TimeAlignCf {root_filename} {root_treename} {time_shift_filename}

#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH2D.h"
#include "TH1D.h"

#include "TMarker.h"

#include <iostream>
#include <fstream>
#include <sstream>

// Make a slice along the x-axis of a 2D histogram
// Histograms must have the same number of bins on the x-axis
// _project_x -------------------------------------------------------
bool MakeSliceX(TH2D* hist2D, TH1D* hist1D, double low, double high, unsigned int &count){
	if(!hist2D || !hist1D){ return false; } // 1D or 2D not initialized
	double range_min = hist2D->GetYaxis()->GetXmin();
	double range_max = hist2D->GetYaxis()->GetXmax();
	int num_bins = hist2D->GetNbinsX();
	if(!(low <= high && low >= range_min && high <= range_max)){ return false; } // Values out of range
	
	int low_bin, high_bin, x_bin, z_bin;
	hist2D->GetBinXYZ(hist2D->FindBin(0.0, low), x_bin, low_bin, z_bin);
	hist2D->GetBinXYZ(hist2D->FindBin(0.0, high), x_bin, high_bin, z_bin);
	
	// Zero the projection histogram and ensure the ranges are the same
	hist1D->SetBins(num_bins, hist2D->GetXaxis()->GetXmin(), hist2D->GetXaxis()->GetXmax());
	
	// Make the slice along the x-axis
	count = 0;
	double bin_content[num_bins];
	
	for(int i = 0; i < num_bins; i++){ bin_content[i] = 0.0; }
	for(int i = low_bin; i <= high_bin; i++){
		for(int j = 1; j <= num_bins; j++){ 
			bin_content[j] += hist2D->GetBinContent(j, i);
		}
	}

	for(int i = 1; i <= num_bins; i++){
		hist1D->SetBinContent(i, bin_content[i]); 
		count += (unsigned int)bin_content[i];
	}
	
	return true;
}

// For compilation
int main(int argc, char *argv[]){
	if(argc < 4){ 
		std::cout << " Error! Invalid number of arguments. Expected 3, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./TimeAlignCf {filename} {treename} {shift}\n";
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
	
	double shifts[42];
	double bar, shift, flash;
	unsigned short num_shifts = 0;

	std::ofstream output;
	if(argc > 3){ output.open(argv[3]); }
	else{ output.open("time_shifts.dat"); }
	if(!output.good()){ 
		if(argc > 3){ std::cout << " Error! Failed to load output file '" << argv[3] << "'\n"; }
		else{ std::cout << " Error! Failed to load output file 'time_shifts.dat'\n"; }
		file->Close();
		return 1; 
	}

	double parameters[3]; 
	parameters[0] = 0.0;
	parameters[1] = 0.0;
	parameters[2] = 0.0;

	double maximum;
	double max_x, max_y;

	std::vector<double> vandle_tof;
	std::vector<int> vandle_loc;
	TBranch *b_vandle_tof, *b_vandle_loc;
	TMarker *marker;
	
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

	TCanvas *can1D = new TCanvas("can1D", "1D Canvas");
	TCanvas *can2D = new TCanvas("can2D", "2D Canvas");
	
	std::vector<double>::iterator tof_iter;
	std::vector<int>::iterator loc_iter;
	
	TH1D *hist1 = new TH1D("hist1", "hist1", 1000, -400, 100); // Histogram for slices
	TH2D *hist2 = new TH2D("hist2", "TOF vs. Vandle Bar", 1000, -400, 100, 42, 0, 42); // Main 2d histogram
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
	can1D->cd();

	// Find bar maxima
	unsigned int entries;
	for(unsigned int i = 0; i < 42; i++){
		std::cout << " Working on bar number " << i << "... ";
		if(MakeSliceX(hist2, hist1, (double)i, (double)i, entries)){
			std::cout << "(" << entries << " entries)\n";
			std::cout << "  Mark neutron peak maximum...\n";
			
			std::stringstream stream;
			stream << "Vandle Bar " << i;
			
			hist1->SetTitle(stream.str().c_str());
			hist1->Draw();
			can1D->Update();
			
			marker = (TMarker*)can1D->WaitPrimitive("TMarker");
			max_x = marker->GetX();
			max_y = marker->GetY();
			marker->Delete();
			
			output << i << "\t" << max_x << std::endl;
		}
		else{ std::cout << "failed!\n"; }
	}
	
	// Cleanup
	can2D->Close();
	file->Close();
	
	delete file;
	
	output.close();
	
	return 0;
}

// For cint
int TimeAlign(int argc, char* argv[]){ return main(argc, argv); }
