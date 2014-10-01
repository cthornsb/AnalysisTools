// CheckAlign.cpp
// C. Thornsberry
// Aug. 26th, 2014
// Check the time alignment for a run
// SYNTAX: ./CheckAlign {filename} {treename} {branchname} {shift_fname} {bins} {low} {high} {#det}

#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

// For compilation
int main(int argc, char* argv[]){
	if(argc < 9){ 
		std::cout << " Error! Invalid number of arguments. Expected 8, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./CheckAlign {filename} {treename} {branchname} {shift_fname} {bins} {low} {high} {#det}\n";
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

	int num_bins, num_dets;
	float low_val, high_val;
	num_bins = atol(argv[5]);
	low_val = atof(argv[6]);
	high_val = atof(argv[7]);
	num_dets = atol(argv[8]);

	std::cout << " Looking for " << num_dets << " detectors\n";
	std::cout << " Using " << num_bins << " bins in range [" << low_val << ", " << high_val << "]\n";
	
	double shifts[num_dets];
	double bar, shift, flash;
	unsigned short num_shifts = 0;
	
	for(unsigned short i = 0; i < num_dets; i++){ shifts[i] = 0.0; }
	
	std::ifstream input;
	if(strcmp(argv[4], "none") != 0){
		input.open(argv[4]);
		if(!input.good()){ 
			std::cout << " Error! Failed to load input file '" << argv[4] << "'\n";
			file->Close();
			return 1; 
		}
		while(true){
			input >> bar >> shift;
			if(num_shifts == 0){ 
				flash = shift; 
				shifts[0] = shift;
			}
			else{ shifts[num_shifts] = flash + (shift - flash); }
			num_shifts++;		
	
			if(num_shifts >= num_dets){ break; }
			if(input.eof()){ break; }
		}
		input.close();
		std::cout << " Finished loading time shifts\n";
	}

	std::vector<double> vandle_tof;
	std::vector<int> vandle_loc;
	TBranch *b_vandle_tof, *b_vandle_loc;
	
	std::stringstream stream; stream << argv[3];
	tree->SetBranchAddress((stream.str() + "_TOF").c_str(), &vandle_tof, &b_vandle_tof);
	tree->SetBranchAddress((stream.str() + "_loc").c_str(), &vandle_loc, &b_vandle_loc);

	if(!b_vandle_tof){
		std::cout << " Error! Failed to load input branch for " << stream.str() + "_TOF\n";
		file->Close();
		return 1;
	}
	if(!b_vandle_loc){
		std::cout << " Error! Failed to load input branch for " << stream.str() + "_loc\n";
		file->Close();
		return 1;
	}

	TCanvas *can2D = new TCanvas("can2D", "2D Canvas");
	TH2D *hist2 = new TH2D("hist2", "TOF vs. Detector", num_bins, low_val-shifts[0], high_val-shifts[0], num_dets, 0, num_dets);
	TH1D *hist1 = new TH1D("hist1", "TOF", num_bins, low_val-shifts[0], high_val-shifts[0]);
	hist1->GetXaxis()->SetTitle("TOF");
	hist2->GetXaxis()->SetTitle("TOF");
	hist2->GetYaxis()->SetTitle("Location");
	hist1->SetStats(false);
	hist2->SetStats(false);
	
	std::vector<double>::iterator tof_iter;
	std::vector<int>::iterator loc_iter;

	// Fill the 2D histogram
	std::cout << " Filling 2D histogram... This may take some time\n";
	unsigned int num_entries = tree->GetEntries();
	if(num_entries > 1000000){ num_entries = 1000000; } // 1 M events is quite sufficient
	
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
		    hist1->Fill(*tof_iter - shifts[*loc_iter]);
			hist2->Fill(*tof_iter - shifts[*loc_iter], *loc_iter);
		}
	}
	std::cout << " Filled " << hist2->GetEntries() << " entries into 2D histogram\n";

	// Plot the 2D histogram
	can2D->Divide(2);
	can2D->cd(1);
	hist2->Draw("COLZ");
	can2D->cd(2);
	hist1->Draw();
	can2D->Update();
	can2D->WaitPrimitive();
	
	// Cleanup
	can2D->Close();
	file->Close();
	
	delete file;
	
	rootapp->Delete();
	
	return 0;
}

// For cint
int CheckAlign(int argc, char* argv[]){ return main(argc, argv); }
