// CheckAlign.cpp
// C. Thornsberry
// Aug. 26th, 2014
// Check the time alignment for a run
// SYNTAX: ./CheckAlign {root_filename} {root_treename} {time_shift_filename}

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

// For compilation
int main(int argc, char* argv[]){
	if(argc < 4){ 
		std::cout << " Error! Invalid number of arguments. Expected 3, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./CheckAlign {filename} {treename} {shift}\n";
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
	
	bool use_shifts = true;
	if(strcmp(argv[3], "none") != 0){ // Load time shifts
		std::ifstream input;
		if(argc > 3){ input.open(argv[3]); }
		else{ input.open("time_shifts.dat"); }
		if(!input.good()){ 
			if(argc > 3){ std::cout << " Error! Failed to load output file '" << argv[3] << "'\n"; }
			else{ std::cout << " Error! Failed to load output file 'time_shifts.dat'\n"; }
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
		
			if(num_shifts >= 42){ break; }
			if(input.eof()){ break; }
		}
		input.close();
		std::cout << " Finished loading time shifts\n";
	}
	else{ use_shifts = false; } // Not using the shifts

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

	TCanvas *can2D = new TCanvas("can2D", "2D Canvas");
	TH2D *hist2 = NULL;
	if(use_shifts){ hist2 = new TH2D("hist2", "TOF vs. Vandle Bar", 500, -400-shifts[0], 100-shifts[0], 42, 0, 42); }
	else{ hist2 = new TH2D("hist2", "TOF vs. Vandle Bar", 1000, -400, 100, 42, 0, 42); }
	hist2->GetXaxis()->SetTitle("TOF (0.5 ns/bin)");
	hist2->GetYaxis()->SetTitle("Vandle Bar Number");
	hist2->SetStats(false);
	
	std::vector<double>::iterator tof_iter;
	std::vector<int>::iterator loc_iter;

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
			if(use_shifts){ hist2->Fill(*tof_iter - shifts[*loc_iter], *loc_iter); }
			else{ hist2->Fill(*tof_iter, *loc_iter); }
		}
	}
	std::cout << " Filled " << hist2->GetEntries() << " entries into 2D histogram\n";

	// Plot the 2D histogram
	can2D->cd();
	hist2->Draw("COLZ");
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
