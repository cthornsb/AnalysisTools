// Cory R. Thornsberry
// Sept. 23rd, 2014
// file: Test.cpp
// Perform PSD on scintillator pulses and separate neutrons from gammas
// Gates on neutron pulses, and outputs the entries to a new root file

#include "Analysis.h"
#include "Loader.h"

#include <iostream>
#include <stdlib.h>
#include <time.h>

#include "TApplication.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TCutG.h"

// Align the liquid_TOF using shift values from a file
// Only call this function after calling PixieLoader::FillEntry(#, false)
void CorrectTOF(LiquidStructure *structure, std::vector<double> &shifts){
	std::vector<double>::iterator tof_iter = structure->liquid_TOF.begin();
	std::vector<int>::iterator loc_iter = structure->liquid_loc.begin();
	for(tof_iter, loc_iter; tof_iter != structure->liquid_TOF.end() && loc_iter != structure->liquid_loc.end(); tof_iter++, loc_iter++){
		*tof_iter -= shifts.at(*loc_iter);
	}
}

// For compilation
int main(int argc, char* argv[]){
	if(argc < 5){
		std::cout << " Error! Invalid number of arguments. Expected 4, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./Test {filename} {treename} {branch} {shift_fname} [method#] [debug]\n";
		return 1;
	}
		
	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	gSystem->Load("libTree");
	
	unsigned short method = 0;
	if(argc >= 6){
		method = atol(argv[5]);
		if(method > 2){
			std::cout << " Encountered undefined method (" << method << "), aborting\n";
			return 1;
		}	
	}
	std::cout << " Using analysis method " << method << "\n";
	
	bool debug = false;
	if(argc > 6 && strcmp(argv[6], "debug") == 0){ 
		debug = true; 
	}

	// TOF shifts
	std::vector<double> shifts;
	double bar, shift, flash;
	bool use_shifts = true;
		
	std::ifstream input;
	if(strcmp(argv[4], "none") != 0){
		input.open(argv[4]);
		if(!input.good()){ 
			std::cout << " Error! Failed to load input file '" << argv[4] << "'\n";
			return 1; 
		}
		while(true){
			input >> bar >> shift;
			if(shifts.size() == 0){ 
				flash = shift; 
				shifts.push_back(shift);
			}
			else{ shifts.push_back(flash + (shift - flash)); }
			if(input.eof()){ break; }
		}
		input.close();
		std::cout << " Finished loading time shifts\n";
	}
	else{ use_shifts = false; }

	// Branch variables
	PixieLoader loader(argv[1], argv[2]);
	if(!loader.IsLoaded()){ 
		std::cout << " Exiting...\n";
		return 1; 
	}
	
	double short_value, long_value;
	unsigned int index;
	
	// This is a memory resident tree for storing S and L
	TTree *mem_tree = new TTree(argv[2], "Pulse analysis tree");
	mem_tree->SetDirectory(0); // Tell root this is a memory resident tree
	mem_tree->Branch("index", &index);
	mem_tree->Branch("S", &short_value);
	mem_tree->Branch("L", &long_value);

	LiquidStructure *out_struct = loader.GetLiquidOut();
	LiquidStructure *liq_struct = loader.GetLiquid();
	LiquidWaveform *liq_wave = loader.GetLiquidWave();
	
	TH2F *hist = new TH2F("hist", "Long vs. Short", 750, 0, 15000, 225, -500, 4000);

	// Get the pulse size
	unsigned int wave_size = 0;
	for(unsigned int i = 0; i < loader.GetEntries(); i++){
		loader.GetEntry(i);
		if(liq_struct->liquid_mult == 0){ continue; }
		else{ 
			wave_size = liq_wave->liquid_wave.size()/liq_struct->liquid_mult;
			break; 
		}
	}
	std::cout << " Using wave size " << wave_size << std::endl;

	// Timing variables
	long time_holder1 = 0;
	double time_holder2 = 0.0;
	clock_t cpu_time = clock();
	time_t real_time;
	time(&real_time);
	unsigned int status_update = 1000000;
	unsigned int num_entries = loader.GetEntries();
	if(debug){ num_entries = 100; }
	else{ num_entries = 1000000; }

	// Analysis variables
	PulseAnalysis *analysis = new PulseAnalysis();
	analysis->SetPulseSize(wave_size);
	if(debug){ 
		std::cout << " DEBUGGING!\n";
		analysis->ToggleDebug(); 
		loader.Test();
	}
	
	std::cout << " Processing " << num_entries << " entries...\n";
	for(unsigned int i = 0; i < num_entries; i++){
		loader.GetEntry(i);
		if(i % status_update == 0){ 
			time_holder1 += (long)(clock()-cpu_time);
			time_holder2 += difftime(time(NULL), real_time);
			cpu_time = clock();
			time(&real_time);
			if(i != 0){ 
				if(time_holder2 < 1){ time_holder2 = 1; } // Prevent zero time remaining
				std::cout << " Entry no. " << i << ", CPU = " << ((float)time_holder1)/CLOCKS_PER_SEC << " s, REAL = " << time_holder2;
				std::cout << " s, REMAIN = " << ((float)(num_entries-i))*(time_holder2/i) << " s\n";
			}
		}
		if(liq_struct->liquid_mult != 1){ continue; } // CHANGED! old: liq_struct->liquid_mult == 0
		else{
			index = i;
			for(unsigned int pulse = 0; pulse < liq_struct->liquid_mult; pulse++){
				analysis->PreProcess(liq_wave->liquid_wave);
				analysis->Process(pulse*wave_size, (pulse+1)*wave_size);
				analysis->PSD_Integration(method, long_value, short_value);
				hist->Fill(long_value, short_value);
				mem_tree->Fill();
			}
		}
	}
	std::cout << "\n Done, found " << mem_tree->GetEntries() << " pulse events\n";
	std::cout << "  Make the graphical cut on the neutrons\n\n";
	
	TCanvas *can = new TCanvas("can", "Canvas");
	can->cd();
	hist->Draw("COLZ");
	TCutG *cut = (TCutG*)can->WaitPrimitive("CUTG");

	// This is the main output file for storing gated neutron signals
	TFile *out_file = new TFile("test.root", "RECREATE");
	TTree *out_tree = new TTree(argv[2], "Pulse analysis tree");
	out_tree->Branch("S", &short_value);
	out_tree->Branch("L", &long_value);
	loader.SetOutputBranches("11001000"); // Turn off unnecessary branches to save space and time
	loader.SetTree(out_tree);
	
	std::cout << " Working...\n";
	for(unsigned int i = 0; i < mem_tree->GetEntries(); i++){
		mem_tree->GetEntry(i);
		if(cut->IsInside(long_value, short_value)){
			loader.FillEntry(index, false); // Load the values, but do not fill the tree			
			CorrectTOF(out_struct, shifts); // Time align the liquid TOFs
			out_tree->Fill(); // Manually fill the tree
		}
	}
	std::cout << "  Done, found " << out_tree->GetEntries() << " pulse events inside the cut\n";
	
	out_file->cd();
	out_tree->Write();
	out_file->Close();

	std::cout << "  Wrote 'test.root' to file\n\n";	
	
	can->Close();
	mem_tree->Delete();
	
	rootapp->Delete();
	delete analysis;
	delete out_file;
	
	return 0;
}
