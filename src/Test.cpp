#include "Analysis.h"
#include "Loader.h"

#include <iostream>
#include <stdlib.h>
#include <time.h>

#include "TApplication.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
	
// For compilation
int main(int argc, char* argv[]){
	if(argc < 4){
		std::cout << " Error! Invalid number of arguments. Expected 4, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./Analysis {filename} {treename} {branch} [method#] [debug]\n";
		return 1;
	}
		
	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	gSystem->Load("libTree");
	
	unsigned short method = 0;
	if(argc >= 5){
		method = atol(argv[4]);
		if(method != 0){
			std::cout << " Encountered undefined method (" << method << "), aborting\n";
			return 1;
		}
		std::cout << " Using analysis method " << method << "\n";		
	}
	
	bool debug = false;
	if(argc > 5 && strcmp(argv[5], "debug") == 0){ 
		debug = true; 
		std::cout << " DEBUGGING...\n";
	}

	// Branch variables
	PixieLoader loader(argv[1], argv[2]);

	/*std::vector<int> wave;
	double TOF;
	unsigned int mult;

	TFile *file = new TFile(argv[1], "READ");
	if(file->IsZombie()){
		std::cout << " Failed to load the input file '" << argv[1] << "'\n";
		return 1;
	}
	TTree *tree = (TTree*)file->Get(argv[2]);
	if(!tree){
		std::cout << " Failed to load the input tree '" << argv[2] << "'\n";
		file->Close();
		return 1;
	}
	tree->SetMakeClass(1);
	
	std::stringstream branch_name;
	branch_name << argv[3];
	TBranch *b_wave, *b_mult, *b_TOF;
	tree->SetBranchAddress((branch_name.str()+"_wave").c_str(), &wave, &b_wave);
	tree->SetBranchAddress((branch_name.str()+"_mult").c_str(), &mult, &b_mult);
	tree->SetBranchAddress((branch_name.str()+"_TOF").c_str(), &TOF, &b_TOF);
	
	if(!b_wave){
		std::cout << " Failed to load the input branch '" << branch_name.str() << "_wave'\n";
		file->Close();
		return 1;
	}
	if(!b_mult){
		std::cout << " Failed to load the input branch '" << branch_name.str() << "_mult'\n";
		file->Close();
		return 1;
	}*/
	
	double short_value, long_value;
	TFile *out_file = new TFile("test.root", "RECREATE");
	TTree *out_tree = new TTree(argv[2], "Pulse analysis tree");
	out_tree->Branch("S", &short_value);
	out_tree->Branch("L", &long_value);
	loader.SetOutputBranches("10001010");
	loader.SetTree(out_tree);

	LiquidStructure *liq_struct = loader.GetLiquid();
	LiquidWaveform *liq_wave = loader.GetLiquidWave();

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

	// Analysis variables
	PulseAnalysis *analysis = new PulseAnalysis();
	double *baseline = new double[wave_size];
	analysis->SetPulseSize(wave_size);
	if(debug){ 
		std::cout << " DEBUGGING!\n";
		analysis->ToggleDebug(); 
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
		if(liq_struct->liquid_mult == 0){ continue; }
		else{
			for(unsigned int pulse = 0; pulse < liq_struct->liquid_mult; pulse++){
				analysis->PreProcess(liq_wave->liquid_wave);
				analysis->Process(pulse*wave_size, (pulse+1)*wave_size);
				analysis->PSD_Integration(3, long_value, short_value);
				loader.FillEntry(i);
			}
		}
	}
	
	std::cout << " Done, filled " << out_tree->GetEntries() << " entries\n";
	
	out_file->cd();
	out_tree->Write();
	out_file->Close();
	
	rootapp->Delete();
	delete analysis;
	delete baseline;
	delete out_file;
	
	return 0;
}
