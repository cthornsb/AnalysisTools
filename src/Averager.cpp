// Cory R. Thornsberry
// Sept. 24th, 2014
// file: Averager.cpp
// Average many scintillator pulses together and view the result

#include "Loader.h"

#include <iostream>
#include <stdlib.h>
#include <vector>

#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TGraph.h"

// For compilation
int main(int argc, char* argv[]){
	if(argc < 4){
		std::cout << " Error! Invalid number of arguments. Expected 3, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./Averager {filename} {treename} {branch} [method#] [debug]\n";
		return 1;
	}
		
	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	gSystem->Load("libTree");
	
	bool debug = false;
	if(argc > 5 && strcmp(argv[5], "debug") == 0){ 
		debug = true; 
	}

	// Branch variables
	PixieLoader loader(argv[1], argv[2]);
	TriggerStructure *trig_struct = loader.GetTrigger();
	TriggerWaveform *trig_wave = loader.GetTriggerWave();

	// Get the pulse size
	unsigned int wave_size = 0;
	for(unsigned int i = 0; i < loader.GetEntries(); i++){
		loader.GetEntry(i);
		if(trig_struct->trigger_mult == 0){ continue; }
		else{ 
			wave_size = trig_wave->trigger_wave.size()/trig_struct->trigger_mult;
			break; 
		}
	}
	std::cout << " Using wave size " << wave_size << std::endl;

	unsigned int num_entries = loader.GetEntries();
	unsigned int count = 0, index = 0;
	
	int sum[wave_size];
	for(unsigned int i = 0; i < wave_size; i++){ sum[i] = 0; }	
	
	std::vector<int>::iterator start, stop, iter;
	std::cout << " Processing " << num_entries << " entries...\n";
	for(unsigned int i = 0; i < num_entries; i++){
		loader.GetEntry(i);
		if(i+1 % 1000000 == 0){ std::cout << "  Entry no. " << i+1 << std::endl; }
		if(trig_struct->trigger_mult == 0){ continue; }
		else{
			//for(unsigned int pulse = 0; pulse < trig_struct->trigger_mult; pulse++){
			for(unsigned int pulse = 0; pulse < 1; pulse++){
				start = trig_wave->trigger_wave.begin() + (pulse * wave_size);
				stop = start + wave_size;
				index = 0;
				for(iter = start; iter != stop; iter++){
					sum[index] += *iter;
					index++;
				}
				count++;
			}
		}
	}

	// Variables for TGraph
	const double *x_val = new double[wave_size];
	const double *y_val = new double[wave_size];
	TCanvas *can = new TCanvas("canvas", "canvas");
	TGraph *graph = new TGraph(wave_size, x_val, y_val);
	graph->SetLineColor(4);
	graph->SetLineWidth(2);
	
	for(unsigned int i = 0; i < wave_size; i++){
		graph->SetPoint(i, i, ((double)sum[i])/count);
	}

	std::cout << " Done, found " << count << " individual pulses\n";

	can->cd();
	graph->Draw();
	can->Update();
	can->WaitPrimitive();

	can->Close();

	return 0;	
}
