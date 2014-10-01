// Gater.cpp
// C. Thornsberry
// Aug. 26th, 2014
// Gate neutron data on the TOF vs. QDC plot
// SYNTAX: ./Gater {filename} {treename} {shift_fname} {#det}

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCutG.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TGraph.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <cmath>

#include "Structures.h"
#include "Loader.h"

#define NEUTRON_MASS 939.565 // Mev/c^2
#define VANDLE_T_RES 1.0 // Time resolution for Vandle (ns)
#define THICKNESS 3.0 // Thickness of Vandle bar (cm)
#define RADIUS 50.0 // Distance from target to Vandle (cm)
#define C 3E10 // Speed of light (cm/s)

#define LEFT 16.2 // Minimum allowed ToF (ns)
#define RIGHT 114.2 // Maximum allowed ToF (ns)

#define ES_LOW 0.0	// Minimum allowed end scintillator energy (arb.)
#define ES_HIGH 1250.0 // Maximum allowed end scintillator energy (arb.)

class LimitFunc{
  private:
	double *x, *y;
	unsigned int num_points;
	
  public:
	LimitFunc(unsigned int num_points_){
		num_points = num_points_;
		x = new double[num_points];
		y = new double[num_points];
	}
	~LimitFunc(){
		delete[] x;
		delete[] y;
	}
	
	void SetPoint(unsigned int index, double x_val, double y_val){
		if(index < num_points){
			x[index] = x_val;
			y[index] = y_val;
		}
	}
	
	// Interpolate along the defined line
	// Return true if the value is in range and false otherwise
	bool Interpolate(double x_val, double &y_val){
		y_val = 0.0;
		for(unsigned int i = 0; i < num_points-1; i++){
			if(x_val >= x[i] && x_val <= x[i+1]){
				if(x_val == x[i]){ y_val = y[i]; }
				else if(x_val == x[i+1]){ y_val = y[i+1]; }
				else{ y_val = y[i] + (x_val - x[i])*(y[i+1] - y[i])/(x[i+1] - x[i]); } // Do interpolation
				return true;
			}
		}
		return false;
	}
};

double Func(double *x, double *par){
	return par[0]/(x[0]*x[0]);
}

double Func(double tof, double coeff){
	return coeff/(tof*tof);
}

// For compilation
int main(int argc, char* argv[]){
	if(argc < 5){
		std::cout << " Error! Invalid number of arguments. Expected 4, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./Gater {filename} {treename} {shift_fname} {#det}\n";
		return 1;
	}

	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	gSystem->Load("libTree");
		
	bool debug = false;
	if(argc > 5){
		if(strcmp(argv[5], "debug") == 0){ 
			debug = true; 
			std::cout << " DEBUGGING...\n";
		}
	}

	// Time shift variables
	float low_val, high_val;
	unsigned int num_dets = (unsigned)atol(argv[4]);

	std::cout << " Looking for " << num_dets << " detectors\n";
	
	double shifts[num_dets];
	double bar, shift, flash;
	unsigned short num_shifts = 0;
	
	for(unsigned short i = 0; i < num_dets; i++){ shifts[i] = 0.0; }

	// Declare the PixieLoader object
	// This object will load the file and all the branches
	PixieLoader loader(argv[1], argv[2]);
	
	// Get pointers to the data stored inside PixieLoader
	TriggerStructure *trig_structure = loader.GetTrigger();
	TriggerWaveform *trig_waveform = loader.GetTriggerWave();
	VandleStructure *van_structure = loader.GetVandle();
	
	std::ifstream input;
	input.open(argv[3]);
	if(!input.good()){ 
		std::cout << " Error! Failed to load input file '" << argv[3] << "'\n";
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

	std::vector<double> TOF, qdc, energy;
	std::vector<unsigned int> loc;
	TFile *out_file = NULL;
	TTree *vandle_tree = NULL;
	if(!debug){
		std::cout << " Opening output file 'gater.root'\n";
		out_file = new TFile("gater.root", "RECREATE");
		vandle_tree = new TTree(argv[2], "Analysis tree");
		loader.SetOutputBranches("11000010");
		loader.SetTree(vandle_tree);
	}

	// Canvas
	TCanvas *can = new TCanvas("can", "canvas");
	can->cd();

	// Histogram
	TH2D *hist = new TH2D("hist", "TOF vs. <E>", 500, -50, 200, 500, 0, 8000);
	hist->GetXaxis()->SetTitle("TOF (0.5 ns/bin)");
	hist->GetYaxis()->SetTitle("QDC (arb. units)");
	hist->SetStats(false);

	std::vector<double>::iterator iter1, iter2;
	std::vector<int>::iterator iter3;
	std::vector<double>::iterator iter4;
	unsigned int count0 = 0;
	unsigned int count1 = 0, count2 = 0;
	unsigned int count3 = 0, count4 = 0;
	
	long time_holder1;
	double time_holder2;
	clock_t cpu_time = clock();
	time_t real_time;
	
	unsigned int num_entries = loader.GetEntries();
		
	TGraph *graph1 = NULL;
	TGraph *graph2 = NULL;
	double qdc_coeff = (NEUTRON_MASS*RADIUS*RADIUS/(2.0*C*C))*(1E18)*(1E3); // Coefficient for maximum neutron qdc (keV/ns^2)
	std::cout << " QDC Coeff = " << qdc_coeff << " keV/ns^2\n";
	if(debug){
		if(num_entries > 1000000){ num_entries = 1000000; } // 1 M events is quite sufficient
		std::cout << " Processing " << num_entries << " entries\n";
		unsigned int num_points = (unsigned int)((RIGHT-LEFT)*2);
		double dE, x[num_points], y1[num_points], y2[num_points];
		
		for(unsigned int i = 0; i < num_points; i++){
			x[i] = 15.0 + 0.5*i; // 0.5 ns steps
			dE = std::sqrt((THICKNESS*THICKNESS/(RADIUS*RADIUS)) + (VANDLE_T_RES*VANDLE_T_RES/(x[i]*x[i])));
			y1[i] = (qdc_coeff/(x[i]*x[i]))*(1 - dE); // lower line
			y2[i] = (qdc_coeff/(x[i]*x[i]))*(1 + dE); // upper line
		}
		
		graph1 = new TGraph(num_points, x, y1); graph1->SetLineColor(2); graph1->SetLineWidth(2);
		graph2 = new TGraph(num_points, x, y2); graph2->SetLineColor(4); graph2->SetLineWidth(2);

		// Fill the histogram to cut on
		for(unsigned int i = 0; i < num_entries; i++){
			loader.GetEntry(i);
			for(iter1 = van_structure->vandle_TOF.begin(), iter2 = van_structure->vandle_qdc.begin(), iter3 = van_structure->vandle_loc.begin();
			  iter1 != van_structure->vandle_TOF.end() && iter2 != van_structure->vandle_qdc.end() && iter3 != van_structure->vandle_loc.end(); iter1++, iter2++, iter3++){
				hist->Fill((*iter1) - shifts[(*iter3)], (*iter2));
				count1++;
			}
		}

		std::cout << " Found " << count1 << " events\n";

		hist->Draw("COLZ");	
		graph1->Draw("SAME");
		graph2->Draw("SAME");
	}
	else{		
		std::cout << " Processing " << num_entries << " entries\n";
		unsigned int num_points = (unsigned int)((RIGHT-LEFT)*2);
		unsigned int position = 0;
		LimitFunc limit(num_points);

		double x;
		for(unsigned int i = 0; i < num_points; i++){ // 0.5 ns steps
			x = 15.0 + 0.5*i;
			limit.SetPoint(i, x, (qdc_coeff/(x*x))*(1 + std::sqrt((THICKNESS*THICKNESS/(RADIUS*RADIUS)) + (VANDLE_T_RES*VANDLE_T_RES/(x*x)))));
		}

		// Fill the histogram to cut on
		double qdc_limit;
		bool valid_event;
		for(unsigned int i = 0; i < num_entries; i++){
			loader.GetEntry(i);
			valid_event = false;
			if(i % 1000000 == 0){ // Timing stuff
				time_holder1 = (long)(clock()-cpu_time);
				time_holder2 = difftime(time(NULL), real_time);
				cpu_time = clock();
				time(&real_time);
				if(i != 0){ 
					if(time_holder2 < 1){ time_holder2 = 1; } // Prevent zero time remaining
					std::cout << " Entry no. " << i << ", CPU = " << ((float)time_holder1)/CLOCKS_PER_SEC << " s, REAL = " << time_holder2 << " s, remain = ";
					std::cout << ((float)(num_entries-i)/1000000)*time_holder2 << " s\n";
				}
			}
			count0 += van_structure->vandle_TOF.size(); 
			for(iter4 = trig_structure->trigger_energy.begin(); iter4 != trig_structure->trigger_energy.end(); iter4++){
				count1++;
				if(*iter4 > ES_LOW && *iter4 <= ES_HIGH){ // Valid recoil energy
					count2++;
					loader.Fill(false);
					loader.GetVandleOut()->Zero();
					for(iter1 = van_structure->vandle_TOF.begin(), iter2 = van_structure->vandle_qdc.begin(), iter3 = van_structure->vandle_loc.begin();
					  iter1 != van_structure->vandle_TOF.end() && iter2 != van_structure->vandle_qdc.end() && iter3 != van_structure->vandle_loc.end(); iter1++, iter2++, iter3++){
					  	count0++;
					  	count3++;
					  	shift = *iter1 - shifts[*iter3];
					  	if(limit.Interpolate(shift, qdc_limit) && (*iter2) <= qdc_limit){
							hist->Fill(shift, (*iter2));
							if(!valid_event){ valid_event = true; }
							loader.GetVandleOut()->Append((*iter3), shift, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, *iter2, 0.0);
							count4++;
						}
					}
					if(valid_event){ vandle_tree->Fill(); }
					break;
				}
			}
		}

		std::cout << " Found " << count2 << " Trigger events within trigger gate (" << 100.0*count2/count1 << "%)\n";
		std::cout << " Found " << count4 << " Vandle events within Vandle gate (" << 100.0*count4/count3 << "%)\n";
		std::cout << " Total number of Vandle events: " << count0 << "\n";
		std::cout << " Total number within both gates: " << count4 << " (" << 100.0*count4/count0 << "%)\n";
		hist->Draw("COLZ");
	}
	
	can->Update();
	can->WaitPrimitive();

	if(!debug){
		std::cout << " Wrote " << vandle_tree->GetEntries() << " entries to the Vandle tree\n";
		out_file->cd();
		vandle_tree->Write();
		out_file->Close();
	}
	else{
		graph1->Delete();
		graph2->Delete();
	}
	
	can->Close();	
	rootapp->Delete();
	
	return 0;
}

// For CINT
//int PulseAnalyzer(int argc, char* argv[]){ return main(argc, argv); }
