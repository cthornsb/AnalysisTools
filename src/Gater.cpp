// Gater.cpp
// C. Thornsberry
// Aug. 26th, 2014
// Gate neutron data on the TOF vs. QDC plot
// SYNTAX: ./Gater {root_filename} {root_treename} {time_shift_filename} [debug]

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

#define NEUTRON_MASS 939.565 // Mev/c^2
#define VANDLE_T_RES 1.0 // ns
#define THICKNESS 3.0 // cm
#define RADIUS 50.0 // cm
#define C 3E10 // cm/s

#define LEFT 16.2 // ns
#define RIGHT 114.2 // ns

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
	if(argc < 4){
		std::cout << " Error! Invalid number of arguments. Expected 3, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./Gater {filename} {treename} {shift} [debug]\n";
		return 1;
	}

	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	gSystem->Load("libTree");
		
	bool debug = false;
	if(argc > 4){
		if(strcmp(argv[4], "debug") == 0){ 
			debug = true; 
			std::cout << " DEBUGGING...\n";
		}
	}

	// Time shift variables
	double shifts[42];
	double bar, shift, flash;
	unsigned short num_shifts = 0;

	// Branch variables
	VandleStructure v_structure;
	TriggerStructure t_structure;
	std::vector<double> input_TOF, input_qdc, input_energy;
	std::vector<unsigned int> input_loc;

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
	
	TBranch *b_TOF, *b_qdc, *b_loc; // Vandle Branches
	TBranch *b_energy; // Trigger Branches
	tree->SetBranchAddress("vandle_TOF", &input_TOF, &b_TOF);
	tree->SetBranchAddress("vandle_qdc", &input_qdc, &b_qdc);
	tree->SetBranchAddress("vandle_loc", &input_loc, &b_loc);
	tree->SetBranchAddress("trigger_energy", &input_energy, &b_energy);
	
	if(!b_TOF){
		std::cout << " Note: Failed to load the input branch 'vandle_TOF'\n";
		file->Close();
		return 1;
	}
	if(!b_qdc){
		std::cout << " Note: Failed to load the input branch 'vandle_qdc'\n";
		file->Close();
		return 1;
	}
	if(!b_loc){
		std::cout << " Note: Failed to load the input branch 'vandle_loc'\n";
		file->Close();
		return 1;
	}
	if(!b_energy){
		std::cout << " Note: Failed to load the input branch 'trigger_energy'\n";
		file->Close();
		return 1;
	}

	std::ifstream input;
	input.open(argv[3]);
	if(!input.good()){ 
		std::cout << " Error! Failed to load time shift file '" << argv[3] << "'\n";
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

	std::vector<double> TOF, qdc, energy;
	std::vector<unsigned int> loc;
	TFile *out_file = NULL;
	TTree *vandle_tree = NULL;
	if(!debug){
		std::cout << " Opening output file 'gater.root'\n";
		out_file = new TFile("gater.root", "RECREATE");

		vandle_tree = new TTree(argv[2], "Analysis tree");
		vandle_tree->Branch("Trigger", &t_structure);
		vandle_tree->Branch("Vandle", &v_structure);		
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
	std::vector<unsigned int>::iterator iter3;
	std::vector<double>::iterator iter4;
	unsigned int count1 = 0, count2 = 0;
	unsigned int count3 = 0, count4 = 0;
	
	long time_holder1;
	double time_holder2;
	clock_t cpu_time = clock();
	time_t real_time;
	
	unsigned int num_entries = tree->GetEntries();
		
	TGraph *graph1 = NULL;
	TGraph *graph2 = NULL;
	double qdc_coeff = (NEUTRON_MASS*RADIUS*RADIUS/(2.0*C*C))*(1E18)*(1E3); // Coefficient for maximum neutron qdc (keV ns^2)
	std::cout << " QDC Coeff = " << qdc_coeff << " keV/ns^2\n";
	if(debug){	
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
			tree->GetEntry(i);
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
			for(iter1 = input_TOF.begin(), iter2 = input_qdc.begin(), iter3 = input_loc.begin();
			  iter1 != input_TOF.end() && iter2 != input_qdc.end() && iter3 != input_loc.end(); iter1++, iter2++, iter3++){
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
			tree->GetEntry(i);
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
			for(iter1 = input_TOF.begin(), iter2 = input_qdc.begin(), iter3 = input_loc.begin(); 
			  iter1 != input_TOF.end() && iter2 != input_qdc.end() && iter3 != input_loc.end(); iter1++, iter2++, iter3++){
			  	shift = *iter1 - shifts[*iter3];
			  	if(limit.Interpolate(shift, qdc_limit) && (*iter2) <= qdc_limit){
			  		v_structure.Append((*iter3), (*iter1) - shifts[(*iter3)], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (*iter2), 0.0);
					hist->Fill((*iter1) - shifts[(*iter3)], (*iter2));
					if(!valid_event){ valid_event = true; }
					count2++;
				}
				count1++;
			}
			count3 += input_energy.size();
			if(valid_event){
				for(iter4 = input_energy.begin(); iter4 != input_energy.end(); iter4++){
					t_structure.Append((*iter4));
					count4++;
				}			
				vandle_tree->Fill();
				v_structure.Zero();
				t_structure.Zero();
			}
		}

		std::cout << " Found " << count1 << " Vandle events and " << count2 << " inside the gate (" << 100.0*count2/count1 << "%)\n";
		std::cout << " Found " << count3 << " Trigger events and " << count4 << " inside the gate (" << 100.0*count4/count3 << "%)\n";
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
	file->Close();
	
	rootapp->Delete();
	
	return 0;
}

// For CINT
//int PulseAnalyzer(int argc, char* argv[]){ return main(argc, argv); }
