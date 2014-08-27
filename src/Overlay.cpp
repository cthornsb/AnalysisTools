// Overlay.cpp
// C. Thornsberry
// Aug. 27th, 2014
// Program to overlay one spectrum over another on the same canvas
// SYNTAX: ./Overlay {root_filename1} {root_treename1} {root_branchname1} {root_filename2} {root_treename2} {root_branchname2} {num_bins} {low} {high} [debug]

#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH1D.h"

#include <iostream>
#include <sstream>
#include <stdlib.h>

int main(int argc, char* argv[]){
	if(argc < 10){ 
		std::cout << " Error! Invalid number of arguments. Expected 9, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./Overlay {filename1} {treename1} {branchname1} {filename2} {treename2} {branchname2} {bins} {low} {high} [debug]\n";
		return 1; 
	}

	bool debug = false;
	if(argc > 10){
		if(strcmp(argv[10], "debug") == 0){ 
			debug = true; 
			std::cout << " DEBUGGING!\n";
		}
	}

	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);

	std::stringstream f1name, f2name, t1name, t2name, b1name, b2name;

	// File 1
	TFile *file1 = new TFile(argv[1], "READ");
	f1name << argv[1];
	if(file1->IsZombie()){
		std::cout << " Error! Failed to load input file '" << f1name.str() << "'\n";
		return 1;
	}

	// Tree 1
	TTree *tree1 = (TTree*)file1->Get(argv[2]);
	t1name << argv[2];
	if(!tree1){
		std::cout << " Error! Failed to load input tree '" << t1name.str() << "'\n";
		file1->Close();
		return 1;
	}
	tree1->SetMakeClass(1);
	
	// File 2
	TFile *file2 = NULL;
	bool f2_is_open = false;
	if(strcmp(argv[1], argv[4]) != 0){
		if(strcmp(argv[4], "same") != 0){ // The spectra are from different files
			file2 = new TFile(argv[4], "READ");
			f2name << argv[4];
			if(file2->IsZombie()){
				std::cout << " Error! Failed to load input file '" << f2name.str() << "'\n";
				file1->Close();
				file2->Close();
				return 1;
			}
			f2_is_open = true;
		}
		else{ f2name << argv[1]; } // The spectra are from the same file
	}
	
	// Tree 2
	TTree *tree2 = NULL;
	bool t2_is_loaded = false;
	if(!f2_is_open){ // Load tree2 from file1
		if(strcmp(argv[2], argv[5]) != 0){ // The spectra are from different trees
			if(strcmp(argv[5], "same") != 0){ // Tree2 name is specified
				tree2 = (TTree*)file1->Get(argv[5]); 
				t2name << argv[5];
			} 
			else{ // Use Tree1 name
				tree2 = (TTree*)file1->Get(argv[2]); 
				t2name << argv[2];
			}
			t2_is_loaded = true;
		}
	}
	else{ // Load tree2 from file2
		if(strcmp(argv[5], "same") != 0){ // Tree2 name is specified
			tree2 = (TTree*)file2->Get(argv[5]); 
			t2name << argv[5];
		} 
		else{ // Use Tree1 name
			tree2 = (TTree*)file2->Get(argv[2]); 
			t2name << argv[2];
		}
		t2_is_loaded = true;
	}	
	
	// Load the branches
	std::vector<double> spec1, spec2;
	TBranch *b_spec1, *b_spec2;

	tree1->SetBranchAddress(argv[3], &spec1, &b_spec1);
	b1name << argv[3];
	
	if(t2_is_loaded){ // Load branch2 from tree2
		if(!tree2){
			std::cout << " Error! Failed to load input tree '" << t2name.str() << "'\n";
			file1->Close();
			if(f2_is_open){ file2->Close(); }
			return 1;
		}	
		tree2->SetMakeClass(1);

		if(strcmp(argv[6], "same") != 0){ 
			tree2->SetBranchAddress(argv[6], &spec2, &b_spec2); 
			b2name << argv[6];
		}
		else{ 
			tree2->SetBranchAddress(argv[3], &spec2, &b_spec2); 
			b2name << argv[3];
		}
	}
	else{ // Load branch2 from tree1
		if(strcmp(argv[6], "same") != 0){ 
			tree1->SetBranchAddress(argv[6], &spec2, &b_spec2); 
			b2name << argv[6];
		}
		else{ // Variable 2 was passed 'same same same'
			std::cout << " Error! Attempted to reset address of branch '" << b1name.str() << "' to variable 2!\n";
			file1->Close();
			if(f2_is_open){ file2->Close(); }
			return 1;
		}
	}
	
	if(!b_spec1){
		std::cout << " Error! Failed to load input branch for '" << b1name.str() << "'\n";
		file1->Close();
		if(f2_is_open){ file2->Close(); }
		return 1;
	}
	if(!b_spec2){
		std::cout << " Error! Failed to load input branch for '" << b2name.str() << "'\n";
		file1->Close();
		if(f2_is_open){ file2->Close(); }
		return 1;
	}
	
	int num_bins = atoi(argv[7]);
	double low_limit = (double)atof(argv[8]);
	double high_limit = (double)atof(argv[9]);


	TH1D *hist1 = new TH1D("hist1", b1name.str().c_str(), num_bins, low_limit, high_limit);
	hist1->GetXaxis()->SetTitle(b1name.str().c_str());
	hist1->GetYaxis()->SetTitle("Counts per bin");
	hist1->SetStats(false);
	hist1->SetLineColor(2);

	TH1D *hist2 = new TH1D("hist2", b2name.str().c_str(), num_bins, low_limit, high_limit);
	hist2->GetXaxis()->SetTitle(b2name.str().c_str());
	hist2->GetYaxis()->SetTitle("Counts per bin");
	hist2->SetStats(false);
	hist2->SetLineColor(4);
	
	std::vector<double>::iterator iter1, iter2;

	// Print info
	std::cout << "  Plotting with " << num_bins << " bins in range (" << low_limit << ", " << high_limit << ")\n";
	if(debug){
		std::cout << "   file2: ";
		if(f2_is_open){ std::cout << "Good\n"; }
		else{ std::cout << "NULL\n"; }
		std::cout << "   tree2: ";
		if(t2_is_loaded){ std::cout << "Good\n"; }
		else{ std::cout << "NULL\n"; }
		std::cout << "   Variable 1:  File = " << f1name.str() << "  Tree = " << t1name.str() << "  Branch = " << b1name.str() << std::endl;
		std::cout << "   Variable 2:  File = " << f2name.str() << "  Tree = " << t2name.str() << "  Branch = " << b2name.str() << std::endl;
	}
	
	unsigned int count1 = 0, count2 = 0;
	if(t2_is_loaded){ // Branches are from different trees
		// Iterate through tree1
		unsigned int num_entries = tree1->GetEntries();
		std::cout << "  Working on tree '" << t1name.str() << "' with " << num_entries << " entries...\n";
		for(unsigned int i = 0; i < num_entries; i++){
			tree1->GetEntry(i);
			if(i % 1000000 == 0 && i != 0){ std::cout << "   Entry number " << i << std::endl; }
			for(iter1 = spec1.begin(); iter1 != spec1.end(); iter1++){
				hist1->Fill((*iter1));
				count1++;
			}
		}

		// Iterate through tree2
		num_entries = tree2->GetEntries();
		std::cout << "  Working on tree '" << t2name.str() << "' with " << num_entries << " entries...\n";
		for(unsigned int i = 0; i < num_entries; i++){
			tree2->GetEntry(i);
			if(i % 1000000 == 0 && i != 0){ std::cout << "   Entry number " << i << std::endl; }
			for(iter2 = spec2.begin(); iter2 != spec2.end(); iter2++){
				hist2->Fill((*iter2));
				count2++;
			}
		}
	}
	else{ // Both branches are from tree1
		// Iterate through tree1
		unsigned int num_entries = tree1->GetEntries();
		std::cout << "  Working on tree '" << t1name.str() << "' with " << num_entries << " entries\n";
		for(unsigned int i = 0; i < num_entries; i++){
			tree1->GetEntry(i);
			if(i % 1000000 == 0 && i != 0){ std::cout << "   Entry number " << i << std::endl; }
			for(iter1 = spec1.begin(); iter1 != spec1.end(); iter1++){ // Vectors may have different lengths
				hist1->Fill((*iter1));
				count1++;
			}
			for(iter2 = spec2.begin(); iter2 != spec2.end(); iter2++){ // Vectors may have different lengths
				hist2->Fill((*iter2));
				count2++;
			}
		}
	}
	std::cout << "  Found " << count1 << " events in '" << b1name.str() << "' and " << count2 << " in '" << b2name.str() << "'\n";	
	
	TCanvas *can = new TCanvas("can", "Overlay Canvas");
	can->cd();
	hist1->Draw();
	hist2->Draw("SAME");
	can->Update();
	can->WaitPrimitive();
	
	file1->Close();
	if(f2_is_open){ file2->Close(); }
	
	can->Close();
	
	rootapp->Delete();
	
	return 0;
}
