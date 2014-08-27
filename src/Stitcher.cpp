// Stitcher.cpp
// C. Thornsberry
// Jul. 30th, 2014
// Stitch multiple root files together into one output file (does not work with newest version of scan code)
// SYNTAX: ./Stitcher

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include <iostream>
#include <sstream>
#include <string>

// Convert arbitrary input to string
template <typename T>
std::string to_str(T input){
	std::stringstream output;
	output << input;
	return output.str();
}

// For CINT
int Stitcher(){ return main(); }

// For compilation
int main(){
	std::string prefix = "/home/cory/MikeF/run26";
	const char* tree_name = "Pixie16";
	int *waveform = new int[124];

	TFile *file;
	TTree *tree;
	TBranch *branch;
	std::string filename;

	TFile *out_file = new TFile("output.root", "RECREATE");
	TTree *out_tree = new TTree("Pixie16", "Tree containing liquid waveforms");
	out_tree->Branch("liquid_wave[124]", waveform, "liquid_wave[124]/I");

	for(unsigned short i = 0; i <= 6; i++){
		filename = prefix;
		if(i == 0){ filename += ".root"; }
		else{ filename += "-" + to_str(i) + ".root"; }

		std::cout << " Working in file " << filename << std::endl;
		file = new TFile(filename.c_str(), "READ");
		tree = (TTree*)file->Get(tree_name);

		if(!tree){ 
			std::cout << "  Failed to load tree " << tree_name << "!\n";
			file->Close();
			delete file;
			continue;
		}

		tree->SetMakeClass(1);
		tree->SetBranchAddress("liquid_wave[124]", waveform, &branch);
		std::cout << "  Stitching " << tree->GetEntries() << " entries... ";

		for(unsigned int j = 0; j < tree->GetEntries(); j++){
			tree->GetEntry(j);
			out_tree->Fill();
		}

		std::cout << "done\n";
		tree->ResetBranchAddresses();
		file->Close();
		delete file;
	}

	std::cout << " Wrote " << out_tree->GetEntries() << " entries to output file\n";

	out_file->cd();
	out_tree->Write();
	out_file->Close();
	delete file;
	delete out_file;
	delete[] waveform;

	return 0;
}
