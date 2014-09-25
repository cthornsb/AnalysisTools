// Stitcher.cpp
// C. Thornsberry
// Jul. 30th, 2014
// Stitch multiple root files together into one output file (does not work with newest version of scan code)
// SYNTAX: ./Stitcher {prefix} {treename} {file1} {file2} ...

#include "Loader.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

// Convert arbitrary input to string
template <typename T>
std::string to_str(T input){
	std::stringstream output;
	output << input;
	return output.str();
}

// Open an output root file without overwriting files with the same name
TFile *OpenOutputFile(std::string prefix){
	unsigned int current_file = 1;
	TFile *output = new TFile((prefix + ".root").c_str(), "CREATE");
	while(!output->IsOpen()){
		current_file++;
		output->Close();
		output = new TFile((prefix + to_str(current_file) + ".root").c_str(), "CREATE");
	}
	if(current_file == 1){ std::cout << " Opened output file '" << prefix << ".root'\n"; }
	else{ std::cout << " Opened output file '" << prefix << to_str(current_file) << ".root'\n"; }
	return output;
}

// For compilation
int main(int argc, char* argv[]){
	if(argc < 4){
		std::cout << " Error! Invalid number of arguments. Expected 3, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./Stitcher {prefix} {treename} {file1} {file2} ...\n";
		return 1;
	}

	std::string prefix = to_str(argv[1]);
	std::string treename = to_str(argv[2]);
	
	std::vector<std::string> filenames;
	for(unsigned int i = 3; i < argc; i++){
		filenames.push_back(to_str(argv[i]));
	}

	PixieLoader loader;

	TFile *outFile = OpenOutputFile(prefix + "stitcher");
	TTree *outTree = new TTree(treename.c_str(), "Stitched TTree");
	loader.BeQuiet();

	unsigned long count = 1;
	std::cout << " Working, please wait\n";
	for(std::vector<std::string>::iterator iter = filenames.begin(); iter != filenames.end(); iter++){
		std::cout << "  Loading file '" << prefix << *iter << "'... ";
		if(loader.LoadFile((prefix + *iter).c_str(), treename.c_str()) == 0){
			std::cout << "fail\n";
			continue;
		}
		std::cout << "pass (" << loader.GetEntries() << " entries)\n";
		if(!loader.IsSet()){ loader.SetTree(outTree); }
		for(unsigned int i = 0; i < loader.GetEntries(); i++){
			if(count % 1000000 == 0){ std::cout << "  Entry no. " << count << std::endl; }
			loader.FillEntry(i);
			count++;
		}
	}

	std::cout << " Wrote " << count-1 << " entries to output file\n";

	outFile->cd();
	outTree->Write();
	outFile->Close();
	delete outFile;

	return 0;
}
