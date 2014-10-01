#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"

//-------------------------------------------------------------------------------------//

std::string filename, treename;
TFile *file;
TTree *tree;
std::vector<TH1F*> hists1d;
std::vector<TH2F*> hists2d;

unsigned int count1d;
unsigned int count2d;

TCanvas *can;

//-------------------------------------------------------------------------------------//

void setFile(std::string filename_){ filename = filename_; }

void setTree(std::string treename_="Pixie16"){ treename = treename_; }

bool load(std::string filename_, std::string treename_="Pixie16"){
	filename = filename_; treename = treename_;
	return load();
}

bool load(){
	if(file && file->IsOpen){ file->Close(); }
	
	// Load the file
	file = new TFile(filename.c_str(), "READ");
	if(!file->IsOpen()){ 
		std::cout << " Warning! Failed to load input file '" << filename << "'\n";
		return false; 
	}
	
	// Load the tree
	tree = (TTree*)file->Get(treename.c_str());
	if(!tree){
		std::cout << " Warning! Failed to load the tree named '" << treename << "'\n";
		return false; 
	}
	
	return true;
}

void ls(){
	gDirectory->ls();
}

TBrowser *browse(){ return new TBrowser; }

//-------------------------------------------------------------------------------------//

void push1D(){
	std::stringstream stream;
	if(count1d < 10){ stream << "hist0"; }
	else{ stream << "hist"; }
	stream << count1d << "_1d";
	
	TH1F *hist = (TH1F*)((TH1F*)gPad->GetPrimitive("")->Clone(stream.str().c_str()));
	hists1d.push_back(hist);
	count1d++;
}

long long plot1D(std::string name, int bins, float low, float high, std::string gate="", std::string draw="", bool push=true){
	if(!file || !file->IsOpen()){ return -1; }
	std::stringstream stream;
	stream << name << ">>(" << bins << "," << low << "," << high << ")";
	
	std::cout << " tree->Draw(\"" << stream.str() << "\", \"" << gate << "\", \"" << draw << "\")\n";
	
	long long return_code = tree->Draw(stream.str().c_str(), gate.c_str(), draw.c_str());
	if(push){ push1D(); }
	return return_code;
}

void list1D(){
	/*unsigned int count = 0;
	for(std::vector<TH1F*>::iterator iter = hists1d.begin(); iter != hists1d.end(); iter++){
		std::cout << " " << count << ": " << iter->GetName() << std::endl;
		count++;
	}*/
	for(unsigned int i = 0; i < hists1d.size(); i++){
		std::cout << " " << i << ": " << hists1d[i]->GetName() << std::endl;
	}
}

void draw1D(unsigned int index, std::string draw=""){
	if(index < hists1d.size()){ hists1d[index]->Draw(draw.c_str()); }
	else{ std::cout << " 1D histogram index is out of range\n"; }
}

TH1F *get1D(unsigned int index){
	if(index < hists1d.size()){ return hists1d[index]; }
	else{ std::cout << " 1D histogram index is out of range\n"; }
	return NULL;
}

bool save1D(TFile *out_file){
	if(!out_file){ return false; }
	out_file->cd();
	/*for(std::vector<TH1F*>::iterator iter = hists1d.begin(); iter != hists1d.end(); iter++){
		iter->Write();
	}*/
	for(unsigned int i = 0; i < hists1d.size(); i++){
		hists1d[i]->Write();
	}
	return true;
}

bool save1D(std::string outname){
	TFile *out_file = new TFile(outname.c_str(), "CREATE");
	if(!out_file->IsOpen()){
		std::cout << " Failed to open output file '" << outname << "'\n";
		return false;
	}
	bool output = save1D(out_file);
	out_file->Close();
	return output;
}

//-------------------------------------------------------------------------------------//

void push2D(){
	std::stringstream stream;
	if(count2d < 10){ stream << "hist0"; }
	else{ stream << "hist"; }
	stream << count2d << "_2d";
	
	TH2F *hist = (TH2F*)((TH2F*)gPad->GetPrimitive("")->Clone(stream.str().c_str()));
	hists2d.push_back(hist);
	count2d++;
}

long long plot2D(std::string xname, int xbins, float xlow, float xhigh, std::string yname, int ybins, float ylow, float yhigh, std::string gate="", std::string draw="", bool push=true){
	if(!file || !file->IsOpen()){ return -1; }
	std::stringstream stream;
	stream << yname << ":" << xname << ">>(";
	stream << xbins << "," << xlow << "," << xhigh << ",";
	stream << ybins << "," << ylow << "," << yhigh << ")";
	
	std::cout << " tree->Draw(\"" << stream.str() << "\", \"" << gate << "\", \"" << draw << "\")\n";
	
	long long return_code = tree->Draw(stream.str().c_str(), gate.c_str(), draw.c_str());
	if(push){ push2D(); }
	return return_code;
}

void list2D(){
	/*unsigned int count = 0;
	for(std::vector<TH2F*>::iterator iter = hists2d.begin(); iter != hists2d.end(); iter++){
		std::cout << " " << count << ": " << iter->GetName() << std::endl;
		count++;
	}*/
	for(unsigned int i = 0; i < hists2d.size(); i++){
		std::cout << " " << i << ": " << hists2d[i]->GetName() << std::endl;
	}
}

void draw2D(unsigned int index, std::string draw=""){
	if(index < hists2d.size()){ hists2d[index]->Draw(draw.c_str()); }
	else{ std::cout << " 2D histogram index is out of range\n"; }
}

TH2F *get2D(unsigned int index){
	if(index < hists2d.size()){ return hists2d[index]; }
	else{ std::cout << " 2D histogram index is out of range\n"; }
	return NULL;
}

bool save2D(TFile *out_file){
	if(!out_file){ return false; }
	out_file->cd();
	/*for(std::vector<TH2F*>::iterator iter = hists2d.begin(); iter != hists2d.end(); iter++){
		iter->Write();
	}*/
	for(unsigned int i = 0; i < hists2d.size(); i++){
		hists2d[i]->Write();
	}
	return true;
}

bool save2D(std::string outname){
	TFile *out_file = new TFile(outname.c_str(), "CREATE");
	if(!out_file->IsOpen()){
		std::cout << " Failed to open output file '" << outname << "'\n";
		return false;
	}
	bool output = save2D(out_file);
	out_file->Close();
	return output;
}

//-------------------------------------------------------------------------------------//

void list(){
	list1D();
	list2D();
}

bool save(TFile *out_file){
	return(save1D(out_file) && save2D(out_file));
}

bool save(std::string outname){
	TFile *out_file = new TFile(outname.c_str(), "CREATE");
	if(!out_file->IsOpen()){
		std::cout << " Failed to open output file '" << outname << "'\n";
		return false;
	}
	bool output = save1D(out_file) && save2D(out_file);
	out_file->Close();
	return output;
}

//-------------------------------------------------------------------------------------//

int RootHelper(){
	std::cout << " Welcome to RootHelper v. 1.0 by C. Thornsberry\n";
	std::cout << "  For a help menu, type 'help()' or 'help(command)'\n";
	filename = "";
	treename = "";
	file = NULL;
	tree = NULL;
	can = new TCanvas("RootHelper v. 1.0");
	
	count1d = 0;
	count2d = 0;
	
	return 0;
}
