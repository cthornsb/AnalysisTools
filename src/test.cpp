#include "Loader.h"

#include "TH2D.h"
#include "TCanvas.h"
#include "TApplication.h"

int main(){
	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);

	PixieLoader loader("../Analysis/run27.root", "Pixie16");
	
	TFile *test_file = new TFile("test.root", "RECREATE");
	TTree *test_tree = new TTree("test_tree", "test_tree");
	loader.SetOutputBranches("11000000");
	loader.SetTree(test_tree);
	
	std::cout << " Plotting " << loader.GetEntries() << " entries\n";
	for(unsigned int i = 0; i < loader.GetEntries(); i++){
		loader.FillEntry(i);
	}
	
	test_file->cd();
	test_tree->Write();
	test_file->Close();
		
	rootapp->Delete();
		
	return 0;
}
