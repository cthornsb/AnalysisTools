// Loader.cpp
// C. Thornsberry
// Aug. 28th, 2014
// Class to load, read, and copy Pixie16 data

#include "Loader.h"

std::string GetBranchName(std::string classname_){
	std::string output = "";
	std::size_t index = classname_.find("Structure");
	if(index == std::string::npos){ // Not structure
		index = classname_.find("Waveform");
		if(index != std::string::npos){ // Waveform
			for(unsigned short i = 0; i < index; i++){ output += classname_[i]; }
		}
	}
	else{ // Structure
		for(unsigned short i = 0; i < index; i++){ output += classname_[i]; }
		output += "Wave";
	}
	return output;
}

bool ReadLinkDef(std::vector<std::string>& bnames_, const char* fname="./dict/LinkDef.h"){
	std::ifstream ldef_file(fname);
	if(!ldef_file.good()){ return false; }
	
	char line[128];
	std::string temp1, temp2;
	std::size_t index;
	
	while(true){
		ldef_file.getline(line, 128);
		if(ldef_file.eof()){ break; }
		
		temp1 = "";
		for(unsigned  int i = 0; i < 128; i++){
			if(line[i] == '\0'){ break; }
			temp1.push_back(line[i]);
		}
		
		index = temp1.find("C++ class");
		if(index != std::string::npos){ // This line has a class def
			temp2 = "";
			for(unsigned short i = index + 10; i < temp1.size(); i++){
				if(temp1[i] == '+' || temp1[i] == ';'){ break; }
				temp2.push_back(temp1[i]);
			}
			bnames_.push_back(GetBranchName(temp2));
		}
	}
	ldef_file.close();
	return true;
}

void PixieLoader::_initialize(){
	/*if(!ReadLinkDef(bnames)){ 
		std::cout << " PixieLoader: Failed to load LinkDef.h\n"; 
		std::cout << " PixieLoader: Falling back to failsafe\n";
		bnames.push_back("Vandle");
	}*/
	file = NULL;
	tree = NULL;
	outTree = NULL;
	for(unsigned short i = 0; i < 8; i++){ 
		branches[i] = NULL; 
		switches[i] = false;
		outSwitches[i] = false;
		branch_check[i] = false;
	}
	file_open = false;
	bad_file = false;
	check_branches = false;
	silence = false;
	is_set = false;
}

unsigned short PixieLoader::_set_addresses(){
	unsigned short count = 0;	
	tree->SetBranchAddress("Trigger", &trig_structure, &branches[0]);
	tree->SetBranchAddress("TriggerWave", &trig_waveform, &branches[1]);
	tree->SetBranchAddress("Logic", &run_structure, &branches[2]);
	tree->SetBranchAddress("LogicWave", &run_waveform, &branches[3]);
	tree->SetBranchAddress("Liquid", &liq_structure, &branches[4]);
	tree->SetBranchAddress("LiquidWave", &liq_waveform, &branches[5]);
	tree->SetBranchAddress("Vandle", &van_structure, &branches[6]);
	tree->SetBranchAddress("VandleWave", &van_waveform, &branches[7]);
	for(unsigned short i = 0; i < 8; i++){ 
		if(branches[i]){ 
			if(check_branches && !branch_check[i]){ // File has different data structure than expected
				if(!silence){ 
					std::cout << " Warning! Encountered unexpected data structure in input file\n";
					std::cout << " Warning! GetEntry has been disabled until this file is unloaded\n";
				}
				bad_file = true;
				return 0;
			}
			switches[i] = true; 
			outSwitches[i] = true;
			count++;
		}
		else{ 
			switches[i] = false; 
			outSwitches[i] = false;
		}
	}
	
	if(!check_branches){ 
		// First time loading a file. Record which branches are active
		// This will prevent loading files with different branches later
		for(unsigned short i = 0; i < 8; i++){
			if(switches[i]){ branch_check[i] = true; }
			else{ branch_check[i] = false; }
		}
		check_branches = true; // Check branches in input file from now on
	}
	
	return count;
}

PixieLoader::PixieLoader(const char* fname, const char* tname){
	_initialize();
	unsigned short num_good_branches = LoadFile(fname, tname);
	if(!silence && num_good_branches > 0){ std::cout << " PixieLoader: Successfully loaded " << num_good_branches << " branches\n"; }
}

PixieLoader::~PixieLoader(){
	if(file_open){ // Close the file if one is open
		file->Close(); 
		delete file;
	}
}

unsigned short PixieLoader::LoadFile(const char* fname, const char* tname){
	if(file_open){ // Close the file if one is open
		file->Close();
		file_open = false;
		bad_file = false;
		delete file;
	}
		
	// Load the file and the tree
	file = new TFile(fname, "READ"); // Open read only
	if(!file->IsOpen()){
		if(!silence){ std::cout << " Warning! Failed to load input file '" << fname << "'\n"; }
		file_open = false;
		file->Close();
		delete file;
		return 0;
	}
	tree = (TTree*)file->Get(tname);
	if(!tree){
		if(!silence){ std::cout << " Warning! Failed to load input tree '" << tname << "'\n"; }
		file_open = false;
		file->Close();
		delete file;
		return 0;
	}
	file_open = true;
	
	// Load the branches
	return _set_addresses();
}

bool PixieLoader::GetEntry(long long entry){
	if(bad_file){ return false; }
	int return_code = tree->GetEntry(entry);
	if(return_code == 0 || return_code == -1){ return false; }
	return true;
}

bool PixieLoader::Fill(bool do_fill){
	if(bad_file || !outTree){ return false; }
	
	// Actually do the filling
	if(outSwitches[0]){ out_trig_structure = trig_structure; }
	if(outSwitches[1]){ out_trig_waveform = trig_waveform; }
	if(outSwitches[2]){ out_run_structure = run_structure; }
	if(outSwitches[3]){ out_run_waveform = run_waveform; }
	if(outSwitches[4]){ out_liq_structure = liq_structure; }
	if(outSwitches[5]){ out_liq_waveform = liq_waveform; }
	if(outSwitches[6]){ out_van_structure = van_structure; }
	if(outSwitches[7]){ out_van_waveform = van_waveform; }
	
	if(do_fill){
		int return_code = outTree->Fill();
		if(return_code == 0 || return_code == -1){ return false; }
	}
	return true;
}

bool PixieLoader::FillEntry(long long entry, bool do_fill){
	if(bad_file){ return false; }
	if(GetEntry(entry)){ return Fill(do_fill); }
	return false;
}

bool PixieLoader::SetTree(TTree *output_tree){
	if(bad_file || !output_tree){ return false; }
	outTree = output_tree;
	if(outSwitches[0]){ outTree->Branch("Trigger", &out_trig_structure); }
	if(outSwitches[1]){ outTree->Branch("TriggerWave", &out_trig_waveform); }
	if(outSwitches[2]){ outTree->Branch("Logic", &out_run_structure); }
	if(outSwitches[3]){ outTree->Branch("LogicWave", &out_run_waveform); }
	if(outSwitches[4]){ outTree->Branch("Liquid", &out_liq_structure); }
	if(outSwitches[5]){ outTree->Branch("LiquidWave", &out_liq_waveform); }
	if(outSwitches[6]){ outTree->Branch("Vandle", &out_van_structure); }
	if(outSwitches[7]){ outTree->Branch("VandleWave", &out_van_waveform); }
	if(!is_set){ is_set = true; }
	return true;
}

void PixieLoader::SetOutputBranches(std::string input){
	if(input.size() < 8){
		if(!silence){ std::cout << " PixieLoader: Failed to set output branches\n"; }
		return;
	}
	unsigned short count = 0;
	for(unsigned short i = 0; i < 8; i++){
		if(switches[i] && input[i] == '1'){ 
			outSwitches[i] = true; 
			count++;
		}
		else{ outSwitches[i] = false; }
	}
	if(!silence){ std::cout << " PixieLoader: Activated " << count << " output branches\n"; }
}

void PixieLoader::Zero(){
	out_trig_structure.Zero();
	out_trig_waveform.Zero();
	out_run_structure.Zero();
	out_run_waveform.Zero();
	out_liq_structure.Zero();
	out_liq_waveform.Zero();
	out_van_structure.Zero();
	out_van_waveform.Zero();
}

void PixieLoader::Test(){
	std::cout << " PixieLoader: TEST...\n";
	std::cout << "  file_open = " << file_open << ", bad_file = " << bad_file << ", check_branches = " << check_branches << ", silence = " << silence << std::endl;
	
	std::cout << "  Input file test... "; if(file){ std::cout << "pass\n"; }else{ std::cout << "fail\n"; }
	std::cout << "  Input tree test... "; if(tree){ std::cout << "pass\n"; }else{ std::cout << "fail\n"; }
	std::cout << "  Output tree test... "; if(outTree){ std::cout << "pass\n"; }else{ std::cout << "fail\n"; }

	std::cout << "  Input branch switches = ";
	for(unsigned short i = 0; i < 8; i++){
		if(switches[i]){ std::cout << "1"; }else{ std::cout << "0"; }
	}
	std::cout << "\n  Output branch switches = ";
	for(unsigned short i = 0; i < 8; i++){
		if(outSwitches[i]){ std::cout << "1"; }else{ std::cout << "0"; }
	}

	/*if(switches[0]){ std::cout << "  Trigger: energy.size()=" << trig_structure->trigger_energy.size() << " mult=" << trig_structure->trigger_mult << std::endl; }
	if(switches[2]){ std::cout << "  Runtime: energy.size()=" << run_structure->rtime_energy.size() << " mult=" << run_structure->rtime_mult << std::endl; }
	if(switches[4]){ std::cout << "  Liquid: TOF.size()=" << liq_structure->liquid_TOF.size() << " loc.size()=" << liq_structure->liquid_loc.size() << " mult=" << liq_structure->liquid_mult << std::endl; }
	if(switches[6]){ std::cout << "  Vandle: TOF.size()=" << van_structure->vandle_TOF.size() << " loc.size()=" << van_structure->vandle_loc.size() << " mult=" << van_structure->vandle_mult << std::endl; }
	if(switches[1]){ std::cout << "  TriggerWave: wave.size()=" << trig_waveform->trigger_wave.size() << std::endl; }
	if(switches[3]){ std::cout << "  RuntimeWave: wave.size()=" << run_waveform->rtime_wave.size() << std::endl; }
	if(switches[5]){ std::cout << "  LiquidWave: wave.size()=" << liq_waveform->liquid_wave.size() << std::endl; }
	if(switches[7]){ std::cout << "  VandleWave: left.size()=" << van_waveform->left_wave.size() << " right.size()=" << van_waveform->right_wave.size() << std::endl; }*/
	std::cout << "\n\n";
}
