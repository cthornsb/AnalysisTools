// Loader.h
// C. Thornsberry
// Aug. 28th, 2014
// Class to load, read, and copy Pixie16 data

#ifndef LOADER_H
#define LOADER_H

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "Structures.h"

// Get the branch name of a pixie branch from the name of its structure class
std::string GetBranchName(std::string classname_);

// Load class names from a root LinkDef file
bool LoadLinkDef(std::vector<std::string> &bnames_, const char* fname="../dict/LinkDef.h");

class PixieLoader{
  private:
	bool file_open;
	bool switches[8], outSwitches[8];
	TFile *file;
	TTree *tree, *outTree;
	TBranch *branches[8];
	std::vector<std::string> bnames;
	
	TriggerStructure out_trig_structure, *trig_structure = NULL;
	TriggerWaveform out_trig_waveform, *trig_waveform = NULL;
	RuntimeStructure out_run_structure, *run_structure = NULL;
	RuntimeWaveform out_run_waveform, *run_waveform = NULL;
	LiquidStructure out_liq_structure, *liq_structure = NULL;
	LiquidWaveform out_liq_waveform, *liq_waveform = NULL;
	VandleStructure out_van_structure, *van_structure = NULL;
	VandleWaveform out_van_waveform, *van_waveform = NULL;
	
	void _initialize();
	
	unsigned short _set_addresses();
  
  public:
	PixieLoader(){ _initialize(); }
	PixieLoader(const char*, const char*);
	~PixieLoader();	
	
	unsigned short LoadFile(const char*, const char*);
	
	TriggerStructure* GetTrigger(){ return trig_structure; }
	TriggerStructure* GetTriggerOut(){ return &out_trig_structure; }
	
	TriggerWaveform* GetTriggerWave(){ return trig_waveform; }
	TriggerWaveform* GetTriggerWaveOut(){ return &out_trig_waveform; }
	
	RuntimeStructure* GetRuntime(){ return run_structure; }
	RuntimeStructure* GetRuntimeOut(){ return &out_run_structure; }
	
	RuntimeWaveform* GetRuntimeWave(){ return run_waveform; }
	RuntimeWaveform* GetRuntimeWaveOut(){ return &out_run_waveform; }
	
	LiquidStructure* GetLiquid(){ return liq_structure; }
	LiquidStructure* GetLiquidOut(){ return &out_liq_structure; }
	
	LiquidWaveform* GetLiquidWave(){ return liq_waveform; }
	LiquidWaveform* GetLiquidWaveOut(){ return &out_liq_waveform; }
	
	VandleStructure* GetVandle(){ return van_structure; }
	VandleStructure* GetVandleOut(){ return &out_van_structure; }
	
	VandleWaveform* GetVandleWave(){ return van_waveform; }
	VandleWaveform* GetVandleWaveOut(){ return &out_van_waveform; }
	
	// Get number of entries in input tree
	long long GetEntries(){ return tree->GetEntries(); }
	
	// Get entry from input tree
	bool GetEntry(long long entry);
	
	// Fill output tree
	bool Fill();
	
	// Get entry and fill output tree
	bool FillEntry(long long entry);
	
	// Set branches on a new tree
	bool SetTree(TTree*);
	
	// Turn output branches on/off
	void SetOutputBranches(std::string);
	
	void Test();
};

#endif
