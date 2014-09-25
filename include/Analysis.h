#ifndef PULSEANALYSIS_H
#define PULSEANALYSIS_H

#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"

struct peak{
	unsigned int left, max, min, cfd;
	double max_value, min_value;
	
	peak(){
		left = 0; max = 0; min = 0; cfd = 0; 
		max_value = -9999.0; min_value = 9999.0;
	}
	
	peak(unsigned int left_, unsigned int max_, unsigned int min_, unsigned int cfd_, double max_value_, double min_value_){
		left = left_; max = max_; min = min_; cfd = cfd_; 
		max_value = max_value_; min_value = min_value_;
	}
	
	std::string print(){
		std::stringstream output;
		output << "left = " << left << ", cfd = " << cfd;
		output << ", max_bin = " << max << ", max_value = " << max_value;
		output << ", min_bin = " << min << ", min_value = " << min_value; 
		return output.str();
	}
};

class PulseAnalysis {
private:
	std::string version; 
	std::string revision;

	unsigned int pulse_size;
	unsigned int full_size;
	double *pulse;

	double deltaT, tempV;
	double PSD, integralS, integralL, sumL, sumS;	
	double linear[4], pulse_temp[2000], summing[2000], MovingAverage[2000];	
	double w,x,y,z;
	int i,j,k,l,m,n;	
	int return_code	;
	bool positive, negative;
	bool debug, can_graph;
	
	peak first, second;
	TCanvas *can;
	TGraph *graph;
	double *x_graph;
	const double *x_val, *y_val;

public:
	PulseAnalysis();
	~PulseAnalysis();
	void GetVersion();	
	int PSD_Integration (int, double&, double&);
	int Baseline_restore (double*, int, int);
	int PSD_Zerocross (int, int, double&);	
	int Parameters (int, double&, double&, double&, double&, double&);	
	int Time_Pickoff (int, int, int, int, double&);
	int Derivative(int);
	int	Integral(double*);
	int PeakFinder (int, int, int, int&, int&);
	int OptimizePSD (int, int, int, int, double*, double*);
	int Half_Integral(double&);
	int Smooth(int, int, double);
	int HPGe(double&);
	
	void SetPulseSize(unsigned int);
	
	void ToggleDebug(){ 
		if(debug){ debug = false; } 
		else{ debug = true; }
	}
	
	unsigned int PreProcess(std::vector<int>&);
	unsigned int Process(unsigned int, unsigned int);
};

#endif
