// PulseAnalyzer.cpp
// C. Thornsberry
// Aug. 25th, 2014
// Load detector pulses and analyze them using various integration methods
// SYNTAX: Viewer {root_filename} {root_treename} {root_branch} {method#} [debug]

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCutG.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <cmath>

const double PIXIE_TIME_RES = 4.0; // In ns
const double ROOT2PI = std::sqrt(2*3.14159);

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

// 1D function to use for pulse fitting
// x[0] = time t in ns
// parameters: 4
//  par[0] = alpha (normalization of pulse 1)
//  par[1] = phi (phase of pulse 1 in ns)
//  par[2] = beta (decay parameter of the 1st pulse exponential in ns)
//  par[3] = gamma (width of the inverted square gaussian of the 1st pulse in ns^4)
double fit_function(double *x, double *par){
	double arg = x[0] - par[1];
	if(arg >= 0.0){ return par[0]*std::exp(-arg/par[2])*(1 - std::exp(-arg*arg*arg*arg/par[3])); }
	return 0.0;
}

// 1D pulse function with delayed double gaussian (bump and valley)
// x[0] = time t in ns
// parameters: 10
//  par[0] = alpha (normalization of pulse 1)
//  par[1] = phi (phase of pulse 1 in ns)
//  par[2] = beta (decay parameter of the 1st pulse exponential in ns)
//  par[3] = gamma (width of the inverted square gaussian of the 1st pulse in ns^4)
//  par[4] = amplitude of delayed "bump"
//  par[5] = mu of "bump"
//  par[6] = sigma of "bump"
//  par[7] = amplitude of "valley"
//  par[8] = mu of "valley" (minimum)
//  par[9] = sigma of "valley"
double fit_function2(double *x, double *par){
	if(x[0] < par[1]){ return 0.0; }
	double arg = x[0] - par[1];
	double sum = par[0]*std::exp(-arg/par[2])*(1 - std::exp(-arg*arg*arg*arg/par[3])); // Pulse
	sum += (par[4]/(par[6]*ROOT2PI))*std::exp(-(x[0]-par[5])*(x[0]-par[5])/(2*par[6]*par[6])); // Bump
	sum += (par[7]/(par[9]*ROOT2PI))*std::exp(-(x[0]-par[8])*(x[0]-par[8])/(2*par[9]*par[9])); // Valley
	return sum;
}

double variance(double *data, unsigned int start, unsigned int stop){
	double mean = 0.0;
	double output = 0.0;
	for(unsigned int i = start; i < stop; i++){
		mean += data[i];
	}
	mean = mean/(stop-start);
	for(unsigned int i = start; i < stop; i++){
		output += (data[i] - mean)*(data[i] - mean);
	}
	return output/(stop-start);
}

double chi_2(double *data, unsigned int start, unsigned int stop, double *pars){
	double var = variance(data, start, stop);
	double output = 0.0, temp;
	double x_value[1] = {0.0};
	for(unsigned int i = start; i < stop; i++){
		x_value[0] = i;
		temp = (data[i] - fit_function(x_value, pars));
		output += temp*temp/var;
	}
	return output;
}

double interpolate(double x, double x1, double y1, double x2, double y2){
	if(x == x1){ return y1; }
	if(x == x2){ return y2; }
	return (y1 + ((y2-y1)/(x2-x1))*(x-x1)); 
}

unsigned int preprocess(double *input, unsigned int start, unsigned int stop, unsigned int size, peak *first, peak *second, bool debug=false){
	// Do baseline correction. Use the range to the left of the pulse.
	// The right side of the pulse may have other pulses
	double base_sum = 0.0;
	unsigned int base_max = start + (unsigned int)(0.1*size); // First 10% of range
	for(unsigned int i = start; i < base_max; i++){
		base_sum += input[i];
	}

	double base_line = base_sum/(base_max-start);
	for(unsigned int i = start; i < stop; i++){
		input[i] = input[i] - base_line;
	}

	// Find the global maximum
	double maximum = -9999.0;
	unsigned int maximum_bin = 0;
	for(unsigned int i = start; i < stop; i++){
		if(input[i] > maximum){ 
			maximum = input[i];
			maximum_bin = i;
		}
	}

	// Peak-finder (there could be more than 1)
	unsigned int leading_edge, peak_bin;
	double back_sub_left = 0.0;
	double back_sub_right = 0.0;
	double slope;
	bool up_slope = false;
	bool down_slope = false;
	bool found_peak = false;
	unsigned short num_peaks = 0;
	for(unsigned int i = start; i < stop-1; i++){
		// Subtract background (10% of maximum)
		back_sub_left = input[i] - 0.1*maximum;
		back_sub_right = input[i+1] - 0.1*maximum;
		if(back_sub_left < 0.0){ back_sub_left = 0.0; }
		if(back_sub_right < 0.0){ back_sub_right = 0.0; }
		
		// Calculate the slope
		slope = back_sub_right - back_sub_left;
		if(!up_slope && slope > 0.0){ // Moving to up_slope marks leading edge of pulse
			up_slope = true; 
			down_slope = false;
			leading_edge = i;
		}
		else if(up_slope && slope < 0.0){ // Moving from up_slope to down_slope marks a pulse peak
			down_slope = true; 
			up_slope = false;
			found_peak = true;
			num_peaks++;
			peak_bin = i;
		}
		else if(down_slope && slope > 0.0){ // Moving from down_slope to non-negative slope marks trailing edge (reset peak-finder)
			up_slope = false;
			down_slope = false;
		}
		
		if(found_peak){
			found_peak = false;
			if(num_peaks == 1){
				first->left = leading_edge;
				first->max = peak_bin;
				first->max_value = input[peak_bin];
			}
			else if(num_peaks == 2){
				second->left = leading_edge;
				second->max = peak_bin;
				second->max_value = input[peak_bin];
				break; // That's enough!
			}
		}
	}

	// Find the local minimum
	if(num_peaks == 1){ // Single pulse, minimum must occur to the right of the pulse
		for(unsigned int i = first->max; i < stop; i++){
			if(input[i] < first->min_value){
				first->min_value = input[i];
				first->min = i;
			}
		}
	}
	else{ // Multiple pulses, minimum must occur between pulse peaks
		for(unsigned int i = first->max; i < second->max; i++){
			if(input[i] < first->min_value){
				first->min_value = input[i];
				first->min = i;
			}
		}
	}

	// Find 80% of the peak maximum
	for(unsigned int i = first->left; i <= first->max; i++){
		if(input[i] >= 0.8*first->max_value){
			if(i > 0){ first->cfd = i-1; }
			else{ first->cfd = i; }
			break;
		}
	}
	if(first->cfd == first->max){ first->cfd--; }
	
	// Print debug information
	if(debug){ 
		std::cout << "Global: baseline = " << base_line << ", maximum_bin = " << maximum_bin;
		std::cout << ", maximum = " << maximum << ", num_peaks = " << num_peaks << std::endl;
		std::cout << "Peak: " << first->print() << std::endl;
	}
	
	return num_peaks;
}

// Integrate a pulse and return the short and long integrals
unsigned int integrate(std::vector<int> &arr, unsigned int pulse_size, std::vector<double> &s_int, std::vector<double> &l_int, 
					   unsigned short method, bool debug, TCanvas *canvas1, TCanvas *canvas2, bool use_in, std::ifstream &infile, std::ofstream &outfile){
	if(debug && (!canvas1 || !canvas2)){ debug = false; }
	unsigned int size = arr.size();
	s_int.clear(); l_int.clear();
	if(size == 0){ return 0; }
	if(size % pulse_size != 0){ return 0; }
	unsigned int num_pulses = size/pulse_size; // Expects size to be divisible by pulse_size	
	
	// Copy the pulse into a new array so the original pulse is unmodified
	double *darr = new double[size];
	unsigned int count = 0;
	for(std::vector<int>::iterator iter = arr.begin(); iter != arr.end(); iter++){
		if(count >= size){ break; }
		darr[count] = (double)(*iter);
		count++;
	}
	
	// Variables for TGraph
	double *x1 = new double[pulse_size];
	double *x2 = new double[pulse_size];
	const double *x_val1 = new double[pulse_size];
	const double *y_val1 = new double[pulse_size];
	const double *x_val2 = new double[pulse_size];
	const double *y_val2 = new double[pulse_size];
	
	TGraph *graph1 = new TGraph(pulse_size, x_val1, y_val1);
	TGraph *graph2 = new TGraph(pulse_size, x_val2, y_val2);
	for(unsigned int i = 0; i < pulse_size; i++){ x1[i] = i; x2[i] = i; }
		
	if(debug){
		graph1->SetLineColor(2); graph1->SetLineWidth(2);
		graph2->SetLineColor(4); graph2->SetLineWidth(2);
		canvas1->cd();
	}
	
	bool write_param;
	if(!use_in && outfile.good()){ write_param = true; }
	else{ write_param = false; }
	
	if(use_in && !infile.good()){ use_in = false; }
	
	unsigned int peak_count;
	unsigned int start, stop;
	peak first_peak;
	peak second_peak;
	double s, l;
	
	//if(num_pulses != 2){ return 0; }
	//for(unsigned int pulse = 1; pulse < num_pulses; pulse++){
	for(unsigned int pulse = 0; pulse < num_pulses; pulse++){
		start = pulse*pulse_size;
		stop = (pulse+1)*pulse_size;

		peak_count = preprocess(darr, start, stop, pulse_size, &first_peak, &second_peak, debug);

		if(debug || method == 5 || method == 6){
			for(unsigned int i = start; i < stop; i++){ 
				graph1->SetPoint(i-start, x1[i-start], darr[i]); 
			}
		}

		if(peak_count == 0){ // Failed to find a peak
			if(debug){ std::cout << "No peaks found, skipping entry\n"; }
			continue;
		}
		if((method >= 4) && peak_count > 1){ // Multiple peaks will throw off the fitting functions
			if(debug){ std::cout << "Multiple peaks, skipping entry\n"; }
			continue;
		}
	
		// Do short and long integrals
		s = 0.0; l = 0.0;			
		if(method == 0){ // Slope inversion
			// To the right of the maximum, we expect a negative slope
			// Stop short integral when the slope goes non-negative
			bool calc_short = true;
			for(unsigned int i = first_peak.cfd; i < first_peak.min; i++){
				// Integrate using trapezoid rule
				if(calc_short){ // Waiting on positive slope to stop calculating short integral
					if(i < first_peak.max || darr[i+1]-darr[i] <= 0.0){ s += (darr[i] + darr[i+1])/2.0; } // Negative slope, still on side of pulse
					else{ calc_short = false; }
				}
				l += (darr[i] + darr[i+1])/2.0;
			}
		} // method == 0
		else if(method == 1){ // Fixed pulse height
			// Stop short integral when we reach 10% of the pulse height
			// Find the 10% of maximum point
			unsigned int ten_percent_bin = 0;
			for(unsigned int i = first_peak.max; i < stop; i++){
				if(darr[i] <= 0.1*first_peak.max_value){ 
					ten_percent_bin = i;
					break;
				}
			}
			
			if(debug){ std::cout << " method = 1, 10p_bin = " << ten_percent_bin << std::endl; }	
			
			for(unsigned int i = first_peak.cfd; i < first_peak.min; i++){
				// Integrate using trapezoid rule with arbitrary bin width of 1
				if(i < ten_percent_bin){ s += (darr[i] + darr[i+1])/2.0; } // Stop calculating short below 10% maximum
				l += (darr[i] + darr[i+1])/2.0;
			}
		} // method == 1
		else if(method == 2){ // Search for decay-side slope variation
			double darr_prime[pulse_size];
			for(unsigned int i = start; i < stop-1; i++){
				darr_prime[i-start] = darr[i+1] - darr[i];
			}
			darr_prime[pulse_size-1] = 0.0;

			// Find the most negative slope, this will be the characteristic slope of the pulse
			double minimum_slope = 9999.0;
			unsigned int minimum_slope_bin = 0, sstop = 0;
			for(unsigned int i = start; i < stop; i++){
				if(darr_prime[i] < minimum_slope){
					minimum_slope = darr_prime[i-start];
					minimum_slope_bin = i-start;
				}
			}

			// Find the bin in which to stop short integration
			for(unsigned int i = minimum_slope_bin; i < stop; i++){
				if(darr_prime[i] >= 0.5*minimum_slope){
					sstop = i;
					break;
				}
			}

			if(debug){
				std::cout << " method = 2, slope_bin = " << minimum_slope_bin << ", min_slope = " << minimum_slope;
				std::cout << ", sstop = " << sstop << ", lstop = " << first_peak.min << std::endl;			
				for(unsigned int i = start; i < stop; i++){ // 1st derivative
					graph2->SetPoint(i-start, x2[i-start], darr_prime[i-start]);
				}
				canvas1->Clear();
				graph2->Draw();
				graph1->Draw("SAME");
				canvas1->Update();
				canvas1->WaitPrimitive();
			}	
			
			// Integrate using trapezoid rule
			for(unsigned int i = first_peak.cfd; i < first_peak.min; i++){
				if(i < sstop){ s += (darr[i] + darr[i+1])/2.0; }
				l += (darr[i] + darr[i+1])/2.0;
			}
			
			if(debug && s == 0){ 
				std::cout << "peak_count = " << peak_count << std::endl;
				std::cout << "slope_bin = " << minimum_slope_bin << ", min_slope = " << minimum_slope;
				std::cout << ", sstop = " << sstop << ", lstop = " << first_peak.min << std::endl;
				std::cout << first_peak.print() << std::endl << std::endl;
			}
		} // method == 2
		else if(method == 3){ // Ratio of "bump" amplitude to pulse amplitude
			s = darr[(first_peak.min % pulse_size) - 3]/first_peak.max_value;
		} // method == 3
		else if(method == 4){ // Use fast fit guessing (S vs. L version) can handle multiple peaks
			double step = (first_peak.min - first_peak.cfd)/100.0; // Integration step
			double x_value[1] = {0.0};

			// Set optimized starting values
			double parameters[4];
			parameters[0] = first_peak.max_value*9.211 + 150.484;			// Normalization of pulse
			parameters[1] = (first_peak.left % pulse_size)*1.087 - 2.359;	// Phase (leading edge of pulse) (ns)
			parameters[2] = 1.7750575;										// Decay constant of exponential (ns)
			parameters[3] = 115.64125;										// Width of inverted square guassian (ns^4)

			if(debug){ 
				std::cout << " method = 4, par[0] = " << parameters[0] << ", par[1] = " << parameters[1] << ", par[2] = ";
				std::cout << parameters[2] << ", par[3] = " << parameters[3] << std::endl;
				
				TF1 *func = new TF1("func",fit_function,0.0,pulse_size,4);
				func->SetLineColor(4); func->SetLineWidth(2);
				func->SetParameters(parameters);

				canvas1->Clear();
				graph1->Draw();
				func->Draw("SAME");
				canvas1->Update();
				canvas1->WaitPrimitive();
				
				func->Delete();			
			}	

			// Find the RMS of the "fit" function
			double temp;
			for(unsigned int j = 0; j <= 100; j++){
				x_value[0] = (first_peak.cfd % pulse_size) + j*step;
				temp = fit_function(x_value, parameters);
				s += temp*temp;
			}
			s = std::sqrt(s / 100.0);
			
			// Find the RMS of the pulse
			for(unsigned int j = first_peak.cfd; j <= first_peak.min; j++){
				l += darr[j]*darr[j];
			}
			l = std::sqrt(l / (first_peak.min - first_peak.cfd));
		} // method == 4
		else if(method == 5){ // Use fast fit guessing (chi^2 version)
			TF1 *func = new TF1("func",fit_function,0.0,pulse_size,4);
			
			// Set optimized starting values
			func->FixParameter(0, first_peak.max_value*9.211 + 150.484); // Normalization of pulse
			func->FixParameter(1, (first_peak.left % pulse_size)*1.087 - 2.359); // Phase (leading edge of pulse) (ns)
			func->FixParameter(2, 1.7750575); // Decay constant of exponential (ns)
			func->FixParameter(3, 115.64125); // Width of inverted square guassian (ns^4)

			if(debug){ 
				func->SetLineColor(4); func->SetLineWidth(2);
				std::cout << " method = 5, par[0] = " << func->GetParameter(0) << ", par[1] = " << func->GetParameter(1) << ", par[2] = ";
				std::cout << func->GetParameter(2) << ", par[3] = " << func->GetParameter(3) << std::endl;
				canvas1->Clear();
				graph1->Draw();
				func->Draw("SAME"); 
				canvas1->Update();
				canvas1->WaitPrimitive();
			}

			TFitResultPtr func_ptr = graph1->Fit(func,"S Q");
			
			double y;
			for(unsigned int j = first_peak.max; j <= first_peak.min; j++){
				y = func->Eval((double)j); 
				s += (darr[j]-y)*(darr[j]-y);
			}
			s = std::sqrt(s / (first_peak.min - first_peak.left));
			l = func_ptr->Chi2();
			
			func->Delete();
		} // method == 5
		else if(method == 6){ // Use root fitting (slow)
			TF1 *func = new TF1("func",fit_function,0.0,pulse_size,4);
			
			double temp;
			if(use_in){
				for(unsigned int j = 0; j < 4; j++){
					infile >> temp;
					func->FixParameter(j, temp);
					if(infile.eof()){ use_in = false; break; }
				}
			}
			else{ // Set optimized starting values
				func->SetParameter(0, first_peak.max_value*9.211 + 150.484); // Normalization of pulse
				func->SetParameter(1, (first_peak.left % pulse_size)*1.087 - 2.359); // Phase (leading edge of pulse) (ns)
				func->SetParameter(2, 1.7750575); // Decay constant of exponential (ns)
				func->SetParameter(3, 115.64125); // Width of inverted square guassian (ns^4)
			}
			TFitResultPtr func_ptr = graph1->Fit(func,"S Q");
			
			if(debug){ 
				func->SetLineColor(4); func->SetLineWidth(2);
				std::cout << " method = 6, par[0] = " << func->GetParameter(0) << ", par[1] = " << func->GetParameter(1) << ", par[2] = ";
				std::cout << func->GetParameter(2) << ", par[3] = " << func->GetParameter(3) << std::endl;
				graph1->Draw();
				func->Draw("SAME");
				canvas1->Update();
				canvas1->WaitPrimitive(); 
			}
			if(write_param){
				outfile << func->GetParameter(0) << "\t" << func->GetParameter(1) << "\t";
				outfile << func->GetParameter(2) << "\t" << func->GetParameter(3) << "\n";
			}

			double y;
			for(unsigned int j = first_peak.max; j <= first_peak.min; j++){
				y = func->Eval((double)j); 
				s += (darr[j]-y)*(darr[j]-y);
			}
			s = std::sqrt(s / (first_peak.min - first_peak.left));
			l = func_ptr->Chi2();
			
			func->Delete();
		} // method == 6
		else if(method == 7){
			TF1 *func = new TF1("func",fit_function2,0.0,pulse_size,10);

			double params[10];
			if(use_in){
				for(unsigned int j = 0; j < 4; j++){
					infile >> params[j];
					func->FixParameter(j, params[j]);
					if(infile.eof()){ use_in = false; break; }
				}
			}
			else{ // Set optimized starting values
				params[0] = first_peak.max_value*9.211 + 150.484; // Normalization of pulse
				params[1] = (first_peak.left % pulse_size)*1.087 - 2.359; // Phase (leading edge of pulse) (ns)
				params[2] = 1.7750575; // Decay constant of exponential (ns)
				params[3] = 115.64125; // Width of inverted square guassian (ns^4)
				for(unsigned int j = 0; j < 4; j++){ 
					func->SetParameter(j, params[j]); // Scintillator pulse. Root is able to fit these really well.
				}
			}
			
			params[4] = 2*darr[(first_peak.min % pulse_size) - 3]; // Amplitude of bump
			params[5] = (first_peak.min % pulse_size) - 3; // Mean of bump
			params[6] = 1.0; // Sigma of bump
			
			params[7] = 2*darr[first_peak.min]; // Amplitude of valley
			params[8] = first_peak.min % pulse_size; // Mean of valley
			params[9] = 1.0; // Sigma of valley

			for(unsigned int j = 4; j < 10; j++){ 
				func->FixParameter(j, params[j]); // Double delayed Gaussian. Root tends to screw these up, so fix their values.
			}
			TFitResultPtr func_ptr = graph1->Fit(func,"S Q");
			
			if(debug){ 
				func->SetLineColor(4); func->SetLineWidth(2);
				std::cout << " method = 7\n  par[0] = " << func->GetParameter(0) << ", par[1] = " << func->GetParameter(1) << ", par[2] = " << func->GetParameter(2) << "  par[3] = " << func->GetParameter(3) << std::endl;
				std::cout << "  par[4] = " << func->GetParameter(4) << ", par[5] = " << func->GetParameter(5) << ", par[6] = " << func->GetParameter(6) << std::endl;
				std::cout << "  par[7] = " << func->GetParameter(7) << ", par[8] = " << func->GetParameter(8) << ", par[9] = " << func->GetParameter(9) << std::endl;
				graph1->Draw("AL");
				func->Draw("SAME");
				canvas1->Update();
				canvas1->WaitPrimitive(); 
			}
			if(write_param){
				outfile << func->GetParameter(0) << "\t" << func->GetParameter(1) << "\t";
				outfile << func->GetParameter(2) << "\t" << func->GetParameter(3) << "\n";
			}
			
			/*// Calculate short and long RMS. Short will be the bare pulse. Long will be the pulse with double gaussian.
			double step = (first_peak.min - first_peak.cfd)/100.0; // Integration step
			double x_value1[1] = {0.0};
			double x_value2[1] = {0.0};
			
			// Integrate the pulse function
			double local_start = (double)(first_peak.cfd % pulse_size);
			for(unsigned int j = 0; j < 100; j++){
				x_value1[0] =  local_start + j*step;
				x_value2[0] = local_start + (j+1)*step;
				s += (x_value2[0] - x_value1[0])*(fit_function(x_value1, params) + fit_function(x_value2, params))/2.0;
			}

			// Integrate the pulse and double gaussian function
			for(unsigned int j = 0; j <= 100; j++){
				x_value1[0] = local_start + j*step;
				x_value2[0] = local_start + (j+1)*step;
				l += (x_value2[0] - x_value1[0])*(fit_function2(x_value1, params) + fit_function2(x_value2, params))/2.0;
			}*/
			
			s = func_ptr->Chi2();
			
			func->Delete();
		} // method == 7
		if(debug){ std::cout << " s = " << s << ", l = " << l << std::endl; }
		s_int.push_back(s);
		l_int.push_back(l);
	} // Over num_pulses
	
	graph1->Delete(); 
	graph2->Delete();
	delete[] darr; delete[] x1; delete[] x2; delete[] x_val1; 
	delete[] y_val1; delete[] x_val2; delete[] y_val2;
	return num_pulses;
}

// For compilation
int main(int argc, char* argv[]){
	std::string method_names[8] = {"slope inversion (S vs. L)", "fixed height (S vs. L)", "decay-side slope variation (S vs. L)", 
								   "pulse height ratio", "fast fit guess (RMS)", "fast fit guess (chi^2)", "full root fitting",
								   "full root fitting w/ delayed gaussians"};

	if(argc < 5){
		std::cout << " Error! Invalid number of arguments. Expected 4, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: Viewer {filename} {treename} {branch} {method#} [debug]\n";
		for(unsigned short i = 0; i < 8; i++){
			std::cout << "  Method " << i << ": " << method_names[i] << std::endl;
		}
		return 1;
	}
		
	// Variables for root graphics
	char* dummy[0]; 
	TApplication* rootapp = new TApplication("rootapp",0,dummy);
	gSystem->Load("libTree");
	
	unsigned short method = atol(argv[4]);
	std::cout << " Using analysis method " << method << ": " << method_names[method] << std::endl;
	if(method > 7){
		std::cout << " Encountered undefined method (" << method << "), aborting\n";
		return 1;
	}
	
	bool debug = false;
	bool use_cut = false;
	if(argc > 5){
		if(strcmp(argv[5], "debug") == 0){ 
			debug = true; 
			std::cout << " DEBUGGING...\n";
		}
		else if(strcmp(argv[5], "cut") == 0){ use_cut = true; }
	}

	// Branch variables
	std::vector<int> wave;
	std::vector<double> energy;
	unsigned int mult;
	unsigned int wave_size = 0;

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
	
	std::stringstream branch_name;
	branch_name << argv[3];
	TBranch *b_wave, *b_mult, *b_energy;
	tree->SetBranchAddress((branch_name.str()+"_wave").c_str(), &wave, &b_wave);
	tree->SetBranchAddress((branch_name.str()+"_energy").c_str(), &energy, &b_energy);
	tree->SetBranchAddress((branch_name.str()+"_mult").c_str(), &mult, &b_mult);
	
	if(!b_wave){
		std::cout << " Failed to load the input branch '" << branch_name.str() << "_wave'\n";
		file->Close();
		return 1;
	}
	if(!b_energy){
		std::cout << " Failed to load the input branch '" << branch_name.str() << "_energy'\n";
		file->Close();
		return 1;
	}
	if(!b_mult){
		std::cout << " Failed to load the input branch '" << branch_name.str() << "_mult'\n";
		file->Close();
		return 1;
	}

	double short_value, long_value, energy_value;
	TFile *out_file = NULL, *cut_file = NULL;
	TTree *out_tree = NULL, *cut_tree = NULL;
	TCutG *cut = NULL;
	if(!debug){
		if(use_cut){
			std::cout << " Loading the cut\n";
			cut_file = new TFile("cuts.root", "READ");
			if(!cut_file->IsZombie()){ 
				cut = (TCutG*)cut_file->Get("cut"); 
				if(!cut){
					std::cout << " Failed to load the cut 'cut'\n";
					use_cut = false; 
					cut_file->Close();
				}
			}
			else{ 
				std::cout << " Failed to open the cut file 'cuts.root'\n";
				use_cut = false; 
				cut_file->Close();
			}
		}

		std::cout << " Opening output file 'analyzer.root'\n";
		out_file = new TFile("analyzer.root", "RECREATE");
		
		// Standard tree
		out_tree = new TTree(argv[2], "Pulse analysis tree");
		out_tree->Branch((branch_name.str() + "_short").c_str(), &short_value);
		out_tree->Branch((branch_name.str() + "_long").c_str(), &long_value);
		out_tree->Branch((branch_name.str() + "_energy").c_str(), &energy_value);
		
		// Cut tree
		if(use_cut){
			cut_tree = new TTree("gated", "Pulse analysis tree (gated)");
			cut_tree->Branch("cut_short", &short_value);
			cut_tree->Branch("cut_long", &long_value);
			cut_tree->Branch("cut_energy", &energy_value);
		}
	}

	// Get the pulse size
	for(unsigned int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		if(mult == 0){ continue; }
		else{ 
			wave_size = wave.size()/mult;
			break; 
		}
	}
	std::cout << " Using wave size " << wave_size << std::endl;

	// Canvas
	TCanvas *can1 = new TCanvas("can1", "canvas");
	TCanvas *can2 = new TCanvas("can2", "canvas");
	can1->cd();

	// Histogram
	TH2D *hist = NULL;
	TH1D *hist2 = NULL;
	unsigned int status_update = 100000;
	if(method == 0){ // slope inversion
		hist = new TH2D("hist", method_names[0].c_str(), 250, 0, 500, 250, 0, 500);
		hist->GetXaxis()->SetTitle("Short (arb. units)");
		hist->GetYaxis()->SetTitle("Long (arb. units)");
	}
	else if(method == 1){ // constant percentage
		hist = new TH2D("hist", method_names[1].c_str(), 250, 0, 500, 250, 0, 500);
		hist->GetXaxis()->SetTitle("Short (arb. units)");
		hist->GetYaxis()->SetTitle("Long (arb. units)");
	}
	else if(method == 2){ // slope variation
		hist = new TH2D("hist", method_names[2].c_str(), 250, 0, 500, 250, 0, 500);
		hist->GetXaxis()->SetTitle("Short (arb. units)");
		hist->GetYaxis()->SetTitle("Long (arb. units)");
	}
	else if(method == 3){ // pulse height ratio
		hist = new TH2D("hist", method_names[3].c_str(), 250, 0, 500, 250, 0, 500);
		hist->GetXaxis()->SetTitle("Short (arb. units)");
		hist->GetYaxis()->SetTitle("Long (arb. units)");
	}
	else if(method == 4){ // fast fit RMS
		hist = new TH2D("hist", method_names[4].c_str(), 200, 0, 100, 200, 0, 100);
		hist->GetXaxis()->SetTitle("Fit RMS (arb. units)");
		hist->GetYaxis()->SetTitle("Real RMS (arb. units)");
	}
	else if(method == 5){ // fast fit chi^2
		hist = new TH2D("hist", method_names[5].c_str(), 400, 0, 8000, 200, 0, 2000);
		hist2 = new TH1D("hist2", method_names[5].c_str(), 200, 0, 50);
		hist->GetXaxis()->SetTitle("Chi^2 (arb. units)");
		hist->GetYaxis()->SetTitle("Energy (arb. units)");
		hist2->GetXaxis()->SetTitle("Absolute Difference (arb. units)");
		hist2->GetYaxis()->SetTitle("Counts per bin");
		status_update = 50000;
	}
	else if(method == 6){ // root pulse fit chi^2
		hist = new TH2D("hist", method_names[6].c_str(), 200, 0, 2000, 200, 0, 2000);
		hist2 = new TH1D("hist2", method_names[6].c_str(), 200, 0, 10);
		hist->GetXaxis()->SetTitle("Chi^2 (arb. units)");
		hist->GetYaxis()->SetTitle("Energy (arb. units)");
		hist2->GetXaxis()->SetTitle("Absolute Difference (arb. units)");
		hist2->GetYaxis()->SetTitle("Counts per bin");
		status_update = 50000;
	}
	else if(method == 7){ // root full pulse fit
		hist = new TH2D("hist", method_names[6].c_str(), 200, 0, 2000, 200, 0, 2000);
		hist2 = new TH1D("hist2", method_names[6].c_str(), 200, 0, 10);
		hist->GetXaxis()->SetTitle("Chi^2 (arb. units)");
		hist->GetYaxis()->SetTitle("Energy (arb. units)");
		hist2->GetXaxis()->SetTitle("Absolute Difference (arb. units)");
		hist2->GetYaxis()->SetTitle("Counts per bin");
		status_update = 10000;
	}
	hist->SetStats(false);

	std::vector<double> short_integral, long_integral;
	std::vector<double>::iterator iter1, iter2, iter3;
	unsigned int count = 0;
	
	std::ifstream in_params("method6param.dat");
	std::ofstream out_params("dumb.dat");
	//std::ifstream in_params("dump.dat");
	//std::ofstream out_params("parameters.dat");
	
	long time_holder1 = 0;
	double time_holder2 = 0.0;
	clock_t cpu_time = clock();
	time_t real_time;
	time(&real_time);
	
	//unsigned int num_entries = 10000;
	unsigned int num_entries = tree->GetEntries();
	if(debug){ num_entries = 100; }
	
	std::cout << " Processing " << num_entries << " entries\n";
	for(unsigned int i = 0; i < num_entries; i++){
		tree->GetEntry(i);
		if(i % status_update == 0){ 
			time_holder1 += (long)(clock()-cpu_time);
			time_holder2 += difftime(time(NULL), real_time);
			cpu_time = clock();
			time(&real_time);
			if(i != 0){ 
				if(time_holder2 < 1){ time_holder2 = 1; } // Prevent zero time remaining
				std::cout << " Entry no. " << i << ", CPU = " << ((float)time_holder1)/CLOCKS_PER_SEC << " s, REAL = " << time_holder2;
				std::cout << " s, REMAIN = " << ((float)(num_entries-i))*(time_holder2/i) << " s\n";
			}
		}
		count += integrate(wave, wave_size, short_integral, long_integral, method, debug, can1, can2, true, in_params, out_params);
		//count += integrate(wave, wave_size, short_integral, long_integral, method, debug, can1, can2, false, in_params, out_params);
		energy_value = *energy.begin();
		for(iter1 = short_integral.begin(), iter2 = long_integral.begin(); iter1 != short_integral.end() && iter2 != long_integral.end(); iter1++, iter2++){
			short_value = (*iter1);
			long_value = (*iter2);
			if(method == 3){ hist->Fill(short_value, energy_value); }
			else if(method < 5){ hist->Fill(short_value, long_value); }
			else{ 
				if(method == 7){ hist->Fill(short_value, energy_value); }
				else{ hist->Fill(long_value, energy_value); }
				hist2->Fill(short_value);
			}
			if(!debug){
				out_tree->Fill();
				if(use_cut){
					if(cut->IsInside(short_value, long_value)){
						cut_tree->Fill(); 
					}
				}
			}
		}
	}
	
	std::cout << " Found " << count << " pulses in " << num_entries << " tree entries\n";
	if(!debug){ 
		std::cout << " Wrote " << out_tree->GetEntries() << " entries to the output tree\n";
		if(use_cut){ std::cout << " Wrote " << cut_tree->GetEntries() << " entries to the gated tree\n"; }
	}

	can1->cd();
	hist->Draw("COLZ");
	can1->Update();
	
	if(hist2){
		can2->cd();
		hist2->Draw();
		can2->Update();
		can2->WaitPrimitive();
	}
	else{ can1->WaitPrimitive(); }

	if(!debug){
		out_file->cd();
		out_tree->Write();
		if(use_cut){ 
			cut_tree->Write(); 
			cut_file->Close();
		}
		hist->Write();
		out_file->Close();
	}
	
	can1->Close();
	if(hist2){ can2->Close(); }
	file->Close();
	in_params.close();
	out_params.close();
	
	rootapp->Delete();
		
	return 0;
}

// For CINT
//int PulseAnalyzer(int argc, char* argv[]){ return main(argc, argv); }
