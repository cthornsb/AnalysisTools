/*************************************************************************
 *
 *  Filename: PulseAnalysis.h 
 *
 *  Description: 
 *
 *	Author(s):
 *     Michael T. Febbraro
 *     
 *
 *  Creation Date: 11/25/2012
 *  Last modified: 9/17/2013 
 *
 * -----------------------------------------------------
 * 	Nuclear Reaction Group
 *  University of Michigan, Ann Arbor, MI, USA 
 *  (c) All Rights Reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "Analysis.h"

/** Constructor  */
PulseAnalysis::PulseAnalysis()
{
	deltaT = 1;
	pulse_size = 0;
	full_size = 0;
	pulse = NULL;
	debug = false;
}

/** Destructor  */
PulseAnalysis::~PulseAnalysis()
{
	if(full_size > 0){ delete pulse; }
}

/** ----------------------------------------------------  
*	Get version
*		-  Print the version and revision date
*          
*  Notes:  
*		Keep this information up to date.
*	----------------------------------------------------
*/
void PulseAnalysis::GetVersion()
{
	std::string  version   = "Pulse analysis 1.0"; 
	std::string  revision  = "9/17/2013";
	
	std::cout << "Package: " << version << std::endl;
	std::cout << "Revision: " << revision << std::endl;
}

/** ----------------------------------------------------  
*	Optimize PSD by charge integration method
*		- Calculates discrimination parameters by 
*		  integration of user defined regions of an
*		  input pulse over a user defined range for evaluation.
*
*	Inputs:
*			pulse - input array
*			start - start index of long integral
*			stop - stop index of both integrals
*			lower - lower index range of short integral
*			upper - upper index range of short integral
*			method - see above...
*	----------------------------------------------------
*/
int PulseAnalysis::OptimizePSD (int start, int stop, int lower, int upper, double *paraL, double *paraS)
{
	if ((upper - lower) < 0) {return -1;}
	if ((upper - lower) < sizeof(paraS)/sizeof(double)) {return -1;}
	
	j = 0;
	for (i = 0; i < pulse_size; i++)
	{
		if(pulse[i] > j) { j = (int)pulse[i]; l= i;}
	}
	
	// Constant fraction discrimination (CFD) at 50% of amplitude
	
	j = l;
	for(i = j - 50; i < j; i++)
	{
		if(pulse[i] < (pulse[l])*0.5) {k = i;}
	}
	
	x = ((0.5*(pulse[l]) - pulse[k - 1])/((pulse[k + 1] 
			- pulse[k - 1])/3))+ ((double)k - 1);
			
			
	if((k - start) > 0 && (k + start) < pulse_size) {
			// Initialization
			integralL = 0;
			for ( j = 0; j < upper - lower; j++){
				 paraS[j] = 0;
			}
			
			// Begin integration using trapezoidal rule 
			for (i = (k - start); i < (k + stop); i++) {
			
				integralL += 0.5*(pulse[i-1] + pulse[i]);
					
				if (i >= (k + lower)) { 
					for ( j = 0; j < upper - lower; j++)
					{
						if (i >= (k + lower + j)) {
							 paraS[j] += 0.5*(pulse[i-1] + pulse[i]);
						}
					}
				}

			}
	}
	
	if (integralL != 0) {
		*paraL = integralL;
		*paraS = integralS;
		return 0;
	} else { return -1;}

}

/** ----------------------------------------------------  
*	Determine PSD by charge integration method
*		- Calculates discrimination parameters by 
*		  integration of user defined regions of an
*		  input pulse.
*
*	Methods:
*			1 - trapezoidal rule
*			2 - composite Simpson's rule
*			3 - rectangular method 
*			4 - Harwell method
*
*	Inputs:
*			pulse - input array
*			start - start index of long integral
*			stop - stop index of both integrals
*			offset - start index of short integral
*			method - see above...
*	----------------------------------------------------
*/
int PulseAnalysis::PSD_Integration (int method, double &paraL, double &paraS)
{	
	int start = first.left;
	int stop = first.min;
	int offset = 5;
	
	j = 0;
	for (i = 0; i < pulse_size; i++)
	{
		if(pulse[i] > j) { j = (int)pulse[i]; l= i;}
	}
	
	// Constant fraction discrimination (CFD) at 50% of amplitude
	
	j = l;
	for(i = j - 50; i < j; i++)
	{
		if(pulse[i] < (pulse[l])*0.5) {k = i;}
	}
	
	x = ((0.5*(pulse[l]) - pulse[k - 1])/((pulse[k + 1] 
			- pulse[k - 1])/3))+ ((double)k - 1);
	
	
	
	if((k - start) > 0 && (k + start) < pulse_size) {
			// Initialization
			integralS = 0; integralL = 0;
			
			// Begin integration using trapezoidal rule 
			if (method == 1) {
				for (i = (k - start); i < (k + stop); i++) {
			
					integralL += 0.5*(pulse[i-1] + pulse[i]);
					
					if (i > (k + offset)) { integralS += 0.5*(pulse[i-1] + pulse[i]);}
				}
			}
			
			// Begin integration using composite Simpson's rule
			if (method == 2) {
			
				sumL = 0;
				sumS = 0;
				integralL = pulse[k - start];
				integralL = pulse[k + offset];
				
				for (i = (k - start) + 1; i < (k + stop) - 2; i+=2) {
			
					sumL += pulse[i];
					
					if (i > (k + offset)) { sumS += pulse[i];}
				}
				
				integralL += 4*sumL;
				integralS += 4*sumS;
				
				sumL = 0;
				sumS = 0;
				for (i = (k - start + 2); i < (k + stop) - 3; i+=2) {
			
					sumL += pulse[i];
					
					if (i > (k + offset)) { sumS += pulse[i];}
				}
				
				integralL += 2*sumL;
				integralS += 2*sumS;
				
				integralL += pulse[k + stop];
				integralS += pulse[k + stop];
				
				integralL = integralL/3;
				integralS = integralS/3;
			}
			
			// Begin integration using rectangular method 
			if (method == 3) {
				for (i = (k - start); i < (k + stop); i++) {
			
					integralL += pulse[i];
					
					if (i > (k + offset)) { integralS += pulse[i];}
					//if (i > (k - start) && i <= (k + offset)) { integralS += pulse[i];}
				}
				
				// Adjust for error due to array indexing (i.e. rounding)
				if (x < k + offset)
				{
					//integralS =+ ((pulse[k + offset] - pulse[k + offset-1])*x*0.5) + ((x - (double)(k + offset))*pulse[k + offset-1]); 
				}
				
			}
			
	} else {return -1;}
	
	if (integralS != 0) {
		paraL = integralL;
		paraS = integralS;
		return 0;
	} else { return -1;}
}

/** ----------------------------------------------------  
*	Restores the baseline of a input pulse
*
*	Methods:
*			1 - SNIP fitting routine
*			2 - baseline averaging
*			 
*	Inputs:
*			pulse - input array
*			iterations - number of iterations
*			method - see above...
*
*	----------------------------------------------------
*/
int PulseAnalysis::Baseline_restore (double *baseline, int iterations, int method)
{
	if (method == 1)
	{
		for (i=0; i< pulse_size; i++) {
			baseline[i] = pulse[i] + 0.3;
		}
		
		// Michael Febbraro 3/2/2013
		std::cout << "Method Disabled!" << std::endl; return -1;
	
		//TSpectrum::Background(baseline, pulse_size, iterations,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,
		//		TSpectrum::kTRUE, TSpectrum::kBackSmoothing3,TSpectrum::kFALSE);

		for (i=0; i< pulse_size; i++) {
			pulse[i] = pulse[i] - baseline[i];
		}
	}
	
	if (method == 2)
	{
		x = 0;
		for (i = 0; i < iterations; i++) {x += pulse[i];}
		x = x/(double)iterations;
		
		for (i = 0; i < pulse_size; i++)
		{
			pulse[i] -= x;
		}
	}
	
	if (method == 3)
	{
		j = 0;
		for (i = 0; i < pulse_size; i++)
		{
			if(pulse[i] > j) { j = (int)pulse[i]; k = i;}
		}
	
		x = 0; y = 0;
		for (i = 0; i < iterations; i++) {x += pulse[i]; y += pulse[pulse_size - i];}
		x = x/(double)iterations; y = y/(double)iterations;
		
		for (i = 0; i < pulse_size; i++)
		{
			if(i >= k) { pulse[i] -= y;}
			else { pulse[i] -= x; }
		}
	}
	
	return 0;
}

/** ----------------------------------------------------  
*	Determine time characteristics of the a pulse
*		- Fits timing regions using linear regression
*		  over a user defined range.
*
*	Inputs:
*			pulse - input array
*			range - range of linear regression for timing
*	----------------------------------------------------
*/
int PulseAnalysis::Parameters (int range, double &CFD, double &amplitude, double &risetime, double &falltime, double &width)
{
	j = 0;
	for (i = 0; i < pulse_size; i++)
	{
		if(pulse[i] > j) { j = (int)pulse[i]; k = i;}
	}
	
	amplitude = pulse[k];
	
	// Constant fraction discrimination (CFD) at 50% of amplitude
	//memset(linear, 0, sizeof(linear));
	
	j = k;
	for(i = j - 50; i < j; i++)
	{
		if(pulse[i] < (amplitude)*0.5) {k = i;}
	}
	
	/*
	for(i = (point_t - (int)(range/2)); i <= (point_t - (int)(range/2) + range); i++) {
		linear[0] += i*pulse[i];
		linear[1] += i;
		linear[2] += pulse[i];
		linear[3] += i*i;
	}
	linear[0] = linear[0]/range;
	linear[1] = linear[1]/range;
	linear[2] = linear[2]/range;
	linear[3] = linear[3]/range;
	
	*CFD = ((0.5*(*amplitude) - pulse[(point_t - (int)(range/2))])/((linear[0] 
		- linear[1]*linear[2])/(linear[3] - linear[1]*linear[1])))+ (point_t - (int)(range/2));
		
	*/
	
	CFD = ((0.5*(amplitude) - pulse[(k - (int)(range/2))])/((pulse[k + (int)(range/2)] 
			- pulse[k - (int)(range/2)])/(double)range))+ ((double)k - ((double)range/2));
		
	// Risetime by 10-90% Method
	j = k;
	for(i = j - 50; i < j; i++)
	{
		if(pulse[i] < (amplitude)*0.9) {k = i;}
	}
	
	//memset(linear, 0, sizeof(linear));
	for(i = (k - (int)(range/2)); i <= (k - (int)(range/2) + range); i++) {
		linear[0] += i*pulse[i];
		linear[1] += i;
		linear[2] += pulse[i];
		linear[3] += i*i;
	}
	linear[0] = linear[0]/range;
	linear[1] = linear[1]/range;
	linear[2] = linear[2]/range;
	linear[3] = linear[3]/range;
	
	risetime = ((0.1*(amplitude) - pulse[(k - (int)(range/2))])/((linear[0] 
		- linear[1]*linear[2])/(linear[3] - linear[1]*linear[1]))) + (k - (int)(range/2));
		
	j = k;
	for(i = j - 50; i < j; i++)
	{
		if(pulse[i] < (amplitude)*0.1) {k = i;}
	}
	
	//memset(linear, 0, sizeof(linear));
	for(i = (k - (int)(range/2)); i <= (k - (int)(range/2) + range); i++) {
		linear[0] += i*pulse[i];
		linear[1] += i;
		linear[2] += pulse[i];
		linear[3] += i*i;
	}
	linear[0] = linear[0]/range;
	linear[1] = linear[1]/range;
	linear[2] = linear[2]/range;
	linear[3] = linear[3]/range;
	
	risetime = risetime - ((0.9*(amplitude) - pulse[(k - (int)(range/2))])/((linear[0] 
		- linear[1]*linear[2])/(linear[3] - linear[1]*linear[1]))) + (k - (int)(range/2));

		
	// NOTE: Still need to do fall time and width...
	falltime = 0;
	width = 0;
		
	return 0;
}

/** ----------------------------------------------------  
*	Determine location and number of peaks in a waveform
*		- Determines the number of peaks by applying a
*		  above-threshold condition at a user defined 
*		  number of std. dev. of the baseline all using
*		  the first derv. of the input waveform.  
*
*	Inputs:
*			pulse - input array
*			sigma - threshold in number of std. deviations
*			range - number of points to deterimine std. dev.
*	Outputs:
*			numPeaks - number of peaks found
*			locPeaks - array of array indexes of peak locations
*
*	----------------------------------------------------
*/
int PulseAnalysis::PeakFinder (int sigma, int range, int method, int &numPeaks, int &value)
{
	for(i = 0; i < 10; ++i) 
	{
		//locPeaks[i] = new int[3];
	}
	
	// Calculate first derivative and average of first derv.
	x = 0;
	for (i = 0; i < pulse_size - 1; i++)
	{
		pulse_temp[i] = pulse[i + 1] - pulse[i];
		x =+ pulse_temp[i];
	}
	
	
	// Determine standard deviation of baseline
	x = x/range;
	z = 0;
	for (i = 0; i < range; i++) { z =+ pow((pulse_temp[i] - x), 2); }
	z = sqrt(z/range);
	
	numPeaks = 0;
	
	if(method == 1) { w = sigma*z; }
	
	if(method == 2) { w = sigma; }
	
	j = 0; k = 0; x = 0; y = 0; m = 0; n = 0;
	positive = 0; negative = 0;
	for (i = 0; i < pulse_size - 5; i++)
	{	
		if( pulse_temp[i] <= w && pulse_temp[i] >= -1*w)
		{	
			if (positive == 1 && negative == 1)
			{
				positive = 0;
				negative = 0;
				//locPeaks[*numPeaks - 1][1] = j;
				//locPeaks[*numPeaks - 1][2] = k;
				//locPeaks[*numPeaks - 1][0] = n;
				value = k;
				
				x = 0; y = 0;
			}
		}	
		else
		{
				
			if (positive == 0 && negative == 0 && i > n + 40) { m++; n = i;}
		
			// Positive inflection point
			if(pulse_temp[i] > 0) 
			{	
				positive = 1;
				if (pulse_temp[i] > x) {x = pulse_temp[i]; j = i;}
			}
			
			// Negative inflection point
			if(pulse_temp[i] < 0) 
			{
				negative = 1;
				if (pulse_temp[i] < y) {y = pulse_temp[i]; k = i;}
			}
			
		}		
		
	}
	
	numPeaks = m;
	
	return 0;
}

/** ----------------------------------------------------  
*	Determine the nth derivative of a waveform
*
*	Inputs:
*			pulse - input array
*			order - order of the dervative (i.e. 1 = first, 2 = second, ...)
*	----------------------------------------------------
*/
int PulseAnalysis::Derivative(int order)
{
	for(j = 1; j <= order; j++)
	{	
		for (i = 0; i < pulse_size - j; i++)
		{
			pulse[i] = pulse[i + 1] - pulse[i];
		}
	}
	return 0;
}

/** ----------------------------------------------------  
*	Determine the integral of a waveform
*
*	Inputs:
*			pulse - input array
*	----------------------------------------------------
*/
int PulseAnalysis::Integral(double *pulse)
{
	for (i = 0; i < pulse_size - 1; i++)
	{
			pulse[i + 1] = pulse[i + 1] + pulse[i];
	}
	return 0;
}

/** ----------------------------------------------------  
*	Output a time pickoff from user defined condition
*
*	Inputs:
*			pulse - input array
*			
*	----------------------------------------------------
*/
int PulseAnalysis::Time_Pickoff (int range, int low, int high , int method, double &CFD)
{
	double max, time;
	int index;
	if (method == 1)
	{	
	max = 0;
	
	for (i = low; i <= high; i++)
	{
		if (pulse[i] > max)
		{
			max = pulse[i];
			index = i;
		}		
	}	

	// Constant fraction discrimination (CFD) at 50% of amplitude
	
	for(i = low; i < index; i++)
	{
		if(pulse[i] < (max*0.5)) {k = i;}
	}

	time = (((max*0.5) - pulse[k])/(pulse[k+1] - pulse[k])) + (double)k;
	
	if (time >= low && time <= high){CFD = time;}
	else {time = -1;}
	
	}
	
	if (method == 2)
	{
	
	max = 0;
	
	// Find first value over 'range'
	for (i = low; i <= high; i++)
	{
		if (pulse[i] > range)
		{
			k = i; break;
		}		
	}	
	
	time = (((((double)(range)) - pulse[k-1]))/(pulse[k] - pulse[k-1])) + (double)k;
	
	if (time >= low && time <= high){CFD = time;}
	else {time = -1;}
	
	
	}
	
	return 0;
}

int PulseAnalysis::PSD_Zerocross (int integration, int differentiation, double &PSD)
{
	j = 0;
	for (i = 0; i < pulse_size; i++)
	{
		if(pulse[i] > j) { j = (int)pulse[i]; k = i;}
	}
	
	// CFD timing pickoff at 50% of amplitude
	//memset(linear, 0, sizeof(linear));
	
	for(i = (k - (int)(3/2)); i <= (k - (int)(3/2) + 3); i++) {
		linear[0] += i*pulse[i];
		linear[1] += i;
		linear[2] += pulse[i];
		linear[3] += i*i;
	}
	linear[0] = linear[0]/3;
	linear[1] = linear[1]/3;
	linear[2] = linear[2]/3;
	linear[3] = linear[3]/3;
	
	PSD = ((0.5*(pulse[k]) - pulse[(k - (int)(3/2))])/((linear[0] 
		- linear[1]*linear[2])/(linear[3] - linear[1]*linear[1])))+ (k - (int)(3/2));
		
	// Shaping amplifier 
	//memset(pulse_temp, 0, sizeof(pulse_temp));
	for(i = (int)integration; i < (pulse_size - integration); i++) {
		for(j = 0; j < integration; j++) {
			if((j + i) >= 0 && (j + i) < pulse_size) {pulse_temp[i] = pulse_temp[i] + pulse[j + i];}
		}
	}
	
	//memset(pulse, 0, sizeof(pulse));
	//for (i = 0; i < (pulse_size - 1); i++) {pulse[i] = pulse_temp[i + 1] - pulse_temp[i];}
	
	
	//PSD = TMath::LocMin(pulse_size,pulse) - TMath::LocMax(pulse_size,pulse);
	
	
	//memset(pulse_temp, 0, sizeof(pulse_temp));
	for (i = 0; i < (pulse_size - 1); i++) {pulse_temp[i] = pulse[i + 1] - pulse[i];}
	
	//for (i=0; i< pulse_size; i++) {
	//	pulse[i] = pulse_temp[i];
	//}
	return 0;
}

/** ----------------------------------------------------  
*	Determine the integral below half the pulse amplitude.  Used for
*   reconstruction of partial pulses.  
*
*	Inputs:
*			pulse - input array
*			integral - integral below half the pulse amplitude
*	----------------------------------------------------
*/
int PulseAnalysis::Half_Integral(double &integral)
{
	integral = 0;
	j = 0;
	for (i = 0; i < pulse_size; i++)
	{
		if(pulse[i] > j) { j = (int)pulse[i]; k = i;}
	}

	for (i = 0; i < pulse_size ; i++)
	{
		if(pulse[i] < (0.5*pulse[k]))
		{
			integral = pulse[i + 1] + pulse[i];
		}
		else
		{
			integral = pulse[i + 1] +pulse[k];
		}
	}
	
	return 0;
}

/** ----------------------------------------------------  
*	Smoothing method for waveforms
* 
*	Methods:
*			1 - Moving Average (option = range of moving average)
*
*	Inputs:
*			pulse - input array
*			iterations - number of iterations to apply smoothing method
*			method - see above
*			option - see above
*	----------------------------------------------------
*/
int PulseAnalysis::Smooth(int iterations, int method, double option)
{
	for(k = 0; k < iterations; k++)
	{
	// Moving average
	if (method == 1)
	{
		for (i = 2; i < pulse_size - 3; i++)
		{
			MovingAverage[i] = 0;
			for (j = i - 2; j <= i + 2; j++) {MovingAverage[i] += pulse[j];}
		}
		for (i =1; i < pulse_size - 3; i++)
		{
			pulse[i] = (MovingAverage[i]/5) ;
		}
	}

	}
	return 0;
}

/** ----------------------------------------------------  
*	HPGe
*	- This method processes signals from HPGe detectors.
*
*	Inputs:
*			pulse - input array
*			amplitude - amplitude of the pulse
*	----------------------------------------------------
*/
int PulseAnalysis::HPGe(double &amplitude){
	return 0;
}

unsigned int PulseAnalysis::PreProcess(std::vector<int> &input){
	// Copy the input vector into an array of doubles for later processing
	if(full_size > 0){ delete pulse; }
	full_size = 0;
	pulse = new double[input.size()];
	for(std::vector<int>::iterator iter = input.begin(); iter != input.end(); iter++){
		pulse[full_size] = *iter;
		full_size++;
	}
	if(debug){ std::cout << "pulse_size = " << pulse_size << ", full_size = " << full_size << ", mult = " << full_size/pulse_size << std::endl; }
	return full_size;
}

unsigned int PulseAnalysis::Process(unsigned int start, unsigned int stop){
	if(full_size == 0 || !pulse){ return 0; }
	
	// Do baseline correction. Use the range to the left of the pulse.
	// The right side of the pulse may have other pulses
	double base_sum = 0.0;
	unsigned int base_max = start + (unsigned int)(0.1*pulse_size); // First 10% of range
	for(unsigned int i = start; i < base_max; i++){
		base_sum += pulse[i];
	}

	double base_line = base_sum/(base_max-start);
	for(unsigned int i = start; i < stop; i++){
		pulse[i] = pulse[i] - base_line;
	}

	// Find the global maximum
	double maximum = -9999.0;
	unsigned int maximum_bin = 0;
	for(unsigned int i = start; i < stop; i++){
		if(pulse[i] > maximum){ 
			maximum = pulse[i];
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
		back_sub_left = pulse[i] - 0.1*maximum;
		back_sub_right = pulse[i+1] - 0.1*maximum;
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
				first.left = leading_edge;
				first.max = peak_bin;
				first.max_value = pulse[peak_bin];
			}
			else if(num_peaks == 2){
				second.left = leading_edge;
				second.max = peak_bin;
				second.max_value = pulse[peak_bin];
				break; // That's enough!
			}
		}
	}

	// Find the local minimum
	if(num_peaks == 1){ // Single pulse, minimum must occur to the right of the pulse
		for(unsigned int i = first.max; i < stop; i++){
			if(pulse[i] < first.min_value){
				first.min_value = pulse[i];
				first.min = i;
			}
		}
	}
	else{ // Multiple pulses, minimum must occur between pulse peaks
		for(unsigned int i = first.max; i < second.max; i++){
			if(pulse[i] < first.min_value){
				first.min_value = pulse[i];
				first.min = i;
			}
		}
	}

	// Find 80% of the peak maximum
	for(unsigned int i = first.left; i <= first.max; i++){
		if(pulse[i] >= 0.8*first.max_value){
			if(i > 0){ first.cfd = i-1; }
			else{ first.cfd = i; }
			break;
		}
	}
	if(first.cfd == first.max){ first.cfd--; }
	
	// Print debug information
	if(debug){ 
		std::cout << "Global: baseline = " << base_line << ", maximum_bin = " << maximum_bin;
		std::cout << ", maximum = " << maximum << ", num_peaks = " << num_peaks << std::endl;
		std::cout << "Peak: " << first.print() << std::endl << std::endl;
	}
	
	return num_peaks;
}
