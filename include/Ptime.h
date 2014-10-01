#ifndef PTIME_H
#define PTIME_H

#include <sstream>
#include <stdlib.h>

struct time3{
	int h, m, s;
	
	time3(){
		set(0, 0, 0);
	}
	
	time3(int h_, int m_, int s_){
		set(h_, m_, s_);
	}
	
	time3 &operator=(const time3 &other){
		h = other.h;
		m = other.m;
		s = other.s;
	}
	
	// Set the time using hours minutes and seconds
	void set(int h_, int m_, int s_){
		h = h_;
		m = m_;
		s = s_;
		correct();
	}
	
	// Set the time using floats
	void set(float h_, float m_, float s_){
		h = (int)h_;
		m = (int)m_;
		s = (int)s_;
		correct();
	}
	
	// Check the time for minutes, seconds >= 60
	void correct(){
		if(s > 59){ m += (s / 60); s = (s % 60); } // correct for seconds >= 60
		if(m > 59){ h += (m / 60); m = (m % 60); } // correct for minutes >= 60
		if(h > 24){ h = (h % 24); } // correct for hours > 24
	}
	
	// Return the time in hours
	float getHours(){
		return(h + m/60.0 + s/3600.0);
	}
	
	// Return the time in minutes
	float getMinutes(){
		return(60.0*h + m + s/60.0);
	}
	
	// Return the time in seconds
	float getSeconds(){
		return(3600.0*h + 60.0*m + s);
	}
	
	// Extract the time from a formatted string (hh:mm:ss)
	float extract(std::string input){
		std::string sub_h = "", sub_m = "", sub_s = "", overflow = "";
		int index = 0;
		for(int i = 0; i < input.size(); i++){
			if(input[i] == ':'){ 
				index++;
				continue; 
			}
			else if(index == 0){ sub_h += input[i]; }
			else if(index == 1){ sub_m += input[i]; }
			else if(index == 2){ sub_s += input[i]; }
			else{ overflow += input[i]; }
		}
		
		h = (int)(atoi(sub_h.c_str()));
		m = (int)(atoi(sub_m.c_str()));
		s = (int)(atoi(sub_s.c_str()));
		correct();
	}
		
	// Return a formated time string
	std::string print(){
		std::stringstream output;
		if(h < 10){ output << "0" << h; }else{ output << h; }
		if(m < 10){ output << ":0" << m; }else{ output << ":" << m; }
		if(s < 10){ output << ":0" << s; }else{ output << ":" << s; }
		
		return output.str();
	}
};

// Convert a time in seconds to a time string with format hh:mm:ss
std::string ConvTime(int myTime){ 
	int hrs = myTime / 3600;
	int min = myTime % 3600;
	int sec = min % 60;
	min = min / 60;	
	
	std::stringstream output;
	if(hrs < 10){ output << "0" << hrs; }else{ output << hrs; }
	if(min < 10){ output << ":0" << min; }else{ output << ":" << min; }
	if(sec < 10){ output << ":0" << sec; }else{ output << ":" << sec; }
	
	return output.str();
}

time3 DiffTime(time3 &start, time3 &stop){
	float t1 = start.getSeconds();
	float t2 = stop.getSeconds();
	float diff;
	if(t1 > t2){ // Time wraps past midnight
		diff = (86400 - t1);
		diff += t2;
	}
	else{ // Times are in the same day
		diff = t2 - t1;
	}
	
	time3 temp;
	temp.set(0.0, 0.0, diff);
	return temp;
}

#endif
