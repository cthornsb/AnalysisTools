// Cory R. Thornsberry
// Oct. 1st, 2014
// file: DiffTime.cpp
// Calculates the difference between two times (separated by a maximum of 23:59:59)

#include <iostream>

#include "Ptime.h"

// For compilation
int main(int argc, char* argv[]){
	if(argc < 3){
		std::cout << " Error! Invalid number of arguments. Expected 2, received " << argc-1 << "\n";
		std::cout << "  SYNTAX: ./DiffTime {start} {stop}\n";
		return 1;
	}
		
	time3 start, stop;
	std::stringstream ss1, ss2;
	
	ss1 << argv[1];
	start.extract(ss1.str());

	ss2 << argv[2];
	stop.extract(ss2.str());

	std::cout << " Start: " << start.print() << std::endl;
	std::cout << " Stop: " << stop.print() << std::endl;
	std::cout << " Diff: " << DiffTime(start, stop).print() << std::endl;

	return 0;
}
