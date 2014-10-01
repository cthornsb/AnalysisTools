// Cory R. Thornsberry
// Oct. 1st, 2014
// file: Qcalc.cpp
// Calculates various beam parameters for the FN tandem at Notre Dame

#include <iostream>
#include <cmath>

#include "Ptime.h"

// For compilation
int main(){
	int option = 0;
	while(true){
		std::cout << " Please select parameter to calculate...\n";
		std::cout << "  1 - Beam current I (Amps)\n";
		std::cout << "  2 - Beam charge Q (Coulombs)\n";
		std::cout << "  3 - Necessary run time dT (Seconds)\n";
		std::cout << "   Enter selection: "; std::cin >> option;
		
		if(option == 0 || option > 3){
			std::cout << "   Invalid selection, try again\n\n";
		}
		else{ break; }
	}

	int S, Q;
	float I;
	time3 dT, start, stop;
	std::cout << "\n Enter integrator scale S: "; std::cin >> S;

	if(option == 1){
		// Calculate beam current in amps (requires Q and dT)
		std::string temp;
		std::cout << " Enter beam charge Q (Coulombs): "; std::cin >> Q;
		std::cout << " Enter run start time (hh:mm:ss): "; std::cin >> temp; start.extract(temp);
		std::cout << " Enter run stop time (hh:mm:ss): "; std::cin >> temp; stop.extract(temp);
		dT = DiffTime(start, stop);
		
		I = Q / ((pow(10.0, 2+S)) * dT.getSeconds());
	}
	else if(option == 2){
		// Calculate beam charge in coulombs (requires I and dT)
		std::string temp;
		std::cout << " Enter beam current I (Amps): "; std::cin >> I;
		std::cout << " Enter run start time (hh:mm:ss): "; std::cin >> temp; start.extract(temp);
		std::cout << " Enter run stop time (hh:mm:ss): "; std::cin >> temp; stop.extract(temp);
		dT = DiffTime(start, stop);
		
		Q = I * (pow(10.0, 2+S)) * dT.getSeconds();
	}
	else if(option == 3){
		// Calculate beam time required in seconds (requires Q and I)
		std::cout << " Enter beam charge Q (Coulombs): "; std::cin >> Q;
		std::cout << " Enter beam current I (Amps): "; std::cin >> I;
		
		dT.set(0.0, 0.0, (float)(Q / ((pow(10.0, 2+S)) * I)));
	}
	else{
		std::cout << " How did you get here???\n";
		Q = 0; I = 0.0;
	}

	std::cout << "\n Beam parameters--------------------\n";
	std::cout << "  Beam Current = " << I << " A\n";
	std::cout << "  Beam Charge = " << Q << " C\n";
	std::cout << "  Time Taken = " << dT.print() << std::endl;

	return 0;
}
