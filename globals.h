#ifndef MY_GLOBALS_H
#define MY_GLOBALS_H

// Const values used in the simulations
const Double_t 
	pi       	= TMath::Pi(),
 	p0       	= 3,
 	r0       	= 6.62,      // fm
 	a        	= 0.542,    // fm
 	sigma    	= 6.5,      // fm^2
 	radiusSq 	= sigma / pi, // fm^2
 	minDisSq 	= 0,          // fm^2
 	nucleons    = 208,
 	simulations = 1e6;
 


// typedef's for avoiding using long type names
typedef std::vector<Double_t>                           vec_1d;
typedef std::vector<std::vector<Double_t>>              vec_2d;
typedef std::vector<std::vector<std::vector<Double_t>>> vec_3d;

// Necessary vectors. Stores the positions multiple simulations of the nucleons for each value of impact parameter
vec_3d xPosSim, yPosSim, collisions;
vec_1d b;

// Functions



#endif
