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

#endif
