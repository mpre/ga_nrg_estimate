ga_nrg_estimate
===============

Estimate energy variation using genetic algorithms

This exercise estimate the energy supply variation for the day n+1 given the
supply variation for the lasts 7 days, today's temperature (in celsius) and
tomorrow's estimate temperature.

energy.txt is the dataset
ga_nrg_estimate.cpp contains the code

You need galib in order to compile this code. Once you have
successfully compiled galib simply use the next command (assuming you have galib
compiled under "ga" subdirectory)

	 g++ -I. -Wall ga_nrg_estimate.cpp ga/libga.a -o exename

and run it.

