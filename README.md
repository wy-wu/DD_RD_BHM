# DD_RD_BHM
The repository contains the following files:

theory_AC.R: This script provides the foundational R functions for calculating the theoretical accuracy for each combination of decision strategy and evaluation criterion. 
	     The primary functions are:
	     DDEC_S: Deterministic decision (DD)’s accuracy under evaluation criterion based on the sub-population mean (EC-S)
	     RDEC_S: Randomized decision (RD)’s accuracy under EC-S
	     DDEC_P: DD’s accuracy under evaluation criterion based on the population mean (EC-P)
	     RDEC_P: RD’s accuracy under EC-P

Figure2.R: This script compares the theoretical accuracy (calculated using theory_AC.R) with simulation-based accuracy and plots the results. 
	   The script explores five distinct scenarios characterized by varying population variances and sub-population variances 


acc.calculate.r: This script is used for the simulated case studies described in the accompanying manuscript. 
		 It performs Section 4 – Simulated Case Studies:
		 4.1 Computes accuracy under the assumption that the true parameter values are known.
		 4.2 Performs decision-making based on observed data and estimated parameters.
		 Note: Three scenarios are evaluated; please refer to the manuscript for details.

example3.R: This script serves as a practical example for running the BHM defined in fully_bayesian_model.stan. It processes data for a specific case study (Scenario 3).

fully_bayesian_model.stan: This file defines the core statistical model—a Bayesian Hierarchical Model (BHM).
