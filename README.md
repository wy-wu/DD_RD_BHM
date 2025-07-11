# Deterministic Versus Randomized Decision Making Under Bayesian Hierarchical Model

## Contents

- **theory_AC.R**  
  Provides foundational R functions for calculating the theoretical accuracy of each combination of decision strategy and evaluation criterion.  
  **Primary functions:**  
  - `DDEC_S`: Deterministic decision (DD) accuracy under evaluation criterion based on the sub-population mean (EC-S)  
  - `RDEC_S`: Randomized decision (RD) accuracy under EC-S  
  - `DDEC_P`: DD accuracy under evaluation criterion based on the population mean (EC-P)  
  - `RDEC_P`: RD accuracy under EC-P  

- **Figure2.R**  
  Compares the theoretical accuracy (from `theory_AC.R`) with simulation-based accuracy and generates visualizations.  
  This script explores five scenarios characterized by varying population variances and sub-population variances as described in paper.

- **acc.calculate.r**  
  Performs simulated case studies as described in the accompanying manuscript.  
  **Key functionality:**  
  - Section 4.1: Computes accuracy assuming true parameter values are known  
  - Section 4.2: Conducts decision-making based on observed data and estimated parameters  
  *Note: Three scenarios are evaluated; see the manuscript for details.*

- **example3.R**  
  A practical example script that demonstrates how to run the Bayesian Hierarchical Model (BHM) defined in `fully_bayesian_model.stan`. This example corresponds to Scenario 3.

- **fully_bayesian_model.stan**  
  Stan code specifying the core Bayesian Hierarchical Model.
