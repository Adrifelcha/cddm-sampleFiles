# cddm-sampleFiles

## Sample code to introduce the "cddm" JAGS module

### root/CDDM_basicNotes.pdf

A brief and didactic presentation of the Circular Drift Diffussion Model.

### root/Functions/

A set of Rscripts to load specific functions build to generate parameter values, simulate data, load and run a JAGS model using the cddm module, and check the posterior samples.

- fun_simData.R - A set of functions to generate and plot the random walk process described by the CDDM and corresponding bivariate data (choices in radians, and RT in seconds).
- fun_genParameters.R - A set of functions to generate random values for the parameters used in the cddm module. This Rscript also contain the functions prepared to compute relevant variable transformations.

### /simple-example/

R scripts:
- toyData.R - An Rscript to load (or generate, as needed) the toy dataset.
- toyExample.R -  A simple Rscript to load the toy dataset and run the cddm jags module

Objects:
- toyData.csv - A dataset with 200 observations (choice in radians and RTs)

### /full-example/

R scripts:
- getData.R - An Rscript to load or generate a simulated dataset from a set of random parameter values.

Objects:
- data.csv  - Simulated dataset
- trueParValues.RData - R object containing the parameter values used in the simulation.

