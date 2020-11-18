# plant_sphingolipid_kinetic_model

This repository contains the files required to generate a kinetic model of the sphingolipid pathway in plants. The original code was published by [Tan et al.](https://www.sciencedirect.com/science/article/abs/pii/S1096717610000996) 

The following modifications have been made to the original code:
(1) Added a module to automate parameter screening and to output models consistent with previously reported network responses to different perturbations.
(2) Amended the process of adding activation reactions. The original code results in mass imbalances and did not account for the activated pathway when calculating the overall flux of the activated reaction.
(3) Introduced an additional layer to sample the initialization of metabolite concentration. This is done to obtain wider coverage of the solution space.

To regenerate the results in the paper, run the code 'MainScript.m'. The kinetic parameters associated with passing models will be in the output matrix 'CuteKvec_passed'. 

To implement other regulatory schemes. The regulatory matrix 'Sreg' which is part of the data structure 'Net' needs to be modified. The matrix 'Sreg' is an mxn matrix where m corresponds to the number of metabolites/components in the network and n is the number of reactions. A regulatory interaction can be added by introducing a non-zero value to Sreg(m',n') where component m' regulates reaction n'. A detailed tutorial on setting up the data structure 'Net' and the different regulatory interactions that can be introduced can be found in the original [work](http://www.seas.ucla.edu/~liaoj/download/Matlab%20module%20for%20the%20web.zip). 
