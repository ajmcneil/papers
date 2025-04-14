# vtarma

Code and data for the paper:

McNeil, A.J. (2021). Modelling volatile time series with v-transforms and copulas. Risks, 9(14). https://www.mdpi.com/2227-9091/9/1/14

### Note

Since publishing the paper I have been able to obtain even better fitting results for the VT-ARMA models by changing to a grid search for the fulcrum parameter. This is theoretically better supported because the profile likelihood can have local maxima. The scripts on this repository now give superior results to the ones published in the paper. The code to exactly reproduce the numbers in the paper is available on request.

To run the code you need the `tscopula` package which is available here: https://github.com/ajmcneil/tscopula.
