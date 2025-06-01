# hopf

To run with min. parameters: 
1) Add to path included folder ADDTOPATH
2) Load into workspace medidas/Default.mat
3) Run Auto_Hopf.m

The average SC and empirical data are saved in Data_Tagliazucchi_full.mat, which will be automatically loaded when running the Hopf model.

Output is included inside struct variable "out".
fval is different from final fitness value because solution is re-simulated to obtain full data from it (simulated tseries for example)

When you have calculated the simulation data through the Hopf model, you can use the functions in the disturbance to conduct disturbance experiments
