
???: nianqiao@purdue.edu

This file explains the MCEM experiments on the linear regression example. 
The experimental setups are similar to that of the RJMCMC experiments.

All scripts are contained in the file **LRMCEM\_unbounded\_helper.R**. The main function is **run\_mcem**, which contains the confidential data process, private query generation process, and implementation of the MCEM algorithm. 


We run our experiments under r/4.3.1. The required packages are **VGAM** and **MASS**.

Our experiments concerns different privacy loss budget on sample size, with $\epsilon_N=0.001,0.01,0.1,1,10,\inf.$
The other parameters are
1. Number of MCEM cycles = 10000;
2. N = 1000
3. privacy loss budget on s = 1;
4. true data generating parameters

```
params_true = list(beta = c(0,-1,1),
                   tau = 1,
                   mu = c(1,-1),
                   phi = diag(p),
                   phi_inv = diag(p))
```

5. hyper parameters for the priors
```
hyperparams = list(p=2,
                   m = rep(0,p+1),
                   V_inv = diag(p+1),
                   a = 2,
                   b = 2,
                   theta = rep(0,p),
                   Sigma_inv = diag(p),
                   d=2,## df needs to be >=p
                   W_inv = diag(p))
```
6. privacy parameters 
```
privacy_params = list(eps=eps,
                      ep_N = epn,
                      deltaa = p^2/2+(2.5)*p+2,  ### bounded should have been p^2+4p+3, not p^2+3p+3...
                      clbounds=c(-5,-5,5,5))
```


We have saved simulation results in the file **summarytable.txt** and we produced the table in the manuscript with script **LRMCEM_analyze.R**. The analysis file uses **dplyr** and **reshape**.
