This code does a joint inversion of receiver functions and surface wave data using a Markov Chain Monte Carlo inversion. You can run a synthetic joint inversion on your computer by running example_github/RUN_example.m. Once you have a synthetic inversion running, you can adapt example_github/parms/bayes_inv_parms.m to your own needs and data. 

This repository was used to run the inversion in Brunsvik et al. (2024). The code is adapted from Eilon et al. (2018). See also Golos et al. (2024), Petruska and Eilon (2022).  

As of 2024/06/11, I am still cleaning and documenting this repository. Expect to make adjustments to the code before it runs on your computer. 

Some dependencies. Look through a0_STARTUP_BAYES.m to understand these.  
- fastBSpline
- Seizmo. Velocity models maybe needed? Maybe also needed for taup. 
- seis_tools-master. Get most recent version from Brennan. 
- EilonmyFUNCTIONS. Get most recent version from Brennan. 
- Some colormaps from Matlab file exchange. 

Dependencies for forward modelling. What you need depends on which data you will use. 
- Mineos for surface wave phase velocity. As of 20240612 I am working to use the Josh Russell version of Mineos.  
- Propmat for receiver functions. 
- Tanimoto HV code for surface wave ellipticity. 
- We partially set up MCMC to work with telewavesim. There are still probably some minor errors. There is an option to use telewavesim instead of propmat, but you should be fine to just use propmat. 
- The anisotropic HK stack repository if you want HK stacks. https://github.com/brennanbrunsvik/hk_anis

Some steps: 
- update_bayes_inv_parms.m will change the parameters of your inversion depending on the stamp name. I recommend using this to develop any custom tests. 
- To speed up splines, go to your fastBSpline folder. Run the compile mex code. 