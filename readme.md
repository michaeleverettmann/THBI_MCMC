This code does a joint inversion of receiver functions and surface wave data using a Markov Chain Monte Carlo inversion. You can run a synthetic joint inversion on your computer by running example_github/RUN_example.m. Once you have a synthetic inversion running, you can adapt example_github/parms/bayes_inv_parms.m to your own needs and data. 

This repository was used to run the inversion in Brunsvik et al. (2024). The code is adapted from Eilon et al. (2018). See also Golos et al. (2024), Petruska and Eilon (2022). See: Brunsvik, B., Eilon, Z., & Lynner, C. (2024). Plate‐scale imaging of eastern US reveals ancient and ongoing continental deformation. Geophysical Research Letters, 51(12), e2024GL109041.

There are dependencies for forward modelling. What you need depends on which data you will use. We are working to make these available as well in Github, but for now, you can contact me for the codes. 
- Mineos for surface wave phase velocity. As of 20240612 I am working to use the Josh Russell version of Mineos.  
- Propmat for receiver functions. 
- Tanimoto HV code for surface wave ellipticity. 
- We partially set up MCMC to work with telewavesim. There are still probably some minor errors. There is an option to use telewavesim instead of propmat, but you should be fine to just use propmat. 
- The anisotropic HK stack repository is available if you want HK stacks. https://github.com/brennanbrunsvik/hk_anis

Some other dependencies. Look through a0_STARTUP_BAYES.m to understand these.  
- fastBSpline
- Seizmo. Velocity models maybe needed? Maybe also needed for taup. 
- seis_tools-master. Get most recent version from Brennan. 
- EilonmyFUNCTIONS. Get most recent version from Brennan. 
- Some colormaps from Matlab file exchange. 

Some steps: 
- Change your paths in a0_STARTUP_BAYES.m. Obtain any dependencies. 
- To speed up splines, go to your fastBSpline folder. Run the compile mex code. 
- Run the synthetic test in example_github. 
- For working with real data, modify "matguts" to load data from your computer. The "ENAM" folder was used for Brunsvik et al. 2024. Use it as a starting place. 
-- There is a step for downloading station information. I don't remember it, but I think it had to do with evdata1_database.m. Add this information to the readme if you find out. 
- The "transfer" folder migrates code from your work machine to HPC. 
- Modify the "ENAM/batch" code to work with your HPC. This has Slurm code for submitting jobs to Linux clusters at UCSB. 
- retrieve_results_example shows how I get results back from the HPC and do some formatting. I ran this from an external drive. Should have ~100 Mb/station.  
- 3d_models takes results from the many stations and does a smoothed inversion for a 3D model that matches the stations model PDFs. 

Some behaviors to expect: 
- In order to speed up IO on systems with many cores (e.g. 40 cores), we read and write to a folder in the ram directly for calling some fortran codes. A folder called "ramDriveTHBI" should mount automatically. 
- update_bayes_inv_parms.m will change the parameters of your inversion depending on the stamp name. I recommend using this to develop any custom tests. 
- run_one_chain.m is obnoxiously complicated because it needs to deal with Mineos failing. Do much testing if you modify this. 