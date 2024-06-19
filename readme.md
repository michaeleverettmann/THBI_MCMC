This code does a joint inversion of receiver functions and surface wave data using a Markov Chain Monte Carlo inversion. You can run a synthetic joint inversion on your computer by running example_github/RUN_example.m. Once you have a synthetic inversion running, you can adapt example_github/parms/bayes_inv_parms.m to your own needs and data. 

This repository was used to run the inversion in Brunsvik et al. (2024). The code is adapted from Eilon et al. (2018). See also Golos et al. (2024), Petruska and Eilon (2022). See: Brunsvik, B., Eilon, Z., & Lynner, C. (2024). Plate‐scale imaging of eastern US reveals ancient and ongoing continental deformation. Geophysical Research Letters, 51(12), e2024GL109041.

There are dependencies for forward modelling. What you need depends on which data you will use. We are working to make these available as well in Github, but email me for any non-public codes you need.
- Mineos for surface wave phase velocity. Email me for this. As of 20240612 I am working to use the Josh Russell version of Mineos, and I will update the MCMC repo once this is done: https://github.com/jbrussell/MINEOS_synthetics.  
- Tanimoto HV code for surface wave ellipticity. We are working to make this public.  
- Propmat for receiver functions: https://github.com/brennanbrunsvik/PropMat
- We partially set up MCMC to work with telewavesim: https://github.com/paudetseis/Telewavesim. If you want to use this instead of Propmat, you will need solve a few more minor errors. There is an option to use telewavesim instead of propmat, but you should be fine to just use propmat. 
- The anisotropic HK stack repository is available if you want HK stacks. https://github.com/brennanbrunsvik/hk_anis

Some other dependencies. You may not need some of these, depending on what you are running. I placed each of these in my ~/MATLAB/ folder (see a0_STARTUP_BAYES.m). You can email me and I will send the whole set of dependencies. 
- borders. https://www.mathworks.com/matlabcentral/fileexchange/50390-borders 
- brewermap. https://github.com/DrosteEffect/BrewerMap https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps 
- fastBSpline. https://www.mathworks.com/matlabcentral/fileexchange/32509-fast-b-spline-class
- Seizmo. Use for the 1-D velocity models, and for Tau-p. https://github.com/g2e/seizmo 
- seis_tools-master. Get most recent version from Brennan. 
- EilonmyFUNCTIONS. Get most recent version from Brennan. 

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