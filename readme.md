# Joint Inversion of Receiver Functions and Surface Waves Using Markov Chain Monte Carlo

This code performs a joint inversion of receiver functions and surface wave data using a Markov Chain Monte Carlo (MCMC) inversion. This is a Trans-dimensional, Hierarchical, Bayesian Inversion (THBI). 

This code accompanies:
- Brunsvik, B., Eilon, Z., & Lynner, C. (2024). Plate‐scale imaging of eastern US reveals ancient and ongoing continental deformation. *Geophysical Research Letters, 51*(12), e2024GL109041.
- The code was modified from Eilon et al. (2018). See also Golos et al. (2024) and Petruska and Eilon (2022).

## Running the example synthetic joint inversion

1. Run the example synthetic joint inversion by executing `example_github/RUN_example.m`.
2. Once the synthetic inversion is running, you can modify the eastern US example, applied to real data, found in the `ENAM` folder.

## Dependencies

### Required depending on your data types

- **Surface wave phase velocity**: Mineos (contact author for access). Working on integration with Josh Russell's public version: [Mineos_synthetics](https://github.com/jbrussell/MINEOS_synthetics).
- **Surface wave ellipticity**: HV code from Tanimoto and Tsuboi (2009, GRL). Public release pending.
- **Receiver functions**: [PropMat](https://github.com/brennanbrunsvik/PropMat). Alternatively, [Telewavesim](https://github.com/paudetseis/Telewavesim) is partially supported.
- **$H-\kappa$ stacks**, including anisotropy: [hk_anis](https://github.com/brennanbrunsvik/hk_anis).

### Other dependencies

These were placed in the `~/MATLAB/` folder. Refer to `a0_STARTUP_BAYES.m` for setup details.

- [borders](https://www.mathworks.com/matlabcentral/fileexchange/50390-borders)
- [brewermap](https://github.com/DrosteEffect/BrewerMap)
- [fastBSpline](https://www.mathworks.com/matlabcentral/fileexchange/32509-fast-b-spline-class)
- [Seizmo](https://github.com/g2e/seizmo) 
- [seis_tools](https://github.com/eilonzach/seis_tools) 
- [myFUNCTIONS](https://github.com/eilonzach/myFUNCTIONS) 

## Setup steps

1. Update paths in `a0_STARTUP_BAYES.m` and obtain necessary dependencies.
1. Compile codes. 
    1. All Fortran forward modelling codes. 
    1. mex code in the `fastBSpline` folder to speed up spline calculations (optional).
1. Run the synthetic test in the `example_github` folder.
1. For real data applications, use the `ENAM` folder as a starting point. Modify "matguts" to load your data from your computer. 
1. For station information, refer to `evdata1_database.m` (add details to the readme if found).
1. Use the "transfer" folder to migrate code from your work machine to HPC. 
1. Modify `ENAM/batch` code (Slurm scripts for UCSB clusters) for your HPC setup.
1. Use `retrieve_results_example` to get results from the HPC and format them (expect ~100 Mb/station).
1. The `3d_models` folder helps create a smoothed 3D model from the station model PDFs.

## Expected behaviors

- IO optimization: For systems with many cores, IO operations need a folder in RAM ("ramDriveTHBI") for efficiency. This should mount automatically.  
- Parameter updates: `update_bayes_inv_parms.m` will modify your inversion parameters depending on your STAMP. 
- Running chains: `run_one_chain.m` is complex because it handles many cases where Mineos fails. Test thoroughly if modified.

For any further assistance or to obtain dependencies, feel free to contact the author.