This code does a joint inversion of receiver functions and surface wave data using a Markov Chain Monte Carlo inversion. You can run a synthetic joint inversion on your computer by running example_github/RUN_example.m. Once you have a synthetic inversion running, you can adapt example_github/parms/bayes_inv_parms.m to your own needs and data. 

This repository was used to run the inversion in Brunsvik et al. (2024). The code is adapted from Eilon et al. (2018). See also Golos et al. (2024), Petruska and Eilon (2022).  

As of 2024/06/11, I am still cleaning and documenting this repository. Expect to make adjustments to the code before it runs on your computer. 

Some steps: 
- update_bayes_inv_parms.m will change the parameters of your inversion depending on the stamp name. I recommend using this to develop any custom tests. 