% This script preps data and paths. Always run on a computer before running any inversions. 
% Before running an inversion on real station data, we need to download
% data AND make sure paths to data are correct. This script (in production,
% 2021.11.15) starts that. It should expand. 
disp('Getting station data and paths. Only need to do this once on a computer unless file is overwritten.')
run('../a0_STARTUP_BAYES.m')
run('evdata1_database.m'); 
disp("Finished getting station data and paths. ")
