function [paths] = getPaths()

% Run the automatically generated script "pahtsAutoGen"
% This creates the paths variables
% It is a structure containing many paths generated in a0_STARTUP_BAYES.m 
% This is one way of bypassing the need for a global variable (which causes
% parfor problems). 
% This approach also means we don't have to pass the variable "paths"
% through a whole ton of functions. 
pathsAutoGen; 

end