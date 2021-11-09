function [paths] = getPaths()
%     arguments 
%         options.pathToMain = '..' 
%     end
%     
% paths = load([options.pathToMain '/misc/paths.mat']); % TODO not sure how to figure out which directory is the THBI directory in the easiest way. 
% paths = paths.paths; 

% Run the automatically generated script "pahtsAutoGen"
% This creates the paths variables
% It is a structure containing many paths generated in a0_STARTUP_BAYES.m 
% This is one way of bypassing the need for a global variable, and for not
% passing "paths" through dozens of functions. 
pathsAutoGen; 

end