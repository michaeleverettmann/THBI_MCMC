function [paths] = getPaths() 
paths = load('../misc/paths.mat'); % TODO not sure how to figure out which directory is the THBI directory in the easiest way. 
paths = paths.paths; 
end