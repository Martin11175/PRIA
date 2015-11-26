% Main Matlab script file
% Author: mh881@york.ac.uk

<<<<<<< HEAD
% Load raw data from CSV
srcData = csvread('trainingData.csv', 1);

% TODO: Solve all unique APs
O = [];
for j = find(srcData(:,1)<0)
    O = [O; srcData(j,1), srcData(j,521), srcData(j,522)];
end
APSolve(O)
=======
% TODO: Load data
>>>>>>> 462b6c16f0e860d651a4c8abfb5c13b892083802

% TODO: Identify set of most useful APs

% TODO: Identify set of most useful locations

% TODO: ERSGA

% TODO: Generate RSS Map

% And beyond...