% Main Matlab script file
% Author: mh881@york.ac.uk

% Load raw data from CSV
srcData = csvread('trainingData.csv', 1);

% TODO: Solve all unique APs
O = [];
for j = find(srcData(:,1)<0)
    O = [O; srcData(j,1), srcData(j,521), srcData(j,522)];
end
APSolve(O)

% TODO: Identify set of most useful APs

% TODO: Identify set of most useful locations

% TODO: ERSGA

% TODO: Generate RSS Map

% And beyond...
