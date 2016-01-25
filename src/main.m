% Main Matlab script file
% Author: mh881@york.ac.uk

% Load raw data from CSV
srcData = csvread('trainingData.csv', 1);

floor = 1;
building = 1;

% Isolate a 2D sub-space (single floor of single building)
rows = ((srcData(:,523) == floor) & (srcData(:,524) == building));
dataSet = srcData(rows, :);

% TODO: Identify set of most useful locations

% TODO: Identify set of most useful APs

% For each access point i
C = zeros(520, 4);
error = zeros(520,1);
for i = 1:520
    % Find which observations include our desired AP
    %J = ((-75 < dataSet(:,i)) & (dataSet(:,i) < 0));
    J = (dataSet(:,i) < 0);
    O = dataSet(J,[i 521 522]);
    
    if size(O,1) > 5
        % Parameterised objective function looking to minimise error with respect
        % to observed RSSI
        objective = @(C) APSolve(C, O);
        
        % Lower and upper bounds for parameter search:
        %   Lat and long (in metres): 50m outside of observation locations
        %   Transmit power: 0 to -50dBm
        %   Path loss: 1.5 to 6.0 from EZ's trials
        % TODO: Is path loss rate from papers comparable? (i.e. from Meters to lat/long)
        lower = [min(O(:,2)) - 50, min(O(:,3)) - 50, -50, 1.5];
        upper = [max(O(:,2)) + 50, max(O(:,3)) + 50, 0, 6.0];

        % Perform simulated annealing, starting from average values
        C0 = [mean(O(:,2)), mean(O(:,3)), -25, 3.0];
        options = optimset('display','off');
        [C(i,:), error(i)] = simulannealbnd(objective, C0, lower, upper, options);
    end
end

% Play 2 seconds of 50Hz sound on completion
sound(sin(linspace(0, 2*50*2*pi, round(2*1000))), 1000)
display(mean(error(error~=0)))

% TODO: ERSGA

% TODO: Generate RSS Map

% And beyond...
