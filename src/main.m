% Main Matlab script file
% Author: mh881@york.ac.uk

% Load raw data from CSV
srcData = csvread('trainingData.csv', 1);

% Progression sound
blip = 0.3*sin(linspace(0, 0.1*430*2*pi, round(0.1*1000)));

floor = 1;
building = 1;
threshold = -75;

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
    J = ((threshold < dataSet(:,i)) & (dataSet(:,i) < 0));
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
    
    % Play progression sound every 100 APs
    if mod(i,100) == 0
       sound(blip, 1000);
       fprintf('i = %d\n',i);
    end
end

% Play sound on completion
beep;
display(mean(error(error~=0)))

% TODO: ERSGA

% TODO: Generate RSS Map

% And beyond...
