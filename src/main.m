% Main Matlab script file
% Author: mh881@york.ac.uk

% Load raw data from CSV
srcData = csvread('UJIndoorLoc/trainingData.csv', 1);
testData = csvread('UJIndoorLoc/validationData.csv', 1);

% Scale known location data to a local space
min_lat = min(srcData(:,521));
min_long = min(srcData(:,522));
srcData(:,521) = srcData(:,521) - (min_lat - bounds);
srcData(:,522) = srcData(:,522) - (min_long - bounds);
testData(:,521) = testData(:,521) - (min_lat - bounds);
testData(:,522) = testData(:,522) - (min_long - bounds);

%-----------------------------EZ-ALGORITHM--------------------------------%

% Relative Gain Estimation Algorithm
G = GroundRGEA(srcData(:, 1:520), srcData(:, 528), srcData(:,521:523));
for d = G'
    deviceMeas = srcData(srcData(:,528) == d(1), 1:520);
    visMeas = (deviceMeas(:,1:520) ~= 100) & (deviceMeas(:,1:520) ~= -100);
    deviceMeas(visMeas) = deviceMeas(visMeas) - d(2);
end

% Iterate over 2D subspaces
floor = 0;
building = 0;
threshold = -70;
bounds = 50;

%for threshold = [-100 -90 -80 -75 -70 -65 -60 -50]

% Isolate a 2D sub-space (single floor of single building)
rows = ((srcData(:,523) == floor) & (srcData(:,524) == building));
dataSet = srcData(rows, :);

% APSelect Algorithm
normalisedData = 1 - (abs(dataSet(:, 1:520)) / 100);
rank = zeros(520, 1);
for i = 1:520 % Rank by number of visible known locations
    % TODO: Omit visible places with threshold too?
    rank(i) = size(find(dataSet(normalisedData(:,i) > 0, 521) > 0), 1);
end
APSelect = HeirarchicalCluster(normalisedData, rank);

% LocSelect Algorithm
% TODO: Test efficacy of LocSelect
%LocSelect = HeirarchicalCluster(normalisedData', zeros(size(normalisedData',1),1));

% For each access point i
APparams = zeros(520, 4);
for i = APSelect
    % Find which observations include our desired AP
    J = ((threshold < dataSet(:,i)) & (dataSet(:,i) < 0));
    O = dataSet(J,[i 521 522]);
    
    if size(O,1) > 4
        % Parameterised objective function looking to minimise error with respect
        % to observed RSSI
        objective = @(C) APSolve(C, O);
        
        % Lower and upper bounds for parameter search:
        %   Lat and long (in metres): boundsm outside of observation locations
        %   Transmit power: 0 to -boundsdBm
        %   Path loss: 1.5 to 6.0 from EZ's trials
        % TODO: In report mention that most access points were hitting the 1.5 bound
        lower = [min(O(:,2)) - bounds, min(O(:,3)) - bounds, -bounds, 1.5];
        upper = [max(O(:,2)) + bounds, max(O(:,3)) + bounds, 0, 6.0];

        % Perform simulated annealing, starting from average values
        C0 = [mean(O(:,2)), mean(O(:,3)), -25, 3.0];
        options = optimset('display','off');
        APparams(i,:) = simulannealbnd(objective, C0, lower, upper, options);
    end
end

% TODO: ERSGA

%------------------------------EVALUATION---------------------------------%

% Isolate a 2D sub-space (single floor of single building)
rows = ((testData(:,523) == floor) & (testData(:,524) == building));
dataSet = testData(rows, :);

% Localise validation data
% TODO: Can we change this to least-squares constraint solving of the form
% C.x = d? (See RGEA solvers)
IPSresults = zeros(size(dataSet,1),2);
RSSIerror = zeros(size(dataSet,1),1);
IPSerror = zeros(size(dataSet,1),1);
n = 1;
for j = dataSet'
    % Find which APs we wish to evaluate
    I = ((threshold < j(1:520)) & (j(1:520) < 0));
    O = j(I);
    C = APparams(I, :);
        
    if size(O,1) > 2
        % Fix parameters to search over
        objective = @(J) LocSolve(J, O, C);
        
        % Bounds on area to search (boundsm outside AP locations)
        lower = [min(C(:,1)) - bounds, min(C(:,1)) - bounds];
        upper = [max(C(:,2)) + bounds, max(C(:,2)) + bounds];
    
        % Perform simulated annealing, starting from average values
        J0 = [mean(C(:,1)), mean(C(:,2))];
        options = optimset('display','off');
        [IPSresults(n,:), RSSIerror(n)] = simulannealbnd(objective, J0, lower, upper, options);
        
        % Calculate localisation error
        IPSerror(n) = sqrt((IPSresults(n,1) - dataSet(n,521))^2 + (IPSresults(n,2) - dataSet(n,522))^2);
    end
    n = n + 1;
end

fprintf('%d: %d/%d localised. %f median IPS error, %f median RSSI error \n', threshold, size(find(IPSerror),1), size(IPSerror,1), median(IPSerror(IPSerror > 0)), median(RSSIerror(RSSIerror > 0)))

%end
beep;
% TODO: Average RSS Map (remember -100 if AP not visible, ignore unclassified APs)
