% Main Matlab script file
% Author: mh881@york.ac.uk

%----------------------------Parameter Setup------------------------------%

% Load raw data from CSV
% Algorithms generally set for 520 APs, should accept 100 (positive
% invisibility marker) for unused access points as a workaround for smaller
% data sets.
if exist('srcData', 'var') == 0
    srcData = csvread('UJIndoorLoc/trainingData.csv', 1);
end
if exist('testData', 'var') == 0
    testData = csvread('UJIndoorLoc/validationData.csv', 1);
end
if exist('chance', 'var') == 0
    chance = 100;
end

% Algorithm types
if exist('RGEA_type', 'var') == 0
    RGEA_type = 'ground';
end
if exist('apselect_type', 'var') == 0
    apselect_type = 'overlap';
end
if exist('locselect_type', 'var') == 0
    locselect_type = 'none';
end
if exist('ips_type', 'var') == 0
    ips_type = 'all';
end

% Scope definition
if exist('bounds', 'var') == 0
    bounds = 50;
end
if exist('floors', 'var') == 0
    floors = 0;
end
if exist('buildings', 'var') == 0
    buildings = 0;
end
if exist('thresholds', 'var') == 0
    thresholds = -100;
end

% Scale known location data to a local space
min_lat = min(srcData(:,521));
min_long = min(srcData(:,522));
srcData(:,521) = srcData(:,521) - (min_lat - bounds);
srcData(:,522) = srcData(:,522) - (min_long - bounds);
testData(:,521) = testData(:,521) - (min_lat - bounds);
testData(:,522) = testData(:,522) - (min_long - bounds);

%-----------------------------EZ-ALGORITHM--------------------------------%

% Relative Gain Estimation Algorithm
if strcmp(RGEA_type,'none') == 0
    tic
    if strcmp(RGEA_type,'rgea') == 1
        G = RGEA(srcData(:, 1:520), srcData(:, 528));
    elseif strcmp(RGEA_type,'simple') == 1
        G = SimpleRGEA(srcData(:, 1:520), srcData(:, 528));
    else
        G = GroundRGEA(srcData(:, 1:520), srcData(:, 528), srcData(:,521:523));
    end
    for d = G'
        deviceMeas = srcData(srcData(:,528) == d(1), 1:520);
        visMeas = (deviceMeas(:,1:520) ~= 100) & (deviceMeas(:,1:520) ~= -100);
        deviceMeas(visMeas) = deviceMeas(visMeas) - d(2);
    end
    toc
end

% Iterate over 2D subspaces
for floor = floors
for building = buildings
for threshold = thresholds

% Isolate a 2D sub-space (single floor of single building)
rows = ((srcData(:,523) == floor) & (srcData(:,524) == building));
dataSet = srcData(rows, :);
normalisedData = 1 - (abs(dataSet(:, 1:520)) / 100);
original_locations = dataSet(:, 521:522);

% Artificial GPS restriction
for j = 1:size(dataSet,1)
    if randi(100) > chance
        dataSet(j, 521:522) = [0 0];
    end
end

% APSelect Algorithm
if strcmp(apselect_type,'none') == 0
    tic
    ap_rank = zeros(520, 1);
    for i = 1:520 % Rank by number of visible known locations
        ap_rank(i) = sum(dataSet(normalisedData(:,i) > 0, 521) > 0);
    end
    if strcmp(apselect_type,'all') == 1
        APSelect = HeirarchicalCluster(normalisedData, ap_rank, 'all');
    elseif strcmp(apselect_type,'max') == 1
        APSelect = HeirarchicalCluster(normalisedData, ap_rank, 'max');
    else
        APSelect = HeirarchicalCluster(normalisedData, ap_rank, 'overlap');
    end
    toc
else
    APSelect = 1:520;
end

% LocSelect Algorithm
if strcmp(locselect_type,'none') == 0
    tic
    % Prefer GPS localised measurements (binary rank)
    loc_rank = (dataSet(:,521) ~= 0 | dataSet(:,522) ~= 0);
    if strcmp(locselect_type,'all') == 1
        LocSelect = HeirarchicalCluster(normalisedData', loc_rank, 'all');
    elseif strcmp(locselect_type,'max') == 1
        LocSelect = HeirarchicalCluster(normalisedData', loc_rank, 'max');
    else
        LocSelect = HeirarchicalCluster(normalisedData', loc_rank, 'overlap');
    end
    toc
else
    LocSelect = 1:size(normalisedData,2);
end
LocSelect = ismember(1:size(normalisedData,2), LocSelect);

% Iterate until no new parameters can be determined
APparams = zeros(520, 4);
new_value_flag = true;
while new_value_flag
    new_value_flag = false;

    % Determine new AP parameters from localised measurements
    for i = APSelect(APparams(APSelect, 4) == 0)
        % Find which observations include our desired AP
        J = ((threshold < dataSet(:,i)) & (dataSet(:,i) < 0) ...
            & ((dataSet(:,521) ~= 0) | (dataSet(:,522) ~= 0)) ...
            & LocSelect);
        O = dataSet(J,[i 521 522]);
        
        if size(O,1) > 4
            new_value_flag = true;
            % Objective function looking to minimise RSSI error
            objective = @(C) APSolve(C, O);
            
            % Lower and upper bounds for parameter search:
            %   Lat and long (in metres): bounds outside of observation locations
            %   Transmit power: 0 to -50dBm
            %   Path loss: 1.5 to 6.0 from EZ's trials
            % TODO: In report mention that most access points were hitting the 1.5 bound
            lower = [min(O(:,2)) - bounds, min(O(:,3)) - bounds, -50, 1.5];
            upper = [max(O(:,2)) + bounds, max(O(:,3)) + bounds, 0, 6.0];
            
            % Perform simulated annealing, starting from average values
            C0 = [mean(O(:,2)), mean(O(:,3)), -25, 3.0];
            options = optimset('display','off');
            fprintf('AP: %d, %d observations...', i, size(O,1));
            tic
            APparams(i,:) = simulannealbnd(objective, C0, lower, upper, options);
            toc
            fprintf('Median localisation error: %fm\n', ...
                median(abs((sqrt((O(:,2) - APparams(i,1)).^2 + (O(:,3) - APparams(i,2)).^2) ...
                - 10.^((APparams(i,3) - O(:,1))./(10*APparams(i,4)))))));
        end
    end
    
    % If new APs are parameterised, check to localise new measurements
    if new_value_flag
        new_value_flag = false;
        
        for j = dataSet((dataSet(:,521) == 0) & (dataSet(:,522) == 0) & LocSelect)
            % Find which APs we wish to evaluate
            I = ((threshold < j(1:520)) & (j(1:520) < 0) ...
                & (APparams(:,4) ~= 0));
            if strcmp(ips_type,'max') == 1
                [~, order] = sort(j(1:520), 'descend');
                I = I & ismember(1:520, order(1:4));
            end
            O = j(I);
            C = APparams(I, :);
            
            if size(O,1) > 2
                new_value_flag = true;
                % Calculate distance
                D = 10.^((C(:,3) - O)./(10*C(:,4)));
                % Fix parameters to search over
                objective = @(J) LocSolve(J, D, C);
                
                % Bounds on area to search (outside AP locations)
                lower = [min(C(:,1)) - bounds, min(C(:,1)) - bounds];
                upper = [max(C(:,2)) + bounds, max(C(:,2)) + bounds];
                
                % Perform simulated annealing, starting from average values
                J0 = [mean(C(:,1)), mean(C(:,2))];
                options = optimset('display','off');
                dataSet(j,521:522) = simulannealbnd(objective, J0, lower, upper, options);
            end
        end
    end
end

%------------------------------EVALUATION---------------------------------%

% Isolate a 2D sub-space (single floor of single building)
rows = ((testData(:,523) == floor) & (testData(:,524) == building));
test_dataSet = testData(rows, :);

% Localise validation data
IPSresults = zeros(size(test_dataSet,1),2);
IPSerror = zeros(size(test_dataSet,1),1);
n = 1;
for j = test_dataSet'
    % Find which APs we wish to evaluate
    I = ((threshold < j(1:520)) & (j(1:520) < 0) & (APparams(:,4) ~= 0));
    if strcmp(ips_type,'max') == 1
        [~, order] = sort(j(1:520), 'descend');
        I = I & ismember(1:520, order(1:4));
    end
    O = j(I);
    C = APparams(I, :);
        
    if size(O,1) > 2
        D = 10.^((C(:,3) - O)./(10*C(:,4)));
        % Fix parameters to search over
        objective = @(J) LocSolve(J, D, C);
        
        % Bounds on area to search (outside AP locations)
        lower = [min(C(:,1)) - bounds, min(C(:,1)) - bounds];
        upper = [max(C(:,2)) + bounds, max(C(:,2)) + bounds];
    
        % Perform simulated annealing, starting from average values
        J0 = [mean(C(:,1)), mean(C(:,2))];
        options = optimset('display','off');
        IPSresults(n,:) = simulannealbnd(objective, J0, lower, upper, options);
        
        % Calculate localisation error
        IPSerror(n) = sqrt((IPSresults(n,1) - dataSet(n,521))^2 + (IPSresults(n,2) - dataSet(n,522))^2);
    end
    n = n + 1;
end

fprintf('%d: %d/%d localised. %f median IPS error\n', ...
    threshold, sum(IPSerror > 0), size(IPSerror,1), median(IPSerror(IPSerror > 0)))

% Create results directory structure
if exist('results', 'dir') == 0
    mkdir('results');
end
building_dir = sprintf('results/building_%d', building);
if exist(building_dir, 'dir') == 0
    mkdir(building_dir);
end
floor_dir = sprintf('results/building_%d/floor_%d', building, floor);
if exist(floor_dir, 'dir') == 0
    mkdir(floor_dir);
end

% Output csv files of results
src_filename = sprintf('results/building_%d/floor_%d/src_rgea:%s_ap:%s_loc:%s_ips:%s_threshold:%d_bounds:%d_chance:%d.csv', ...
    building, floor, RGEA_type, apselect_type, locselect_type, ips_type, threshold, bounds, chance);
csvwrite(src_filename, [original_locations dataSet(:, 521:522) LocSelect]);

tst_filename = sprintf('results/building_%d/floor_%d/tst_rgea:%s_ap:%s_loc:%s_ips:%s_threshold:%d_bounds_chance:%d:%d.csv', ...
    building, floor, RGEA_type, apselect_type, locselect_type, ips_type, threshold, bounds, chance);
csvwrite(tst_filename, [test_dataSet(:,521:522) IPSresults]);

ap_filename = sprintf('results/building_%d/floor_%d/ap_rgea:%s_ap:%s_loc:%s_ips:%s_threshold:%d_bounds_chance:%d:%d.csv', ...
    building, floor, RGEA_type, apselect_type, locselect_type, ips_type, threshold, bounds, chance);
csvwrite(ap_filename, APparams);
end
end
end
