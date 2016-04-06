% Relative Gain Estimation Test
% Author: mh881@york.ac.uk

% Load raw data from CSV
srcData = csvread('UJIndoorLoc/trainingData.csv', 1);

% Scale known location data to a local space
min_lat = min(srcData(:,521));
min_long = min(srcData(:,522));
srcData(:,521) = srcData(:,521) - min_lat;
srcData(:,522) = srcData(:,522) - min_long;

% Prepare summary files
if exist('ez_gain_summary.csv', 'file') == 0
    fprintf(fopen('ez_gain_summary.csv', 'w'), ...
        'floor,building,prox_threshold,min_AP_overlap,gain_error\n');
end
if exist('h_gain_summary.csv', 'file') == 0
    fprintf(fopen('h_gain_summary.csv', 'w'), ...
        'floor,building,num_strongest_APs,min_AP_strength,min_strong_overlap,gain_error\n');
end

% Parameters to test
% EZ
prox_threshold = [ 2 3 4 5 7 10 ];
min_AP_overlap = [ 1 2 3 4 5 ];
% Heuristic
min_AP_strength = [ -100 -90 -80 -70 ];
num_strongest = [ 2 3 4 5 7 10 ];
min_strong_overlap = [ 2 3 4 5 6 ];

% Test total accuracy -----------------------------------------------------
J = srcData(:, 1:520);
K = srcData(:, 528);
D = sort(unique(K));
X = srcData(:, 521:523);

% Calculate base line estimates
g = GroundRGEA(J,K,X);

% Run algortihms and write summaries
for prox = prox_threshold
    for overlap = min_AP_overlap
        if prox == 2 || prox == 3 || prox == 4 || ...
            ((prox == 5) && (overlap == 1 || overlap == 2 || overlap == 3))
            continue;
        end
        fprintf('Complete | EZ: %d %d\n', prox, overlap);
        ez = RGEA(J,K,prox,overlap);
        
        dlmwrite('ez_gain_summary.csv', ...
            [-1, -1, prox, overlap, sum(abs(g(:,2) - ez(:,2)))], ...
            'delimiter', ',', '-append');
    end
end

for strength = min_AP_strength
    for num_str = num_strongest
        for str_overlap = min_strong_overlap
            fprintf('Complete | H: %d %d %d\n', num_str, strength, str_overlap);
            h = SimpleRGEA(J,K, strength, num_str, str_overlap);

            dlmwrite('h_gain_summary.csv', ...
                [-1, -1, num_str, strength, str_overlap, sum(abs(g(:,2) - h(:,2)))], ...
                'delimiter', ',', '-append');
        end
    end
end


% Test building-wide ------------------------------------------------------
for building = [ 0 1 2 ]
    
    dataSet = (srcData(:,524) == building);
    
    % Set data up for test
    J = srcData(dataSet, 1:520);
    K = srcData(dataSet, 528);
    D = sort(unique(K));
    X = srcData(dataSet, 521:523);
    
    % Calculate base line estimates
    g = GroundRGEA(J,K,X);
    
    % Run algortihms and write summaries
    for prox = prox_threshold
        for overlap = min_AP_overlap
            fprintf('Building %d | EZ: %d %d\n', building, prox, overlap);
            ez = RGEA(J,K,prox,overlap);
            
            dlmwrite('ez_gain_summary.csv', ...
                [-1, -1, prox, overlap, sum(abs(g(:,2) - ez(:,2)))], ...
                'delimiter', ',', '-append');
        end
    end
    
    for strength = min_AP_strength
        for num_str = num_strongest
            for str_overlap = min_strong_overlap
                fprintf('Building %d | H: %d %d %d\n', building, num_str, strength, str_overlap);
                h = SimpleRGEA(J,K, strength, num_str, str_overlap);
                
                dlmwrite('h_gain_summary.csv', ...
                    [-1, -1, num_str, strength, str_overlap, sum(abs(g(:,2) - h(:,2)))], ...
                    'delimiter', ',', '-append');
            end
        end
    end
end

% Test individual floor accuracy ------------------------------------------
for building = [ 0 1 2 ]
for floor = [ 0 1 2 3 ]
    
    dataSet = ((srcData(:,523) == floor) & (srcData(:,524) == building));
    
    % Set data up for test
    J = srcData(dataSet, 1:520);
    K = srcData(dataSet, 528);
    D = sort(unique(K));
    X = srcData(dataSet, 521:523);
    
    % Calculate base line estimates
    g = GroundRGEA(J,K,X);
    
    % Run algortihms and write summaries
    for prox = prox_threshold
        for overlap = min_AP_overlap
            fprintf('F(%d) B(%d) | EZ: %d %d\n', building, floor, prox, overlap);
            ez = RGEA(J,K,prox,overlap);
            
            dlmwrite('ez_gain_summary.csv', ...
                [-1, -1, prox, overlap, sum(abs(g(:,2) - ez(:,2)))], ...
                'delimiter', ',', '-append');
        end
    end
    
    for strength = min_AP_strength
        for num_str = num_strongest
            for str_overlap = min_strong_overlap
                fprintf('F(%d) B(%d) | H: %d %d %d\n', building, floor, num_str, strength, str_overlap);
                h = SimpleRGEA(J,K, strength, num_str, str_overlap);
                
                dlmwrite('h_gain_summary.csv', ...
                    [-1, -1, num_str, strength, str_overlap, sum(abs(g(:,2) - h(:,2)))], ...
                    'delimiter', ',', '-append');
            end
        end
    end

end
end