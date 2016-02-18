function [ G ] = RGEA( J, K )
%RGEA EZ Relative Gain Estimation Algorithm
%   J [in] - Matrix of RSSI measurements (APs as columns)
%   K [in] - Vector of device IDs relating rows to device
%   G [out] - Pairs of device IDs (1) to estimated gain (2)

prox_threshold = 3;
min_AP_overlap = 2;

D = sort(unique(K)); % List of device IDs
J(J == 100) = -100; % Replace positive invisibility markers to prevent skew

% Calculate relative gain between pairs of devices
deltaG = zeros(size(D,1));
sigma_deltaG = zeros(size(D,1)); % Uncertainty (estimated standard deviation)
% TODO: parfor to increase performance?
for i = 1:(size(D,1) - 1)
    for j = (i+1):size(D,1)
        
        % Test pairs of measurements for proximity by similarity
        avg_diff = zeros(sum(K == D(i)), sum(K == D(j)));
        
        % TODO: Difference matrix between vectors = bsxfun(@minus, x, x')
        count_m = 1;
        for m = J(K == D(i), :)'
            count_n = 1;
            for n = J(K == D(j), :)'
                % Only compare APs visible at at least one location
                vis = (m > -100) | (n > -100);
                if(sum(vis) > min_AP_overlap && mean(abs(m(vis) - n(vis))) < prox_threshold)
                    avg_diff(count_m, count_n) = mean(m(vis) - n(vis));
                    count_n = count_n + 1;
                end
            end
            count_m = count_m + 1;
        end
        
        % Compare all proximate locations
        prox = abs(avg_diff) < prox_threshold & (avg_diff ~= 0);
        
        if sum(sum(prox)) > 10
            deltaG(i, j) = mean(avg_diff(prox));
            deltaG(j, i) = -deltaG(i, j);
            sigma_deltaG(i, j) = (1 / sum(sum(prox))) * sqrt(sum((avg_diff(prox) - deltaG(i,j)).^2));
            sigma_deltaG(j, i) = sigma_deltaG(i, j);
        end
        fprintf('%f | %d:%d / %d\n', deltaG(i,j), i, j, size(D,1))
    end
end

% TODO: Determine sets of devices that can be estimated relative to each other
% Currently assuming all devices can be see at least one proximate location

% Objective function aiming to determine individual gain values from estimated relative gain
objective = @(gain) GainSolve(gain, deltaG, sigma_deltaG);

% Gain should be between +/- 20 dB
lower = zeros(size(D,1), 1) - 20;
upper = zeros(size(D,1), 1) + 20;

% Perform simulated annealing, starting from average values
G0 = zeros(size(D,1), 1);
options = optimset('display','off');
G = simulannealbnd(objective, G0, lower, upper, options);

G = [D G];

end

