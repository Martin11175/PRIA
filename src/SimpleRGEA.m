function [ G ] = SimpleRGEA( J, K, X )
%RGEA Simplified Relative Gain Estimation Algorithm using location data
%   J [in] - Matrix of RSSI measurements (APs as columns)
%   K [in] - Vector of device IDs relating rows to device
%   X [in] - Matrix of locations for RSSI measurements
%       [latitude, longitude, floor]
%   G [out] - Estimated gain by device ID

D = sort(unique(K)); % List of device IDs
J(J == 100) = -100; % Replace positive invisibility markers to prevent skew

% Calculate relative gain between pairs of devices
deltaG = zeros(size(D,1));
sigma_deltaG = zeros(size(D,1)); % Uncertainty (estimated standard deviation)
% TODO: Replace i and j with nchoosek
for i = 1:(size(D,1) - 1)
    k_1 = J(K == D(i), :);
    x_1 = X(K == D(i), :);
    for j = (i+1):size(D,1)
        k_2 = J(K == D(j), :);
        x_2 = X(K == D(j), :);
        
        % Test pairs of measurements for proximity by similarity
        avg_diff = zeros(size(k_1, 1), size(k_2, 1));
        
        % Proximate on the same floor within 1m
        prox = (abs(bsxfun(@minus, x_1(:,1), x_2(:,1)')) < 1 ...
                & abs(bsxfun(@minus, x_1(:,2), x_2(:,2)')) < 1 ...
                & bsxfun(@eq, x_1(:,3), x_2(:,3)'));
        [m, n] = find(prox);
        
        if size(m, 1) > 0
            for p = [m n]'
                % Only compare visible APs
                vis = (k_1(p(1),:) > -100) & (k_2(p(2),:) > -100);
                if sum(vis) > 0
                    avg_diff(p(1), p(2)) = mean(k_1(p(1), vis) - k_2(p(2), vis));
                end
            end
        
            deltaG(i, j) = mean(avg_diff(prox));
            deltaG(j, i) = -deltaG(i, j);
            sigma_deltaG(i, j) = (1 / size(m, 1)) * sqrt(sum((avg_diff(prox) - deltaG(i,j)).^2));
            sigma_deltaG(j, i) = sigma_deltaG(i, j);
        end
        fprintf('%d:%d / %d\n', i, j, size(D,1))
    end
end

% Print out device observation connectedness
for i = 1:size(deltaG,1);
    fprintf('%d: ', i);
    for m = find(deltaG(i,:))
        fprintf('%d ',m);
    end
    fprintf('\n');
end

% Objective function aiming to determine individual gain values from estimated relative gain
objective = @(gain) GainSolve(gain, deltaG, sigma_deltaG);

% Gain should be between +/- 20 dB
lower = zeros(size(D,1), 1) - 20;
upper = zeros(size(D,1), 1) + 20;

% Perform simulated annealing, starting from average values
G0 = zeros(size(D,1), 1);
options = optimset('display','off');
G = simulannealbnd(objective, G0, lower, upper, options);

end
