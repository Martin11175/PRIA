function [ G ] = RGEA( J, K )
%RGEA Relative Gain Estimation Algorithm
%   J [in] - Matrix of RSSI measurements (APs as columns)
%   K [in] - Vector of device IDs relating rows to device
%   G [out] - Estimated gain by device ID

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
        
        count_m = 1;
        for m = J(K == D(i), :)'
            count_n = 1;
            for n = J(K == D(j), :)'
                % Only compare APs visible at at least one location
                vis = (m > -100) | (n > -100);
                avg_diff(count_m, count_n) = mean(m(vis) - n(vis));
                count_n = count_n + 1;
            end
            count_m = count_m + 1;
        end
        
        % Compare all proximate locations
        prox = abs(avg_diff) < 3;
        
        deltaG(i, j) = mean(avg_diff(prox));
        deltaG(j, i) = -deltaG(i, j);
        sigma_deltaG(i, j) = (1 / sum(sum(prox))) * sqrt(sum((avg_diff(prox) - deltaG(i,j)).^2));
        sigma_deltaG(j, i) = sigma_deltaG(i, j);
        fprintf('%d:%d / %d\n', i, j, size(D,1))
    end
end

% Calculate gain by device
G = zeros(size(D,1), 1);


end

