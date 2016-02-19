function [ G ] = SimpleRGEA( J, K )
%SimpleRGEA Simple Relative Gain Estimation Algorithm using strongest RSSIs
%   J [in] - Matrix of RSSI measurements (APs as columns)
%   K [in] - Vector of device IDs relating rows to device
%   G [out] - Pairs of device IDs (1) to estimated gain (2)

D = sort(unique(K)); % List of device IDs
J(J == 100) = -100; % Replace positive invisibility markers to prevent skew

% TODO: Isolate strongest RSSI measurements per location
[~, order] = sort(J(1,:));

% Calculate relative gain between pairs of devices
deltaG = zeros(size(D,1));
sigma_deltaG = zeros(size(D,1)); % Uncertainty (estimated standard deviation)

for i = nchoosek(1:size(D,1), 2)'
    % Test pairs of measurements for proximity by similarity
    avg_diff = zeros(sum(K == D(i(1))), sum(K == D(i(2))));
    
    count_m = 1;
    for m = J(K == D(i(1)), :)'
        count_n = 1;
        for n = J(K == D(i(2)), :)'
            % Only compare APs visible at at least one location
            vis = (m > -100) | (n > -100);
            if(sum(vis) > min_AP_overlap && mean(abs(m(vis) - n(vis))) < prox_threshold)
                avg_diff(count_m, count_n) = mean(m(vis) - n(vis));
            end
            count_n = count_n + 1;
        end
        count_m = count_m + 1;
    end
    
    % Compare all proximate locations
    prox = abs(avg_diff) < prox_threshold & (avg_diff ~= 0);
    
    if sum(sum(prox)) > 10
        deltaG(i(1), i(2)) = mean(avg_diff(prox));
        sigma_deltaG(i(1), i(2)) = (1 / sum(sum(prox))) * sqrt(sum((avg_diff(prox) - deltaG(i(1), i(2))).^2));
    end
    fprintf('%f | %d:%d / %d\n', deltaG(i(1), i(2)), i(1), i(2), size(D,1))
end

% Solve least mean squares set of simultaneous equations
[i, j] = find(deltaG);
C = zeros(size(i,1), size(D,1));
d = zeros(size(i,1),1);

% Prepare system of simultaneous equations weighted by estimated standard deviation
% TODO: Determine sets of devices that can be estimated relative to each
% other and fix a starting point for each cluster.
for k = 1:size(i,1)
    % TODO: Scale weight, sigma^-1 increases exponentially as reaching peak
    % of bell curve (relatively little confidence improvement)!
    weight = 1 / sigma_deltaG(i(k), j(k));
    C(k, [i(k) j(k)]) = 1 * weight;
    d(k) = deltaG(i(k),j(k)) * weight;
end

G = [D (C \ d)];

end
