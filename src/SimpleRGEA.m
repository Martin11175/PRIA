function [ G ] = SimpleRGEA( J, K )
%SimpleRGEA Simple Relative Gain Estimation Algorithm using strongest RSSIs
%   J [in] - Matrix of RSSI measurements (APs as columns)
%   K [in] - Vector of device IDs relating rows to device
%   G [out] - Pairs of device IDs (1) to estimated gain (2)

prox_threshold = 3;
min_AP_overlap = 2;
num_strongest_APs = 4;
min_AP_strength = -70;

D = sort(unique(K)); % List of device IDs
J(J == 100) = -100; % Replace positive invisibility markers

% Isolate strongest RSSI measurements per location
for i = 1:size(J,1)
    [~, order] = sort(J(i,:));
    J(i, order(1:(size(J,2) - num_strongest_APs))) = 0;
end
J(J < min_AP_strength) = 0;
J = sparse(J); % TODO: Is sparse faster?

% Calculate relative gain between pairs of devices
deltaG = zeros(size(D,1));
sigma_deltaG = zeros(size(D,1)); % Uncertainty (estimated standard deviation)

for i = nchoosek(1:size(D,1), 2)'
    % Test pairs of measurements for proximity by similarity
    diff = zeros(sum(K == D(i(1))), sum(K == D(i(2))), size(J,2));
    m = J(K == D(i(1)), :);
    n = J(K == D(i(2)), :);
    
    for k = 1:520
        diff(:,:,k) = bsxfun(@minus, m(:,k), n(:,k)');
    end
    
    tic
    % Compare all proximate locations
    num_compared = sum(diff ~= 0, 3);
    sum_abs_diff = sum(abs(diff),3);
    avg_abs_diff = (sum_abs_diff ./ num_compared);
    prox = (avg_abs_diff < prox_threshold) & (num_compared > min_AP_overlap);
    
    if sum(sum(prox)) > 10
        avg_diff = (sum(diff,3) ./ num_compared);
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
