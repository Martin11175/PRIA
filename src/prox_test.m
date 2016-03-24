% Proximity metric test script
% Author: mh881@york.ac.uk

% Load raw data from CSV
srcData = csvread('UJIndoorLoc/trainingData.csv', 1);

% Scale known location data to a local space
bounds = 50;
min_lat = min(srcData(:,521));
min_long = min(srcData(:,522));
srcData(:,521) = srcData(:,521) - (min_lat - bounds);
srcData(:,522) = srcData(:,522) - (min_long - bounds);

% Set data up for test
J = srcData(:, 1:520);
K = srcData(:, 528);
D = sort(unique(K));
X = srcData(:,521:523);

% Ground truth params
g_dist = 1;
% EZ Params
ez_prox_threshold = 3;
min_AP_overlap = 2;
% Heuristic Params
min_AP_strength = -100;
num_strongest_APs = 4;
min_strong_overlap = 3;

% Data prep
J(J == 100) = -100; % Replace positive invisibility markers
J(J < min_AP_strength) = -100;
max_J = zeros(size(J,1), num_strongest_APs);
for i = 1:size(J,1)
    [~, order] = sort(J(i,:), 'descend');
    max_J(i,:) = order(1:num_strongest_APs);
end

ez_result = zeros(size(D,1),3);
h_result = zeros(size(D,1),3);

% Test device by device
count_d = 1;
for i = nchoosek(1:size(D,1), 2)'
    x_1 = X(K == D(i(1)), :);
    x_2 = X(K == D(i(2)), :);
    max_1 = max_J(K == D(i(1)), :);
    max_2 = max_J(K == D(i(2)), :);
    k_1 = J(K == D(i(1)), :);
    k_2 = J(K == D(i(2)), :);
    
    % Ground truth proximate on the same floor within 1m
    tic
    g_prox = (abs(bsxfun(@minus, x_1(:,1), x_2(:,1)')) < g_dist ...
        & abs(bsxfun(@minus, x_1(:,2), x_2(:,2)')) < g_dist ...
        & bsxfun(@eq, x_1(:,3), x_2(:,3)'));
    [g_m, g_n] = find(g_prox);
    toc
    
    % EZ proximity by similarity
    tic
    ez_avg_diff = zeros(sum(K == D(i(1))), sum(K == D(i(2))));
    ez_count_m = 1;
    for ez_m = J(K == D(i(1)), :)'
        ez_count_n = 1;
        for ez_n = J(K == D(i(2)), :)'
            % Only compare APs visible at at least one location
            vis = (ez_m > -100) | (ez_n > -100);
            if(sum(vis) > min_AP_overlap && mean(abs(ez_m(vis) - ez_n(vis))) < ez_prox_threshold)
                ez_avg_diff(ez_count_m, ez_count_n) = mean(ez_m(vis) - ez_n(vis));
            end
            ez_count_n = ez_count_n + 1;
        end
        ez_count_m = ez_count_m + 1;
    end
    ez_prox = abs(ez_avg_diff) < ez_prox_threshold & (ez_avg_diff ~= 0);
    [ez_m, ez_n] = find(ez_prox);
    toc
    
    % Heuristic proximity by strongest APs
    tic
    h_prox = zeros(size(max_1, 1), size(max_2, 1));
    for h_m = 1:size(max_1,1)
        if sum(k_1(h_m, max_1(h_m,:))) > (num_strongest_APs * min_AP_strength)
        for h_n = 1:size(max_2,1)
            if sum(ismember(max_1(h_m,:), max_2(h_n,:))) > min_strong_overlap
                if sum(ismember([h_m h_n], [g_m g_n], 'rows')) == 0
                    dist_diff = sqrt((x_1(h_m,1) - x_2(h_n,1))^2 + (x_1(h_m,2) - x_2(h_n,2))^2);
                    fprintf('False positive distance: %f, %d floors', dist_diff, x_1(h_m,3) - x_2(h_n,3));
                end
                h_prox(h_m, h_n) = 1;
            end
        end
        end
    end
    [h_m, h_n] = find(h_prox);
    toc
    
    % TODO: For each point in h_prox, limit to those with avg_diff of 3dB
    
    ez_correct = sum(ismember([g_m, g_n], [ez_m, ez_n], 'rows'));
    h_correct = sum(ismember([g_m, g_n], [h_m, h_n], 'rows'));
    ez_missed = size(g_m,1) - ez_correct;
    h_missed = size(g_m,1) - h_correct;
    ez_false_positive = size(ez_m,1) - ez_correct;
    h_false_positive = size(h_m,1) - h_correct;
    
    ez_result(count_d,:) = [ez_correct, ez_missed, ez_false_positive];
    h_result(count_d,:) = [h_correct, h_missed, h_false_positive];
    [ez_result(count_d,:) h_result(count_d,:)]
    count_d = count_d + 1;
end