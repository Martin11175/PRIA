function [ Z ] = HeirarchicalCluster( X, rank )
%HEIRARCHICALCLUSTER Clustering as specified by EZ for APs and locations
%   X [in]  - Normalised data matrix, clusters by columns
%   rank [in] - Ranking of all columns indicating preference when selecting 
%       for cluster representative (higher rank = higher preference)
%   Z [out] - Array of columns elected to represent clusters

Z = 1:size(X,2); % Mapping of all columns to clusters
max_similarity = 1;
threshold_similarity = 0.9; % Point at which to cluster

% Iterate until no more similar clusters
while max_similarity > threshold_similarity
    max_similarity = threshold_similarity; % Reset maximum observed similarity between clusters
    C = sort(unique(Z)); % Clusters to be evaluated
    
    % Calculate similarity and mark maximum on each pass
    for i = nchoosek(unique(C), 2)'
            % Similarity between clusters is average of all inter-cluster pairs
            % TODO: Difference matrix between vectors = bsxfun(@minus, x, x')
            pair_sim = zeros(sum(Z == i(1)), sum(Z == i(2)));
            count_x = 1;
            for x = X(:, Z == i(1))
                count_y = 1;
                for y = X(:, Z == i(2))
                    pair_sim(count_x, count_y) = 1 - (1 / size(X,2) * sum(abs(x - y)));
                    count_y = count_y + 1;
                end
                count_x = count_x + 1;
            end
            similarity = mean2(pair_sim);
            
            if similarity > max_similarity
                max_similarity = similarity;
                cluster = i; % Pair of clusters to merge
            end
            
        % Nothing can be more similar than exact (speed boost)
        if max_similarity == 1
            break;
        end
    end
    
    % Merge most similar cluster
    if max_similarity > threshold_similarity
        Z(Z == cluster(2)) = cluster(1);
        %fprintf('Merging %d / %d. %d clusters remaining. Similarity: %f\n', cluster(1), cluster(2), size(unique(Z), 2), max_similarity);
    end
end

% Select representative column for each cluster
for c = sort(unique(Z))
    % First select based on the max rank for the cluster
    rep = find(rank == max(rank(Z == c)) & (Z == c)');
    
    % Check for similarity to cluster if multiple with equivalent rank
    if size(rep,1) > 1
        best_avg_similarity = 0;
        
        % Calculate average similarity to the rest of the cluster
        for i = rep'
            avg_similarity = 0;
            for j = find(Z == c)
                if j ~= i
                    avg_similarity = avg_similarity + (1 - (1 / size(X,2) * sum(abs(X(:,i) - X(:,j)))));
                end
            end
            avg_similarity = avg_similarity / (sum(Z == c) - 1);
            
            if avg_similarity > best_avg_similarity
                best_avg_similarity = avg_similarity;
                rep = i;
            end
        end
    end
    
    % Replace cluster representative
    Z(Z == c) = rep;
end

Z = sort(unique(Z));

end

