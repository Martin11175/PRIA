function [ Z ] = HeirarchicalCluster( X )
%HEIRARCHICALCLUSTER Clustering as specified by EZ for APs and locations
%   X [in]  - Normalised data matrix, clusters by columns
%   Z [out] - Array of columns elected to represent clusters

Z = 1:size(X,1); % Mapping of all columns to clusters
max_similarity = 1;

% Iterate until no more similar clusters
while max_similarity > 0.9
    max_similarity = 0.9; % Maximum observed similarity between clusters
    C = sort(unique(Z)); % Clusters to be evaluated
    
    % Calculate similarity and mark maximum on each pass
    for i = 1:(size(C,1) - 1)
        for j = i:size(C,1) % Don't re-evaluate pairs
            similarity = 1 - (1/size(X,2) * abs(symsum(X(i,k) - X(j,k), k, 1, size(X,2))));
            if similarity > max_similarity
                max_similarity = similarity;
                cluster_from = [i j]; % Pair of clusters to merge
            end
        end
    end
    
    % Only compute once the most similar pair has been determined
    if max_similarity > 0.9
        % Decide which column should represent the cluster
        % (First by known location visibility, then by similarity)
        curr_cluster = X(Z == cluster_from(1) | Z == cluster_from(2));
        rep = 0; % New representing column
        
        % TODO: Implement known location check (remember to avoid including location data in similarity checks!)
        
        if rep == 0
            best_avg_similarity = 0;
        
            % Calculate average similarity to the rest of the cluster
            for i = curr_cluster
                avg_similarity = 0;
                for j = curr_cluster(curr_cluster ~= i)
                    avg_similarity = avg_similarity + (1 - (1/size(X,2) * abs(symsum(X(i,k) - X(j,k), k, 1, size(X,2)))));
                end
                avg_similarity = avg_similarity / (size(curr_cluster,1) - 1);
            
                if avg_similarity > best_avg_similarity
                    rep = i;
                end
            end
        end
        
        % Replace cluster representative
        Z(curr_cluster) = rep;
    end
end

Z = sort(unique(Z));

end

