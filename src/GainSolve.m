function [ e ] = GainSolve( G, deltaG, sigma_deltaG )
%GAINSOLVE Objective function for estimating gain from relative gain
%   G [in] - Vector of gain values by device
%   deltaG [in] - Matrix of relative gains between devices
%   sigma_deltaG [in] - Matrix of estimated standard deviation for relative
%       gains (weight for MSE)
%   e [out] - MSE between individual and relative gain estimates

error_matrix = zeros(size(deltaG,1), size(deltaG,2));

for i = 1:size(deltaG,1)
    for j = 1:size(deltaG,2)
        if deltaG(i,j) ~= 0
            % TODO: 1/sigma instead? Currently weighs greater deviation as higher
            % TODO: Tune the 100? Setting larger values will cause the
            %   simulated annealing to search longer?
            error_matrix(i,j) = 100 * sigma_deltaG(i,j) * (G(i) - G(j) - deltaG(i,j))^2;
        end
    end
end

e = mean2(error_matrix);

% TODO: Tune Occam's razor?
e = e + (0.1 * mean(abs(G)));

end

