function [ e ] = LocSolve( J, D, C )
%LocSolve Objective function for calculating location of an observation using distance
%   J [in] - Access point location (longitude, latitude)
%   D [in] - Vector of calculated distances to APs
%   C [in] - Matrix of AP locations (c_longitude, c_lattitude)
%   e [out] - Mean Squared Error between calculated and localised distance
%
%   Must be able to see at least 3 APs to uniquely trilaterate

%d_ij = sqrt((x_jx - c_ix)^2 + (x_jy - c_iy)^2)

error = zeros(size(D,1),1);

for i = 1:size(D,1)
    error(i) = D(i) - sqrt((J(1) - C(i,1))^2 + (J(2) - C(i,2))^2);
end

e = mean(error.^2);

end

