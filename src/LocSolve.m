function [ e ] = LocSolve( J, O, C )
%LocSolve Objective function for calculating location of an observation using RSSIs
%   J [in] - Access point location (longitude, latitude)
%   O [in] - Vector of observed RSSI values
%   C [in] - Matrix of AP parameters:
%       c longtitude, c latitude, transmit power, path loss rate
%   e [out] - Error between observed and estimated RSSI
%
%   Must be able to see at least 3 APs to uniquely solve.

%p_ij = P_i - (10* gamma_i)*log(d_ij) + R

%d_ij = sqrt((x_jx - c_ix)^2 + (x_jy - c_iy)^2)

error = zeros(size(O,1),1);

for i = 1:size(O,1)
    error(i) = O(i) - C(i,3) + ((10 * C(i,4)) * log(sqrt((J(1) - C(i,1))^2 + (J(2) - C(i,2))^2)));
end

e = median(sort(abs(error)));

end

