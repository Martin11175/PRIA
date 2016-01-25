function [ J ] = APSolve( C, O )
%APSolve Objective function for calculating parameters of Access Point using observed RSSI
%   O [in]  - Matrix of RSSI values from known locations relevant to
%       to this Access Point (RSSI, X_long, X_lat)
%   C [in] - Vector of AP parameters:
%       c longtitude, c latitude, transmit power, path loss rate
%   J [out] - Error between observed and estimated RSSI
%
%   Requires at least 5 fixed locations (rows in O) to solve all 4 AP parameters
%   If not enough fixed locations are visible, use APRandomInit

assert(size(O,1) > 4, 'Must have at least 5 fixed locations for unique solution');

%p_ij = P_i - (10* gamma_i)*log(d_ij) + R

%d_ij = sqrt((x_jx - c_ix)^2 + (x_jy - c_iy)^2)

error = zeros(size(0,1));

for j = 1:size(O,1)
    error(j) = O(j,1) - C(3) + ((10 * C(4)) * log(sqrt((O(j,2) - C(1))^2 + (O(j,3) - C(2))^2)));
end

J = median(sort(abs(error)));

end

