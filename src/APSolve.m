function [ C ] = APSolve( O )
%APSolve Calculate unique parameters of Access Point using observed RSSI
%   O [in]  - Matrix of RSSI values from known locations relevant to
%       to this Access Point (RSSI, X_long, X_lat)
%   C [out] - Vector of calculated AP parameters:
%       c longtitude, c latitude, transmit power, path loss rate
%
%   Requires at least 5 fixed locations (rows in O) to solve all 4 AP parameters
%   If not enough fixed locations are visible, use APRandomInit

assert(size(O,1) > 4, 'Must have at least 5 fixed locations for unique solution');

%p_ij = P_i - (10* gamma_i)*log(d_ij) + R

%d_ij = sqrt((x_jx - c_ix)^2 + (x_jy - c_iy)^2)

syms c_ix c_iy P_i gamma_i;
equations = [];

for j = 1:size(O,1)
    eqn = P_i - ((10 * gamma_i) * log(sqrt((O(j,2) - c_ix)^2 + (O(j,3) - c_iy)^2))) == O(j,1);
    equations = [equations, eqn];
end
C = solve(equations, [c_ix, c_iy, P_i, gamma_i]);
    
end

