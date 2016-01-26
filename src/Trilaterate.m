function [ X_long, X_lat ] = Trilaterate( O, C )
%Trilaterate Trilaterate a location using AP parameters and observed RSSI
%   O [in]       - Vector of observed RSSI values from the location
%       being trilaterated
%   C [in]       - Matrix of known Access Point parameters (APs as rows)
%   X_long [out] - Calculated longtitude of location
%   X_lat [out]  - Calculated latitude of location

% d_ij = 10^((P_i - p_ij)/(10*gamma_i))


end

