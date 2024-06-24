function eff = nodal_efficiency(D)
%NODAL_EFFICIENCY
% Input: D, distance matrix (output from distance_wei)
% Output: eff, vector of nodal efficiencies
% Author: Haatef Pourmotabbed

DD = 1./D;
DD(1:length(DD)+1:end) = 0; % clear main diagonal
eff = sum(DD,2)./(length(DD)-1);

end