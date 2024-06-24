function [syn,ac,S,s] = syn_wu(W)
%SYN_WU     Spectral measure of network synchronizability
%
%      Inputs:     W,     (symmetric) weighted undirected adjacency matrix
%
%      Outputs:    syn,   synchronizability (ratio of the second smallest eignvalue to the largest eigenvalue of the Laplacian matrix)
%                  ac,    algebraic connectivity (second smallest eigenvalue of the Laplacian matrix)
%                  s,     column vector of nodal strengths (for each node, the sum of all of its weights)
%                  S,     network strength (sum of all of the nodal strengths)
%
%      Author:     Haatef Pourmotabbed

n = length(W);
W(1:n+1:end) = 0; % set weights on the main diagonal (self-connections) to zero

s = sum(W,2);
S = sum(s);
L = diag(s) - W; % combinatorial Laplacian (Kirchhoff) matrix
[~,E,~] = eig(L);
e = sort(diag(E),'descend');
ac = e(end-1);
syn = ac/e(1);

end