function Uy = normalize(Ux)
% normalize each col of Ux
% Args:
%   Ux  : input matrix
% Returns:
%   Uy  :  matrix with normalized cols
Uy = Ux * diag(1./sqrt(sum(Ux.^2,1)));
end