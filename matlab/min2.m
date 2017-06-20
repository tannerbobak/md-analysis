% Absolute minimum of a 2-D numeric matrix.
%
% Parameters:
%   A: Input 2-D matrix
% Returns:
%   B: Absolute minimum of matrix.
function B=min2(A)
    B = min(min(A));
end