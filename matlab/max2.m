% Absolute maximum of a 2-D numeric matrix.
%
% Parameters:
%   A: Input 2-D matrix
% Returns:
%   B: Absolute maximum of matrix.
function B=max2(A)
    B = max(max(A));
end