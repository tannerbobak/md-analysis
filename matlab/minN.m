% Absolute minimum of a n-dimensional numeric matrix.
%
% Parameters:
%   A: A numeric matrix.
%   n: Number of dimensions that 'A' has.
% Returns:
%   B: Absolute min of the matrix.
function B=minN(A,n)
    B = A;
    for i = 1:n
        B = min(B);
    end
end