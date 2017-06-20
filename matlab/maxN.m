% Absolute maximum of a n-dimensional numeric matrix.
%
% Parameters:
%   A: A numeric matrix.
%   n: Number of dimensions that 'A' has.
% Returns:
%   B: Absolute max of the matrix.
function B=maxN(A,n)
    B = A;
    for i = 1:n
        B = max(B);
    end
end