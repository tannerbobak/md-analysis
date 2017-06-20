% Gets a matrix of the mean inter-residue distance between each pair of
% residues. The mean is longitudinal, i.e. it is the mean across frames,
% not within a single frame.
function B=getMeanMatrix(A)
    B = zeros(size(A,1),size(A,2));
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            B(i,j) = mean(A(i,j,:));
        end
    end
end