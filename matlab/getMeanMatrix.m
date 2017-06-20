% Gets a matrix of the mean inter-residue distance between each pair of
% residues. The mean is longitudinal, i.e. it is the mean across frames,
% not within a single frame.
%
% Parameters:
%   A: n_res x n_res x n_frames 3-D numeric matrix representing 
%   inter-residue distances on each frame.
% Returns:
%   B: n_res x n_res 2-D numeric matrix representing the mean inter-residue
%   distance for each residue pair across all frames.
function B=getMeanMatrix(A)
    B = zeros(size(A,1),size(A,2));
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            B(i,j) = mean(A(i,j,:));
        end
    end
end