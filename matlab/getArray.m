% Takes a 2D array imported from python output and turns it into a 3D
% distance matrix with a page corresponding to each frame.
%
% Parameters:
%   data: n_frames*n_res x n_res numeric matrix of interacting residue data
%   imported from Python
%   n_res: number of residues that are interating with the peptide.
% Returns:
%   arr: n_res x n_res x n_frames 3-D numeric matrix representing
%   inter-residue distances on each frame.
function arr=getArray(data, n_res)
    frames = size(data,1)/n_res;
    arr = zeros(n_res, n_res, frames);
    for frame = 1:frames
        sIdx = (frame-1)*n_res+1;
        arr(:,:,frame) = data(sIdx:(sIdx+n_res-1),:);
    end
end