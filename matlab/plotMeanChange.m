% Creates a plot showing how the difference between the mean inter-residue
% distance and the current inter-residue distance evolves with time.
%
% Parameters:
%   arr: 3-D numeric array of inter-residue distances for each frame.
%   nstxout_compressed: Number of steps between writing to .xtc, taken from
%   the .mdp file.
%   dt: Time step (ps) for the simulation, taken from the .mdp file.
%   stride: The number of frames skipped between samples when analyzing the
%   trajectory in Python.
% Returns: None
function plotMeanChange(arr,nstxout_compressed,dt,stride)
    close all
    figure(1);
    n = size(arr,3);
    picosecondsPerFrame = dt*nstxout_compressed*stride;
    means = zeros(n,1);
    totMean = mean2(getMeanMatrix(arr));
    for frame = 1:n
        means(frame) = totMean - mean2(arr(:,:,frame));
    end
    t = ((1:n)*picosecondsPerFrame)/1000;
    plot(t, means);
    title('Inter-residue distance difference vs Time');
    xlabel('Time Elapsed (ns)');
    % It is plotting 'totMean - mean2(arr(:,:,frame))' on the y-axis
    ylabel({'Difference between overall mean inter-residue distance and'; ...
            'current mean inter-residue distance (nm)'});
end