% Plots a matrix of the difference between the inter-residue distance and
% the mean inter-residue distance for each residue, and animates across all
% frames.
%
% Parameters:
%   arr: 3-D numeric array of inter-residue distances for each frame.
%   labels: Residue numbers that are interacting.
%   nstxout_compressed: Number of steps between writing to .xtc, taken from
%   the .mdp file.
%   dt: Time step (ps) for the simulation, taken from the .mdp file.
%   stride: The number of frames skipped between samples when analyzing the
%   trajectory in Python.
%   fps: Framerate of the animation, at least in theory. In practice MATLAB
%   is too slow, or rather this algorithm is too slow.
%   startTime: The time (in ns) to start the animation from.
%   repeat: Whether to loop animation or not.
% Returns: None
function plotMe(arr, labels, nstxout_compressed, dt, stride, fps, ...
    startTime, repeat)
    close all
    figure(1)
    arr = arr - getMeanMatrix(arr);
    picosecondsPerFrame = dt*nstxout_compressed*stride;
    count = 1;
    if(startTime ~= 0)
        count = startTime*1000/picosecondsPerFrame;
    end
    ptime = 1/fps;
    
    % Do resizing beforehand
    disp('Resizing array, hold on to your RAM!!!!');
    n_frames = size(arr,3);
    scale = 100;
    frames = cell(n_frames,1);
    for i = 1:n_frames
        frames{i} = imresize(arr(:,:,i), scale, 'nearest');
    end
    
%     Normalization and stuff, don't use anymore.
%     arr = (arr - minN(arr,3))/(maxN(arr,3) - minN(arr,3));
    
%     for i = 1:size(arr,1)
%         for j = 1:size(arr,2)
%             if (i >= j) 
%                 arr(i,j,:) = minN(arr,3); 
%             end
%         end
%     end
%     rot90(arr,3);
%     clear i j
    
    while true
        
        colormap jet
        image(frames{count}, 'CDataMapping','scaled');
        xticks(0:scale:scale*size(arr,2) + scale/2)
        yticks(0:scale:scale*size(arr,1) + scale/2)
        xlim([0,scale*size(arr,2)])
        ylim([0,scale*size(arr,1)])
        xticklabels(string(labels));
        yticklabels(string(labels));
        pbaspect([1,1,1])
        t = sprintf('%.2f',count*picosecondsPerFrame/1000);
        title(['Frame #', num2str(count), ' (', t, ' ns)']);
        ylabel('Residue Number');
        xlabel('Residue Number');
        caxis([minN(arr,3),maxN(arr,3)]);
        cb = colorbar;
        ylabel(cb, 'RMSD distance from mean (nm)');
        
        count = count+1;
        if count > size(arr,3)
            if(repeat) break;
            else count = 1; end
        end
        if(fps ~= -1)
            pause(ptime)
        end
    end
end