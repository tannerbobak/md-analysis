% Makes an AVI file out of the plot shown in plotMe() and saves to disk.
% This way it can be viewed fast!
%
% Parameters:
%   arr: 3-D numeric array of inter-residue distances for each frame.
%   labels: Residue numbers that are interacting.
%   nstxout_compressed: Number of steps between writing to .xtc, taken from
%   the .mdp file.
%   dt: Time step (ps) for the simulation, taken from the .mdp file.
%   stride: The number of frames skipped between samples when analyzing the
%   trajectory in Python.
%   fps: Framerate of the animation.
%   startTime: The time (in ns) to start the animation from.
% Returns: None
function movieMaker(arr, labels, nstxout_compressed, dt, stride, fps, startTime)
    [file, path] = uiputfile('*.avi');
    v = VideoWriter(strcat(path,file));
    v.FrameRate = fps;
    open(v);
    
    close all
    grid on
    figure('Position',[10,10,1000,1000]);
    arr = arr - getMeanMatrix(arr);
    picosecondsPerFrame = dt*nstxout_compressed*stride;
    count = 1;
    if(startTime ~= 0)
        count = startTime*1000/picosecondsPerFrame;
    end
    
    for f = count:size(arr,3)
        scale = 100;
        colormap jet
        image(imresize(arr(:,:,f), scale, 'nearest'), 'CDataMapping','scaled');
        xticks(0:scale:scale*size(arr,2))
        yticks(0:scale:scale*size(arr,1))
        xlim([0,scale*size(arr,2)])
        ylim([0,scale*size(arr,1)])
        xticklabels(string(labels));
        yticklabels(string(labels));
        pbaspect([1,1,1])
        t = sprintf('%.2f',f*picosecondsPerFrame/1000);
        title(['Frame #', num2str(f), ' (', t, ' ns)']);
        ylabel('Residue Number');
        xlabel('Residue Number');
        caxis([minN(arr,3),maxN(arr,3)]);
        cb = colorbar;
        ylabel(cb, 'RMSD distance from mean (nm)');
        
        frame = getframe(gcf);
        writeVideo(v, frame);
        disp(strcat(['Wrote frame #', num2str(f), ' of ', num2str(size(arr,3)-count+1), '...']));
    end
    
    grid off
    close(v);
    
    disp('Done!');
end