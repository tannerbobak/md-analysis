function plotMe(arr,labels,nstxout_compressed,dt,stride, fps, startTime)
    close all
    hold on
    figure(1)
    picosecondsPerFrame = dt*nstxout_compressed*stride;
    count = 1;
    if(startTime ~= 0)
        count = startTime*1000/picosecondsPerFrame;
    end
    ptime = 1/fps;
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
        
        scale = 100;
        colormap jet
        image(imresize(arr(:,:,count), scale, 'nearest'), 'CDataMapping','scaled');
        xticks(0:scale:scale*size(arr,2))
        yticks(0:scale:scale*size(arr,1))
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
%             break;
            count = 1;
        end
        if(fps ~= -1)
            pause(ptime)
        end
    end
    hold off
end