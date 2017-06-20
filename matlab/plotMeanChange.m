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
    ylabel({'Difference between overall mean inter-residue distance and'; ...
            'current mean inter-residue distance (nm)'});
end