mp4list = dir('*.mp4');
ffmpegPath = 'C:\Users\shires\Documents\GitHub\jkWhisker\tracking_windows\ffmpeg';
for i = 1 
    fn = mp4list(i).name;
%     renamefn = ['temp', fn];
%     FileRename(fn, renamefn);
    outfntmp = ['temp10', fn];
    outfn = ['10', fn];
    v = VideoReader(fn);
    vout = VideoWriter(outfntmp,'MPEG-4');
    vout.Quality = 100;    
    open(vout);
    while hasFrame(v)
        temp = readFrame(v);
        if length(size(temp)) > 2
            temp = squeeze(temp(:,:,1));
        end        
%         temp = adapthisteq(temp);
        writeVideo(vout, temp(:,1:250));
    end
    close(vout)
    mp4CorrectionCMD = sprintf('%s -i %s -b:v 800k -codec:v mpeg4 -c:a copy %s', ffmpegPath, outfntmp, outfn);
    system(mp4CorrectionCMD);
end