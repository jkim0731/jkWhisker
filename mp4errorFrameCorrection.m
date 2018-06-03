% fnlist = {'203','211','221','230','250','278','294','302', '311','320','340','358','378', '387','395', '413','422'};
% fnlist = dir('*.mp4');
ffmpegPath = 'C:\Users\shires\Documents\GitHub\jkWhisker\tracking_windows\ffmpeg';

for fi = 1 : length(fnlist)
% for fi = 1 
    fn = [fnlist{fi}, '.mp4'];
    vr = VideoReader(fn);
    nof = fix(vr.FrameRate*vr.Duration);
    img = zeros(vr.Height, vr.Width, 1, nof, 'uint8');
    for i = 1 : nof
        temp = readFrame(vr);
        img(:,:,1,i) = temp(:,:,1);
    end    
    ts = squeeze(mean(mean(img)));
    idx = find(ts < mean([max(ts), min(ts)]));
    if ~isempty(idx)
        img(:,:,1,idx) = img(:,:,1,idx-1);
        
        fntemp = [fnlist{fi}, 'temp.mp4'];
        fntest = [fnlist{fi}, 'test.mp4'];
        FileRename(fn,fntemp);
        vw = VideoWriter(fntest, 'MPEG-4');    
        vw.Quality = 100;
        open(vw)
        writeVideo(vw,img)
        close(vw)    
        mp4CorrectionCMD = sprintf('%s -i %s -b:v 5M -vcodec libx264 -c:a copy %s', ffmpegPath, fntest, fn);

        system(mp4CorrectionCMD);
    end
end

%%
fntest = '24.mp4';
fntest2 = '24test.mp4';
mp4CorrectionCMD = sprintf('%s -i %s -b:v 5M -vcodec libx264 -c:a copy %s', ffmpegPath, fntest, fntest2);    
system(mp4CorrectionCMD);
