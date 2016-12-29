function avi2tiff(fn)
%
%
%
% DHO, 12/08.
%
% fn = 'JF4793_091407_S1_trial64.avi';
outfn = [fn(1:(end-3)) 'tif'];
m = aviread(fn);
ai = aviinfo(fn);

numFrames = ai.NumFrames;

s = ['NumFrames: ' int2str(ai.NumFrames) ', FPS: ' int2str(ai.FramesPerSecond) ...
    ', FileModDate: ' ai.FileModDate ...
    ', Width: ' int2str(ai.Width) ', Height: ' int2str(ai.Height)];

for k=1:numFrames
    imwrite(m(k).cdata,outfn,'TIFF','Compression','none',...
        'WriteMode','append','Description',s);
%     if mod(k,100)
%         disp(['Frame: ' int2str(k)]);
%     end
end

