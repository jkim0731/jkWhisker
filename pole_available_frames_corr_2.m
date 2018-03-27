function pole_up_frames = pole_available_frames_corr_2(fn)

% implemented in WhiskerSignalTrial_2pad 2018/03/06 JK
% modified for bartracker 2018/03/26 JK

%% Manually set factors. Should be changed individually for different settings
width_factor = 1;
height_factor = 0.7;
%%

v = VideoReader(fn);

nof = fix(v.FrameRate*v.Duration); % number of frames
h = floor(v.height*height_factor);
w =  floor(v.width*width_factor);
pole_img_ts = zeros(h,w,nof);
for i = 1 : nof
    temp = readFrame(v);
    pole_img_ts(:,:,i) = adapthisteq(temp(1:h,1:w,1),'NumTiles',[5,20]);    
end
ref_img = mean(pole_img_ts(:,:,floor(nof/2 - nof/10):floor(nof/2+nof/10)),3);
[pole_edge_ref, edge_thresh] = edge(imgaussfilt(ref_img,3),'Prewitt','nothinning');

im_corr = zeros(nof,1);
for i = 1 : nof
    pole_edge_ts = edge(imgaussfilt(pole_img_ts(:,:,i),3),'Prewitt', 'nothinning', edge_thresh);
    im_corr(i) = sum(sum(xcorr2(single(pole_edge_ts), single(pole_edge_ref))));
end                
corr_upper = im_corr(im_corr > (max(im_corr)+min(im_corr))/2);

threshold = prctile(corr_upper(floor(length(corr_upper)/5) : length(corr_upper) - floor(length(corr_upper)/5)),20);
pole_up_frames = find(im_corr >= threshold,1,'first') : find(im_corr >= threshold,1,'last');

