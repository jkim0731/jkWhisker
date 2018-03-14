function pole_available_frames_corr(mouseName,sessionName,videoloc,varargin)

% implemented in WhiskerSignalTrial_2pad 2018/03/06 JK

%% Manually set factors. Should be changed individually for different settings
width_factor = 1;
height_factor = 0.7;

%% Setup whisker array builder 
% mouseName = 'JK025'
% sessionName ='S03'
% videoloc = 'Y:\Whiskernas\JK_temp\whisker\tracked\';
% optional = 'noskip'
optional = 'skip';
if nargin > 3
    optional = varargin{1};
elseif nargin > 4
    error('Too many input arguments')
end

d = ([videoloc filesep mouseName sessionName filesep])
cd(d)

paf_fn = [mouseName sessionName 'pole_available_frames.mat'];
if exist(paf_fn,'file') && strcmp(optional, 'skip')
    disp([paf_fn ' already exists. Skipping.'])    
    return
elseif exist(paf_fn,'file') && strcmp(optional, 'noskip')
    disp([paf_fn ' already exists. Overriding.'])    
end
%%
vlist = dir('*.mp4');
pole_up_frames = struct;
parfor vind = 1:length(vlist)
    tic
    [~,fn,~] = fileparts(vlist(vind).name);
    if exist([fn,'_WT.mat'],'file') % Let's run this only on those with WT files (otherwise there has been an error in the video file, e.g., interruption)
        %%
        v = VideoReader(vlist(vind).name);

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
        
        pole_up_frames(vind).name = str2double(fn);
        pole_up_frames(vind).frames = find(im_corr >= threshold,1,'first') : find(im_corr >= threshold,1,'last');
        fprintf('%s %s %d/%d done.\n', mouseName, sessionName, vind, length(vlist));            
    else
        pole_up_frames(vind).frames = NaN;
        fprintf('Whisker tracking error: %s %s %d/%d \n', mouseName, sessionName, vind, length(vlist));
    end
end
savefn = sprintf('pole_up_frames_%s%s.mat', mouseName, sessionName);
save(savefn,'pole_up_frames')
