function pole_available_frames(mouseName,sessionName,videoloc,varargin)


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
firsts = zeros(length(vlist),1);
lasts = zeros(length(vlist),1);
vnames = zeros(length(vlist),1);
for vind = 1 : length(vlist)       
    tic
    [~,fn,~] = fileparts(vlist(vind).name);
    if exist([fn,'_WT.mat'],'file') % Let's run this only on those with WT files (otherwise there has been an error in the video file, e.g., interruption)
        v = VideoReader(vlist(vind).name);

        nof = fix(v.FrameRate*v.Duration); % number of frames
        vv = zeros(v.height,v.width,nof,'uint8');
        pole_shadow = zeros(nof,1);
        for i = 1 : nof
            temp = readFrame(v);
            vv(:,:,i) = temp(:,:,1);    
            pole_shadow(i) = mean(mean(temp(1:round(v.height*height_factor),1:round(v.width*width_factor),1)));    
        end
        % figure, plot(1:length(pole_shadow),pole_shadow), hold on, plot(1:length(pole_shadow), ones(1,length(pole_shadow))*(mean(pole_shadow)), 'r-');

        pole_candidate = find(pole_shadow < mean(pole_shadow));
        bg_frames = find(pole_shadow > min([pole_shadow(1:100); pole_shadow(end-100:end)]));
        bg = imgaussfilt(mean(vv(:,:,bg_frames),3),2);

        % %%
        pole_im = zeros(round(v.height*height_factor),round(v.width*width_factor),length(pole_candidate));
        pole_bw = zeros(round(v.height*height_factor),round(v.width*width_factor),length(pole_candidate),'logical');
        edge_front = zeros(length(pole_candidate),1);
        for i = 1 : length(pole_candidate)
            pole_im(:,:,i) = bg(1:round(v.height*height_factor),1:round(v.width*width_factor)) - imgaussfilt(double(vv(1:round(v.height*height_factor),1:round(v.width*width_factor),pole_candidate(i))),2);
            pole_bw(:,:,i) = pole_im(:,:,i) > 3*std(std(pole_im(:,:,i)));
            a = regionprops(pole_bw(:,:,i),'Extrema','Area');
            maxArea = a(1).Area;
            maxInd = 1;
            if length(a) > 1
                for j = 2 : length(a)
                    if a(j).Area > maxArea
                        maxArea = a(j).Area;
                        maxInd = j;
                    end
                end
            end    
            edge_front(i) = a(maxInd).Extrema(5,2);
        end            

        [c, p] = histcounts(edge_front,min(edge_front):0.25:max(edge_front));
        [~,maxInd] = max(c);
        maxPos = p(maxInd-1);

        [~,fnum] = fileparts(vlist(vind).name);
        vnames(vind) = str2double(fnum);
        firsts(vind) = find(edge_front >= maxPos,1,'first') + pole_candidate(1) - 1;
        lasts(vind) = find(edge_front >= maxPos,1,'last') + pole_candidate(1) - 1;
        %%
    %     figure, plot(1:length(edge_front),edge_front), hold on,
    %     plot(first,edge_front(first),'r.','MarkerSize',20)
    %     plot(last,edge_front(last),'r.','MarkerSize',20)
        %%
    %     implay(pole_bw)
        fprintf('%s %s %d/%d done.\n', mouseName, sessionName, vind, length(vlist));            
    else
        firsts(vind) = NaN; lasts(vind) = NaN;
        fprintf('Whisker tracking error: %s %s %d/%d \n', mouseName, sessionName, vind, length(vlist));
    end
end
save('firsts_n_lasts.mat','firsts','lasts', 'vnames')
