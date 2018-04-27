function [nof, poleUpFrames, poleMovingFrames, poleAxesUp, poleAxesMoving, topPix, barPos] = pole_edge_detection(videoFn, angle, radius)

% Automatic detection of pole_edge from video. Both front and top view.
% Currently only for 2pad 2017 JK

% Outputs:
% (1) total number of frames
% (2) pole up frames
% (3) pole moving frames
% (4) poleAxis when up
% (5) poleAxis durin moving
% (6) pixel value (bottom right) of the front pole where it has fully risen
% (to calculate fitting for ap position estimation later)

% Updates:
% 2018/03/06 Added automatic pole_available_frames 
% 2018/04/11 Using imbinarize, calculating bottom-left tip pixel point for
% pole_position estimation and calculating pole available timepoints
%% Fixed parameters
% targeting right top 1/4 of the whole FOV for top-view pole tracking
wFactorTop = 0.5;
hFactorTop = 0.5;
% targeting left top ~1/5.6 of the whole FOV for front-view pole tracking
wFactorFront = 0.3;
hFactorFront = 0.55;

%% Initialization
if isnumeric(videoFn)
    v = VideoReader([num2str(videoFn),'.mp4']);
else
    if length(videoFn) > 4 && strcmp(videoFn(end-3:end),'.mp4')
        v = VideoReader(videoFn);
    else
        v = VideoReader([videoFn,'.mp4']);
    end
end
targetWidth = round(v.Width*wFactorTop):v.Width-10;
targetHeight = 10:round(v.Height*hFactorTop);
rowSub = ones(length(targetWidth),1)*targetHeight;
rowSub = rowSub(:);
colSub = repmat(targetWidth,[1,length(targetHeight)]);
topTargetInd = sub2ind([v.Height, v.Width], rowSub, colSub');

targetWidth = 10:round(v.Width*wFactorFront);
targetHeight = 10:round(v.Height*hFactorFront);
rowSub = ones(length(targetWidth),1)*targetHeight;
rowSub = rowSub(:);
colSub = repmat(targetWidth,[1,length(targetHeight)]);
frontTargetInd = sub2ind([v.Height, v.Width], rowSub, colSub');

nof = fix(v.FrameRate*v.Duration);

topPix = NaN(nof,2); % bottom-right of top-view pole
topPixforcheck = NaN(nof,2); % right-bottom of top-view pole to check if the pole was in the FOV
frontPix = NaN(nof,2); % left-bottom of front-view pole

%% Gathering frame-by-frame information
for i = 1 : nof
    temp = readFrame(v);
    if  length(size(temp)) > 2 % temporary solution for having RGB-like mp4 file 2018/03/16 JK
        temp = temp(:,:,1);
    end
    btemp = 1 - imbinarize(uint8(temp), 'adaptive','ForegroundPolarity','dark','Sensitivity',0.1);
    btemp(:,1:20) = deal(0);
    btemp(:,end-10:end) = deal(0);
    bcc = bwconncomp(btemp);
    candid = find(cellfun(@(x) length(intersect(x,topTargetInd)), bcc.PixelIdxList));
    if ~isempty(candid) 
        btempTop = zeros(size(btemp),'logical');
        for j = 1 : length(candid)
            btempTop(bcc.PixelIdxList{candid(j)}) = 1;
        end
        bccTop = bwconncomp(btempTop);

        if bccTop.NumObjects > 1
            [~,bccind] = max(cellfun(@(x) length(x), bccTop.PixelIdxList));
        else
            bccind = 1;
        end
        btempTop = zeros(size(btempTop),'logical');
        btempTop(bccTop.PixelIdxList{bccind}) = 1;
        bccTop = bwconncomp(btempTop);
%         if length(bccTop.PixelIdxList{1}) > 10
            s = regionprops(bccTop,'Extrema');        
            topPix(i,:) = (floor(s.Extrema(5,:)) + floor(s.Extrema(6,:)))/2;
            topPixforcheck(i,:) = floor(s.Extrema(4,:));
%         end
    end    
    
    candid = find(cellfun(@(x) length(intersect(x,frontTargetInd)), bcc.PixelIdxList));
    if ~isempty(candid) 
        btempFront = zeros(size(btemp),'logical');
        for j = 1 : length(candid)
            btempFront(bcc.PixelIdxList{candid(j)}) = 1;
        end
        bccFront = bwconncomp(btempFront);

        if bccFront.NumObjects > 1
            [~,bccind] = max(cellfun(@(x) length(x), bccFront.PixelIdxList));
        else
            bccind = 1;
        end
        btempFront = zeros(size(btempFront),'logical');
        btempFront(bccFront.PixelIdxList{bccind}) = 1;
        bccFront = bwconncomp(btempFront);    
%         if length(bccFront.PixelIdxList{1}) > 50
            s = regionprops(bccFront,'Extrema');
            frontPix(i,:) = floor(s.Extrema(7,:));
%         end
    end        
end

poleUpPix = mode(topPix(:,1));
poleUpFrames = find(topPix(:,1) <= poleUpPix + 1, 1, 'first') : find(topPix(:,1) <= poleUpPix + 1, 1, 'last'); % for just in case where pixel values are noisy
poleMovingFrames = setdiff(     find(topPix(:,1)),    union(poleUpFrames,  union( find(isnan(topPix(:,1))), find(isnan(frontPix(:,1))) )  )     );
poleAxesMoving = cell(length(poleMovingFrames),2);
%% Binarized image from averaged pole up images (~ 10 frames, evenly distributed across pole up frames)
% And calculate pole edge slopes
frames = poleUpFrames(5): round((poleUpFrames(end-5) - poleUpFrames(5))/10) : poleUpFrames(end-5);
vavg = zeros(v.height,v.width);
for i = 1 : length(frames)
    v.CurrentTime = frames(i)/v.FrameRate;
    temp = readFrame(v);
    if  length(size(temp)) > 2 % temporary solution for having RGB-like mp4 file 2018/03/16 JK
        temp = temp(:,:,1);
    end
    vavg = vavg + double(temp)/length(frames);    
end
btemp = 1 - imbinarize(uint8(vavg), 'adaptive','ForegroundPolarity','dark','Sensitivity',0.1);
btemp(:,1:20) = deal(0);
btemp(:,end-10:end) = deal(0);
bcc = bwconncomp(btemp);

candid = find(cellfun(@(x) length(intersect(x,topTargetInd)), bcc.PixelIdxList));
if ~isempty(candid) 
    topPole = zeros(size(btemp),'logical');
    for j = 1 : length(candid)
        topPole(bcc.PixelIdxList{candid(j)}) = 1;
    end
    bccTop = bwconncomp(topPole);

    if bccTop.NumObjects > 1
        [~,bccind] = max(cellfun(@(x) length(x), bccTop.PixelIdxList));
    else
        bccind = 1;
    end
    topPole = zeros(size(topPole),'logical');
    topPole(bccTop.PixelIdxList{bccind}) = 1;
    bccTop = bwconncomp(topPole);
    s = regionprops(bccTop,'Extrema');
    if angle~=90 % if it's NOT 90 degrees
        topPole(:,floor(s.Extrema(4,1))-7:end) = 0; % So, remove right-most 8 columns of the image from the pole. To remove tip (or kink) noise.
        bccTop = bwconncomp(topPole);

        if bccTop.NumObjects > 1
            maxYval = zeros(bccTop.NumObjects,1);
            for i = 1 : bccTop.NumObjects
                [yval, ~] = ind2sub([v.height, v.width],bccTop.PixelIdxList{i});
                maxYval = max(yval);
            end
            [~, maxInd] = max(maxYval);
            bccind = maxInd;
        else
            bccind = 1;
        end
        topPole = zeros(size(topPole),'logical');
        topPole(bccTop.PixelIdxList{bccind}) = 1;
        bccTop = bwconncomp(topPole);
        s = regionprops(bccTop,'Extrema');        
        topSlope = (s.Extrema(4,2) - s.Extrema(5,2))/(s.Extrema(4,1) - s.Extrema(5,1));
        
    else % if it's 90 degrees, calculate slope based on the pole movement
        % 5 frames before and after pole up
        frames = [poleUpFrames(1)-5:poleUpFrames(3),poleUpFrames(end-2):poleUpFrames(end)+5];
        p = polyfit(topPix(frames,1),topPix(frames,2),1); % linear fitting. p(1) is going to be the slope
        topSlope = p(1);
    end
    
end

candid = find(cellfun(@(x) length(intersect(x,frontTargetInd)), bcc.PixelIdxList));
if ~isempty(candid) 
    frontPole = zeros(size(btemp),'logical');
    for j = 1 : length(candid)
        frontPole(bcc.PixelIdxList{candid(j)}) = 1;
    end
    bccFront = bwconncomp(frontPole);

    if bccFront.NumObjects > 1
        [~,bccind] = max(cellfun(@(x) length(x), bccFront.PixelIdxList));
    else
        bccind = 1;
    end
    frontPole = zeros(size(frontPole),'logical');
    frontPole(bccFront.PixelIdxList{bccind}) = 1;
    bccFront = bwconncomp(frontPole);
    s = regionprops(bccFront,'Extrema');
    frontPole(:,floor(s.Extrema(4,1))-2:end) = 0; % Remove right-most 3 columns of the image from the pole. 
    bccFront = bwconncomp(frontPole);
    s = regionprops(bccFront,'Extrema');
    frontSlope = (s.Extrema(4,2) - s.Extrema(7,2))/(s.Extrema(4,1) - s.Extrema(7,1));    
end

%% Calculte pole Up axes and adjust topPix, frontPix, and moving frames
% Based on pixel values (left-bottom for front-view and bottom-right for top-view) and the slopes calculated above

% for front
q = linspace(1,v.width*0.4);
middleUpFrame = poleUpFrames(round(length(poleUpFrames)/2));
poleAxesUp{2} = [frontPix(middleUpFrame,2) + q * frontSlope; q];
ind = find(frontPix(:,1)>0);
if ~isempty(ind)
    for i = 1 : length(ind)
        frontPix(ind(i),:) = [NaN, NaN];
    end
end

% for top
q = linspace(v.width,v.width*0.4);
originY = topPix(middleUpFrame,2) + (v.width - topPix(middleUpFrame,1)) * topSlope; % originX is at v.width
% if angle ~= 90
%     ind = find(topPixforcheck(:,2) < originY);
% else
    ind = find(topPix(:,2) < floor(originY));
% end

%%
% figure, imshow(topPole), hold on, plot(v.width,originY,'r.', 'MarkerSize', 20)


%%

if ~isempty(ind)
    for i = 1 : length(ind)
        topPix(ind(i), :) = [NaN, NaN];
    end
end

poleAxesUp{1} = [originY + (q-q(1)) * topSlope; q];

% adjust pole moving frames
poleMovingFrames = setdiff(   poleMovingFrames,    union( find(isnan(topPix(:,1))), find(isnan(frontPix(:,1))) )   );

%% calculate bar position when the angle is 90 degrees
if angle == 90
    polePresentFrames = union(poleMovingFrames, poleUpFrames); % union is sorted in default
    barPos = zeros(length(polePresentFrames),3);
    for i = 1 : size(barPos,1)
        barPos(i,1) = polePresentFrames(i);
        barPos(i,2) = topPix(polePresentFrames(i),1);
        barPos(i,3) = topPix(polePresentFrames(i),2) - radius;
    end
else
    barPos = [];
end

%% Calculte axes during pole movement
% for front
q = linspace(1,v.width*0.4);
for i = 1 : length(poleMovingFrames)
    poleAxesMoving{i,2} = [frontPix(poleMovingFrames(i),2)+ q * frontSlope; q];
end

% for top
q = linspace(v.width,v.width*0.4);
for i = 1 : length(poleMovingFrames)
    originY = topPix(poleMovingFrames(i),2) + (v.width - topPix(poleMovingFrames(i),1)) * topSlope;
    poleAxesMoving{i,1} = [originY + (q-q(1)) * topSlope; q];
end

end