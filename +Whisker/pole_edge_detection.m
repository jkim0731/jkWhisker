function [nof, poleUpFrames, poleMovingFrames, poleAxesUp, poleAxesMoving, topPix, barPos, binvavg] = pole_edge_detection(videoFn, angle, radius)

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
hFactorTop = 0.6;
% targeting left top ~1/5.6 of the whole FOV for front-view pole tracking
wFactorFront = 0.3;
hFactorFront = 0.7;

excludeHFactor = 0.8; % anything has pixel value under 0.8 height should be ignored (because it is the face)

% hard-coded padding for some front pole image error
topKinkPad = 7;
topTipPad = 2;
frontTipPad = 2;
frontLinkPad = 20;
% topLinkPad = 10;
% topExPad = 40; % sometimes top pole is divided into two, and linker gets chosen because of bulky body. To solve this, pad 0 at the top (only for top pole part)
% Instead, choose the lower one when there are multiple objects on the top view
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

excludeHeight = round(v.Height * excludeHFactor):v.Height;
excludeWidth = 1:v.Width;
rowSub = ones(length(excludeWidth),1)*excludeHeight;
rowSub = rowSub(:);
colSub = repmat(excludeWidth,[1,length(excludeHeight)]);
excludeHInd = sub2ind([v.Height, v.Width], rowSub, colSub');

nof = fix(v.FrameRate*v.Duration);

topPix = NaN(nof,2); % bottom-right of top-view pole
topPixforcheck = NaN(nof,2); % right-bottom of top-view pole to check if the pole was in the FOV
frontPix = NaN(nof,2); % left-bottom of front-view pole

%% Gathering frame-by-frame information
for i = 1 : nof
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % For debugging
%     %%
%     i=400;
%     v.CurrentTime = i/v.FrameRate;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp = readFrame(v);
    if  length(size(temp)) > 2 % temporary solution for having RGB-like mp4 file 2018/03/16 JK
        temp = temp(:,:,1);
    end
    btemp = 1 - imbinarize(uint8(temp), 'adaptive','ForegroundPolarity','dark','Sensitivity',0.1);
    btemp(:,1:frontLinkPad) = deal(0);
%     btemp(:,end-topLinkPad:end) = deal(0);
    bcc = bwconncomp(btemp);
    candid = find(cellfun(@(x) length(intersect(x,topTargetInd)), bcc.PixelIdxList));

%     % for debuggin
%     [topi,topj] = ind2sub(size(btemp),topTargetInd);
%     topTargetArea = zeros(size(btemp),'logical');
%     for iii = 1:length(topi)
%         topTargetArea(topi(iii),topj(iii)) = 1;
%     end
%     figure, imshowpair(topTargetArea, btemp)
    
    if ~isempty(candid) 
        btempTop = zeros(size(btemp),'logical');
        for j = 1 : length(candid)
            if isempty(intersect(bcc.PixelIdxList{candid(j)},excludeHInd))
                btempTop(bcc.PixelIdxList{candid(j)}) = 1;
            end
        end
%         btempTop(1:topExPad,:) = deal(0);
        
        bccTop = bwconncomp(btempTop);
        if bccTop.NumObjects
            if bccTop.NumObjects > 1
                lowestPix = zeros(bccTop.NumObjects,1);
                for lpi = 1 : bccTop.NumObjects
                    [yval,~] = ind2sub(size(btempTop),bccTop.PixelIdxList{lpi});
                    lowestPix(lpi) = max(yval);
                end
%                 [~,bccind] = max(cellfun(@(x) length(x), bccTop.PixelIdxList));
                [~,bccind] = max(lowestPix);
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
    end    
    
    candid = find(cellfun(@(x) length(intersect(x,frontTargetInd)), bcc.PixelIdxList));
    if ~isempty(candid) 
        btempFront = zeros(size(btemp),'logical');
        for j = 1 : length(candid)
            if isempty(intersect(bcc.PixelIdxList{candid(j)},excludeHInd))
                btempFront(bcc.PixelIdxList{candid(j)}) = 1;
            end
        end
        bccFront = bwconncomp(btempFront);

        if bccFront.NumObjects
            if bccFront.NumObjects > 1
%                 rightMostPix = zeros(bccFront.NumObjects,1);
%                 for rmpi = 1 : bccFront.NumObjects
%                     [~,xval] = ind2sub(size(btempFront),bccTop.PixelIdxList{rmpi});
%                     rightMostPix(rmpi) = max(xval);
%                 end
                [~,bccind] = max(cellfun(@(x) length(x), bccFront.PixelIdxList));
%                 [~,bccind] = max(rightMostPix);
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
    
%     % for debugging
%     figure, imshow(btempTop)
    
end
%%
poleFrames = find(~isnan(topPix(:,1)));
if ~isempty(find(diff(poleFrames)-1,1))
    chunks = [1; find(diff(poleFrames)-1)+1; length(poleFrames)];
    [~,chunkInd] = max(diff(chunks));      
    for i = 1 : length(chunks)-1
        if i ~= chunkInd
            topPix(poleFrames(chunks(i):chunks(i+1)-1),:) = deal(NaN);
        end
    end
end
%%        
[n, edges] = histcounts(topPix(:,1), floor(min(topPix(:,1))):ceil(max(topPix(:,1))));
edgeInd = find(n > max(n)/2, 1, 'last');
poleUpPix = edges(edgeInd+1);
% poleUpPix = mode(topPix(:,1));

poleUpFrames = find(topPix(:,1) <= poleUpPix, 1, 'first') : find(topPix(:,1) <= poleUpPix, 1, 'last'); % for just in case where pixel values are noisy
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
binvavg = btemp;
btemp(:,1:frontLinkPad) = deal(0);
% btemp(:,end-topLinkPad:end) = deal(0);
bcc = bwconncomp(btemp);

candid = find(cellfun(@(x) length(intersect(x,topTargetInd)), bcc.PixelIdxList));
topPole = zeros(size(btemp),'logical');
for j = 1 : length(candid)
    if isempty(intersect(bcc.PixelIdxList{candid(j)},excludeHInd))
        topPole(bcc.PixelIdxList{candid(j)}) = 1;
    end
end
% topPole(1:topExPad,:) = deal(0);
bccTop = bwconncomp(topPole);
if bccTop.NumObjects 
    [~,bccind] = max(cellfun(@(x) length(x), bccTop.PixelIdxList));
else
    btempToptemp = btemp;
    btempToptemp(floor(v.Height*0.7):end,:) = deal(0);
    btempToptemp(:,1:100) = deal(0);
    bccToptemp = bwconncomp(btempToptemp);

    candid = find(cellfun(@(x) length(intersect(x,topTargetInd)), bccToptemp.PixelIdxList));
    topPole = zeros(size(btempToptemp),'logical');
    for j = 1 : length(candid)
        if isempty(intersect(bccToptemp.PixelIdxList{candid(j)},excludeHInd))
            topPole(bccToptemp.PixelIdxList{candid(j)}) = 1;
        end
    end
    bccTop = bwconncomp(topPole);
    if bccTop.NumObjects 
        [~,bccind] = max(cellfun(@(x) length(x), bccTop.PixelIdxList));
    else
        figure, imshow(btemp), axis image, axis off, title(['Video ', videoFn])
        error(['Error in ', videoFn])
    end
end
topPole = zeros(size(btemp),'logical');
topPole(bccTop.PixelIdxList{bccind}) = 1;
s = regionprops(topPole,'Extrema');
if angle < 90
    topPole(:,1:floor(s.Extrema(6,1))+topTipPad) = 0; % to remove dirty pole tip
    topPole(:,floor(s.Extrema(4,1))-topKinkPad:end) = 0; % to remove kink region of the pole
elseif angle > 90
    topPole(:,floor(s.Extrema(4,1))-topTipPad:end) = 0; % to remove dirty pole tip
    topPole(:,1:floor(s.Extrema(6,1))+topKinkPad) = 0; % to remove kink region of the pole
end
% the above will help better estimate the pole edge
bccTop = bwconncomp(topPole);

if bccTop.NumObjects % Sometimes there can be no pole because of paw movements and other reasons
    if bccTop.NumObjects > 1
        lowestPix = zeros(bccTop.NumObjects,1);
        for lpi = 1 : bccTop.NumObjects
            [yval,~] = ind2sub(size(btempTop),bccTop.PixelIdxList{lpi});
            lowestPix(lpi) = max(yval);
        end
%         [~,bccind] = max(cellfun(@(x) length(x), bccTop.PixelIdxList));
        [~,bccind] = max(lowestPix);
    else
        bccind = 1;
    end

    topPole = zeros(size(topPole),'logical');
    topPole(bccTop.PixelIdxList{bccind}) = 1;
    s = regionprops(topPole,'Extrema');
    if angle~=90 % if it's NOT 90 degrees
%         topSlope = (s.Extrema(4,2) - mean([s.Extrema(5,2), s.Extrema(6,2)]))/(s.Extrema(4,1) - mean([s.Extrema(5,1),s.Extrema(6,1)]));
        topSlope = (s.Extrema(4,2) - s.Extrema(5,2))/(s.Extrema(4,1) - s.Extrema(5,1));
    else % if it's 90 degrees, calculate slope based on the pole movement
        % 5 frames before and after pole up       
        try
            frames = [poleUpFrames(1)-5:poleUpFrames(3),poleUpFrames(end-2):poleUpFrames(end)+5];
            p = polyfit(topPix(frames,1),topPix(frames,2),1); % linear fitting. p(1) is going to be the slope
            topSlope = p(1);
        catch
            error(['topPix error in ', videoFn])
        end
    end   
    topExtrema = (floor(s.Extrema(5,:)) + floor(s.Extrema(6,:)))/2;    

    %%
    candid = find(cellfun(@(x) length(intersect(x,frontTargetInd)), bcc.PixelIdxList));
    frontPole = zeros(size(btemp),'logical');
    for j = 1 : length(candid)
        if isempty(intersect(bcc.PixelIdxList{candid(j)},excludeHInd))
            frontPole(bcc.PixelIdxList{candid(j)}) = 1;
        end
    end
    bccFront = bwconncomp(frontPole);
    if bccFront.NumObjects > 0
        [~,bccind] = max(cellfun(@(x) length(x), bccFront.PixelIdxList));
    else 
        btempFronttemp = btemp;
        btempFronttemp(floor(v.Height*0.7):end,:) = deal(0);
        btempFronttemp(:,100:end) = deal(0);
        bccFronttemp = bwconncomp(btempFronttemp);

        candid = find(cellfun(@(x) length(intersect(x,frontTargetInd)), bccFronttemp.PixelIdxList));
        frontPole = zeros(size(btemp),'logical');
        for j = 1 : length(candid)
            if isempty(intersect(bccFronttemp.PixelIdxList{candid(j)},excludeHInd))
                frontPole(bccFronttemp.PixelIdxList{candid(j)}) = 1;
            end
        end
        bccFront = bwconncomp(frontPole);
        if bccFront.NumObjects > 0
            [~,bccind] = max(cellfun(@(x) length(x), bccFront.PixelIdxList));
        else
            figure, imshow(btemp), axis image, axis off, title(videoFn)
            error(['No front pole from video ', videoFn])
        end
    end
    frontPole = zeros(size(btemp),'logical');
    frontPole(bccFront.PixelIdxList{bccind}) = 1;
    s = regionprops(frontPole,'Extrema'); 
    frontPole(:,floor(s.Extrema(4,1))-frontTipPad:end) = 0; % Remove right-most 3 columns of the image from the pole.   
    bccFront = bwconncomp(frontPole);

    if bccFront.NumObjects
        [~,bccind] = max(cellfun(@(x) length(x), bccFront.PixelIdxList));
   
        frontPole = zeros(size(frontPole),'logical');
        frontPole(bccFront.PixelIdxList{bccind}) = 1;
        s = regionprops(frontPole,'Extrema'); 
        bccFront = bwconncomp(frontPole);
        if bccFront.NumObjects             
            frontSlope = (s.Extrema(4,2) - s.Extrema(7,2))/(s.Extrema(4,1) - s.Extrema(7,1));    
            frontExtrema = floor(s.Extrema(7,:));

            %% Calculte pole Up axes and adjust topPix, frontPix, and moving frames
            % Based on pixel values (left-bottom for front-view and bottom-right for top-view) and the slopes calculated above

            % for front
            q = linspace(-frontLinkPad,v.width*0.4-frontLinkPad);
            poleAxesUp{2} = [q + frontLinkPad; frontExtrema(1,2) + q * frontSlope]; % order is changed from [y;x] to [x;y] 2018/06/13 JK
%             poleAxesUp{2} = [frontExtrema(1,2) + q * frontSlope; q + frontLinkPad;];
            % for top
            q = linspace(v.width,v.width*0.4);
            originY = topExtrema(1,2) + (v.width - topExtrema(1,1)) * topSlope; % originX is at v.width
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

            poleAxesUp{1} = [q; originY + (q-q(1)) * topSlope]; % order is changed from [y;x] to [x;y] 2018/06/13 JK
%             poleAxesUp{1} = [originY + (q-q(1)) * topSlope; q]; 
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
                poleAxesMoving{i,2} = [q; frontPix(poleMovingFrames(i),2)+ q * frontSlope]; % order is changed from [y;x] to [x;y] 2018/06/13 JK
            end

            % for top
            q = linspace(v.width,v.width*0.4);
            for i = 1 : length(poleMovingFrames)
                originY = topPix(poleMovingFrames(i),2) + (v.width - topPix(poleMovingFrames(i),1)) * topSlope;
                poleAxesMoving{i,1} = [q; originY + (q-q(1)) * topSlope]; % order is changed from [y;x] to [x;y] 2018/06/13 JK
            end
        else
            barPos = [];
            poleAxesUp{2} = [];
            poleAxesUp{1} = [];
            poleAxesMoving = [];
            poleUpFrames = [];
            poleMovingFrames = [];
            topPix = [];
        end
    else
        barPos = [];
        poleAxesUp{2} = [];
        poleAxesUp{1} = [];
        poleAxesMoving = [];
        poleUpFrames = [];
        poleMovingFrames = [];
        topPix = [];
    end
else
    barPos = [];
    poleAxesUp{1} = [];
    poleAxesUp{2} = [];
    poleAxesMoving = [];
    poleUpFrames = [];
    poleMovingFrames = [];
    topPix = [];
end

end