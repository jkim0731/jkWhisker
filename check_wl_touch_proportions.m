wla = Whisker.WhiskerTrialLite_2padArray([whiskerDir,mouseName,sessionName]);
validInd = find(cellfun(@(x) ~strcmp(x.trialType, 'oo'),wla.trials));
ntInd = find(cellfun(@(x) isempty(union(x.protractionTouchFrames, x.retractionTouchFrames)), wla.trials(validInd)));
tns = cellfun(@(x) x.trialNum, wla.trials(validInd(ntInd)));

%%

height = size(ws.binvavg,1);
width = size(ws.binvavg,2);
topViewPole = zeros([height, width],'logical');
targetH = 0.8; targetW = 0.8; targetArea = zeros(size(ws.binvavg), 'logical'); targetArea(1:round(height*targetH), round((1-targetW)*width):end) = deal(1);
tempBinvavg = ws.binvavg;
tempBinvavg(:,1:round(width*(1-targetW))) = deal(0);
tempBinvavg(round(height*targetH):end, :) = deal(0);
targetInd = find(targetArea);                    
bw = bwconncomp(tempBinvavg);
bwInd = find(cellfun(@(x) length(intersect(x,targetInd)),bw.PixelIdxList));
if isempty(bwInd)
    error(['no top-view pole detected at mouse ', obj.mouseName, ' session ', obj.sessionName, ' trial # ', num2str(obj.trialNum)]);
end
for tempi = 1 : length(bwInd)
    topViewPole(bw.PixelIdxList{bwInd(tempi)}) = deal(1);
end
                    
                    

obviousNoTouchFrames = zeros(wl.nof,1,'logical');
for tempi = 1 : length(wl.poleUpFrames)
    trackerFrameTop = find(ws.trackerFrames{1} == wl.poleUpFrames(tempi),1);
    if ~isempty(trackerFrameTop)
        x = polyval(ws.polyFits{1}{1}(trackerFrameTop,:),linspace(0,1.3)); % stretch the whisker fitting outwardly 30 % for cases where whisker tracing is cut off because of the pole
        y = polyval(ws.polyFits{1}{2}(trackerFrameTop,:),linspace(0,1.3));
        xyi = find(y >= 1 & y <= height & x >= 1 & x <= width); % take only those within image dimension
        x = x(xyi); y = y(xyi);
        whiskerBW = zeros([height, width],'logical');
        whiskerBW(sub2ind(size(whiskerBW),round(y),round(x))) = deal(1);
        whiskerPoleIntersection = intersect(find(whiskerBW), find(topViewPole));
        if isempty(whiskerPoleIntersection)
            obviousNoTouchFrames(wl.poleUpFrames(tempi)) = 1;
        end
    end
end

height = size(ws.binvavg,1);
                    width = size(ws.binvavg,2);
                    topViewPole = zeros([height, width],'logical');
                    targetH = 0.8; targetW = 0.8; targetArea = zeros(size(ws.binvavg), 'logical'); targetArea(1:round(height*targetH), round((1-targetW)*width):end) = deal(1);
                    tempBinvavg = ws.binvavg;
                    tempBinvavg(:,1:round(width*(1-targetW))) = deal(0);
                    tempBinvavg(round(height*targetH):end, :) = deal(0);
                    targetInd = find(targetArea);                    
                    bw = bwconncomp(tempBinvavg);
                    bwInd = find(cellfun(@(x) length(intersect(x,targetInd)),bw.PixelIdxList));
                    if isempty(bwInd)
                        error(['no top-view pole detected at mouse ', obj.mouseName, ' session ', obj.sessionName, ' trial # ', num2str(obj.trialNum)]);
                    end
                    for tempi = 1 : length(bwInd)
                        topViewPole(bw.PixelIdxList{bwInd(tempi)}) = deal(1);
                    end