function nft = buildWTandWST(mouseName, sessionName, d, bSession, ppm)
    nft = []; % network fail time
    curr_d = pwd;
    whisker_d = [d, mouseName, sessionName, filesep];   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    cd(whisker_d)
    delete *_WT.mat
    delete *_WST.mat
    delete *_errorWST.mat
    delete *_WL_2pad.mat
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    load([whisker_d, mouseName, sessionName, '_post.mat']) % for maskx and masky. Ignoring includef from this one. It will be re-defined.
    
    if ~isempty(bSession) % try only the ones with behavior session

        filelist=dir([whisker_d '*.measurements']);
        dirTrialNums=zeros(1,size(filelist,1));

        % %%
        % Assign the trial numbers to existing .measurements files in the directory
        % NOTE : This assumes that the .measurements files have leading numbers
        % corresponding to trial number in string positions 1:end-13 of the file
        % name. These index numbers may need to be changed to match up to the
        % numerical code of the trial number.  (2016/09/05 JK)

        for i=1:length(filelist)
            dirTrialNums(i)=str2double(filelist(i).name(1:end-13)); % extract out the trial number from each measurements file present in directory
        end
        trialNums = sort(dirTrialNums);
        trialNums = trialNums(~isnan(trialNums));
        trialNums = intersect(trialNums,bSession.trialNums); % try only the ones with behavior trials

        includef=cell(size(trialNums,1),1);
        for i = 1: length(trialNums)
            includef{i} = num2str(trialNums(i));
        end
    end
%%
while true
    try
        if size(maskx{1},1) > size(maskx{1},2)
            Whisker.makeAllDirectory_WhiskerTrial_2pad(whisker_d,[0 1],'mask', {[maskx{1}';masky{1}'],[maskx{2}';masky{2}']},...
                'trial_nums',trialNums,'include_files',includef,...
                'barRadius',3,'faceSideInImage', 'bottom', 'framePeriodInSec',1/310,...
                'imagePixelDimsXY',[width height],'pxPerMm',ppm,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward', 'behavior', bSession);
        else
            Whisker.makeAllDirectory_WhiskerTrial_2pad(whisker_d,[0 1],'mask', {[maskx{1};masky{1}],[maskx{2};masky{2}]},...
                'trial_nums',trialNums,'include_files',includef,...
                'barRadius',3,'faceSideInImage', 'bottom', 'framePeriodInSec',1/310,...
                'imagePixelDimsXY',[width height],'pxPerMm',ppm,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward', 'behavior', bSession);
        end
        break
    catch
        disp('It''s in pause because of network failure')
        nft = [nft; clock]; 
    end
end

while true
    try
        Whisker.makeAllDirectory_WhiskerSignalTrial_2pad(whisker_d,'include_files',includef,'polyRoiInPix',[ppm 6*ppm]);
        break
    catch
        disp('It''s in pause because of network failure')
        nft = [nft; clock]; 
    end 
end

% assigning frame-by-frame AP positions
wstList = dir('*_WST.mat');
wstNums = zeros(length(wstList),1);
for i = 1 : length(wstNums)
    wstNums(i) = str2double(wstList(i).name(1:end-8));
end

angles = unique(cellfun(@(x) x.servoAngle, bSession.trials));
distances = unique(cellfun(@(x) x.motorDistance, bSession.trials));
for i = 1 : length(angles)
    for j = 1 : length(distances)
        trialInds = find(cellfun(@(x) x.servoAngle == angles(i) && x.motorDistance == distances(j), bSession.trials));
        trialNums = zeros(length(trialInds),1);
        for k = 1 : length(trialNums)
            trialNums(k) = bSession.trials{trialInds(k)}.trialNum;
        end    
        trialNums = intersect(trialNums,wstNums);
        if ~isempty(trialNums)
            trialFns = cell(length(trialNums),1);
            for k = 1 : length(trialNums)    
                trialFns{k} = num2str(trialNums(k));
            end

            wsArray = Whisker.WhiskerSignalTrialArray_2pad(whisker_d,'include_files',trialFns);
    
            apPositions = zeros(length(wsArray),1);
            polePixVals = zeros(length(wsArray),2);
            for k = 1 : length(wsArray)
                apPositions(k) = wsArray.trials{k}.apUpPosition;
                polePixVals(k,:) = wsArray.trials{k}.topPix(wsArray.trials{k}.poleUpFrames(round(length(wsArray.trials{k}.poleUpFrames)/2)),:); % value at the center of poleUpFrames
            end    
%             Calculate linear fit between euclidean distance between pixel values and differences in pole position values
            [~, baseInd] = max(apPositions);
            pixelDistances = sum((polePixVals - polePixVals(baseInd)).^2,2).^0.5;
            positionDiff = apPositions - apPositions(baseInd);
            
            poleImP = polyfit(polePixVals(:,1), polePixVals(:,2), 1);
            mirrorAngle = atand(-poleImP(1));
            
            p = polyfit(pixelDistances , positionDiff, 1); % linear fitting
            slope = p(1);        
            for k = 1 : length(wsArray)
                apPosition = NaN(wsArray.trials{k}.nof,1);
                apPosition(wsArray.trials{k}.poleUpFrames) = wsArray.trials{k}.apUpPosition;
                for m = 1 : length(wsArray.trials{k}.poleMovingFrames)
                    apPosition(wsArray.trials{k}.poleMovingFrames(m)) = wsArray.trials{k}.apUpPosition + slope * sqrt(sum((wsArray.trials{k}.topPix(wsArray.trials{k}.poleMovingFrames(m))...
                                                                                                                 - wsArray.trials{k}.topPix(wsArray.trials{k}.poleUpFrames(round(length(wsArray.trials{k}.poleUpFrames)/2)))).^2));
                end
                load([wsArray.trials{k}.trackerFileName, '_WST.mat']) % loading ws
                ws.apPosition = apPosition;
                ws.mirrorAngle = mirrorAngle;
                save([wsArray.trials{k}.trackerFileName, '_WST.mat'], 'ws') % saving ws            
            end
        end
    end
end




cd(curr_d)

end