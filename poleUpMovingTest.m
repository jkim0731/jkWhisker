close all
% sessions = {[4,19,22],[3,16,17],[3,21,22],[1,17,18,91],[7],[2],[1,22:25],[3]};
mouse = 'JK025';
session = 'S18';
dirBase = 'E:\WhiskerVideo\';
dirName = [dirBase, mouse, session];
cd(dirName)
flist = dir('*_WT.mat');
fnums = zeros(length(flist),1);
for i = 1 : length(flist)
    fnums(i) = str2double(strtok(flist(i).name, '_'));
end
fnums = sort(fnums);
wtArray = cell(length(fnums),1);
for i = 1 : length(wtArray)
    disp(['Loading ', num2str(fnums(i)), '_WT, ', num2str(i) , '/', num2str(length(wtArray))])
    fname = sprintf('%d_WT.mat',fnums(i));
    load(fname); % loading w
    wtArray{i} = w;
end
poleUpLength = cellfun(@(x) length(x.poleUpFrames), wtArray);
poleMovingLength = cellfun(@(x) length(x.poleMovingFrames), wtArray);
noPoleLength = cellfun(@(x) x.nof - length(x.poleUpFrames) - length(x.poleMovingFrames), wtArray);
%
figure,
oonum = find(cellfun(@(x) strcmp(x.trialType,'oo'), wtArray));
subplot(311), plot(poleUpLength), hold on, plot(oonum,poleUpLength(oonum), 'r.')
subplot(312), plot(poleMovingLength), hold on, plot(oonum,poleMovingLength(oonum), 'r.')
subplot(313), plot(noPoleLength), hold on, plot(oonum,noPoleLength(oonum), 'r.')
%
figure,
poleAxesUpX = [];
poleAxesUpY = [];
poleAxesUpX90 = [];
poleAxesUpY90 = [];
angles = unique(cellfun(@(x) x.angle, wtArray));
rds = unique(cellfun(@(x) x.radialDistance, wtArray));
rds = rds(find(rds));
for ai = 1 : length(angles)
    for ri = 1 : length(rds)
        tn = find(cellfun(@(x) x.angle == angles(ai) && x.radialDistance == rds(ri), wtArray));
        imshow(wtArray{tn(10)}.binvavg), hold on, 
        for j = 1 : length(tn)
            plot(wtArray{tn(j)}.poleAxesUp{1}(1,:), wtArray{tn(j)}.poleAxesUp{1}(2,:), 'r.'), 
            plot(wtArray{tn(j)}.poleAxesUp{2}(1,:), wtArray{tn(j)}.poleAxesUp{2}(2,:), 'b.')
        end
        if wtArray{tn(j)}.angle == 90
            poleAxesUpX90 = wtArray{tn(1)}.poleAxesUp{1}(1,:);
            poleAxesUpY90 = wtArray{tn(1)}.poleAxesUp{1}(2,:);
        else
            poleAxesUpX = [poleAxesUpX; wtArray{tn(1)}.poleAxesUp{1}(1,:)];
            poleAxesUpY = [poleAxesUpY; wtArray{tn(1)}.poleAxesUp{1}(2,:)];
        end
        waitforbuttonpress
%         hold off, imshow(wtArray{tn(10)}.binvavg), hold on, 
%         for j = 1 : length(tn)
%             plot(wtArray{tn(j)}.poleAxesUp{1}(1,:), wtArray{tn(j)}.poleAxesUp{1}(2,:), 'r.'), 
%             plot(wtArray{tn(j)}.poleAxesUp{2}(1,:), wtArray{tn(j)}.poleAxesUp{2}(2,:), 'b.')
%         end
%         waitforbuttonpress
%         hold off
    end
end
figure, hold on
for i = 1 : size(poleAxesUpX,1)
    plot(poleAxesUpX(i,:), poleAxesUpY(i,:), 'k-')
end
plot(poleAxesUpX90, poleAxesUpY90, 'r-')
%%

for i = 1 : length(wsArray.trials)
    for j = 1 : length(wsArray.trials{i}.poleMovingFrames)
        if sum(wsArray.trials{i}.poleAxesUp{1}(1,:) - wsArray.trials{i}.poleAxesMoving{1}(1,:)) + sum(wsArray.trials{i}.poleAxesUp{1}(2,:) - wsArray.trials{i}.poleAxesMoving{1}(2,:)) == 0
            disp('right')
        end
    end
end

%%
% [6,7,25,26,48,81,86,88,89,98,111,121,125,136,142,146,169,171,175,178,189,196,198,211,212,220,222,226,232,234,243,253,255,267,269,274,289,296,309,313,323,326,338,355,357,365,374,381,382,387,395,403,414,423]
%%
figure
plotnum = setdiff(1:length(poleUpLength),oonum);
subplot(311), plot(poleUpLength(plotnum))
subplot(312), plot(poleMovingLength(plotnum))
subplot(313), plot(noPoleLength(plotnum))

%%

olnum = 85; % outlier number

figure, 
if ~isempty(wtArray{olnum}.topPix)
    subplot(311)
    plot(wtArray{olnum}.topPix(:,1)), 
end
if ~isempty(wtArray{olnum}.frontPix)
    subplot(312)
    plot(wtArray{olnum}.frontPix(:,2)), 
end

title([wtArray{olnum}.trackerFileName, '  ', wtArray{olnum}.trialType])
subplot(313), imshow(wtArray{olnum}.binvavg), hold on, 
plot(wtArray{olnum}.poleAxesUp{1}(1,:), wtArray{olnum}.poleAxesUp{1}(2,:), 'r.'), 
plot(wtArray{olnum}.poleAxesUp{2}(1,:), wtArray{olnum}.poleAxesUp{2}(2,:), 'b.')
%%
% wtArray{olnum}.
%%
figure, plot(wsArray.trials{olnum}.topPix(:,2)), title(wsArray.trials{olnum}.trackerFileName)
%%
figure, histogram(abs(diff(wsArray.trials{olnum}.topPix(:,1))))
%%
prctile(wsArray.trials{olnum}.topPix(:,1),10)