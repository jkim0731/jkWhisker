close all
% sessions = {[4,19,22],[3,16,17],[3,21,22],[1,17,18,91],[7],[2],[1,22:25],[3]};
mouse = 'JK041';
session = 'S03';
dirBase = 'E:\WhiskerVideo\';
dirName = [dirBase, mouse, session];
cd(dirName)
flist = dir('*_touch_hp.mat');
load(flist(1).name)
wsArray = Whisker.WhiskerSignalTrialArray_2pad(dirName);
wtArray = cell(length(wsArray.trials),1);
for i = 1 : length(wtArray)
    disp(['Loading ', wsArray.trials{i}.trackerFileName, '_WT, ', num2str(i) , '/', num2str(length(wtArray))])
    wtArray{i} = load([wsArray.trials{i}.trackerFileName, '_WT.mat']);
end
poleUpLength = cellfun(@(x) length(x.poleUpFrames), wsArray.trials);
poleMovingLength = cellfun(@(x) length(x.poleMovingFrames), wsArray.trials);
noPoleLength = cellfun(@(x) x.nof - length(x.poleUpFrames) - length(x.poleMovingFrames), wsArray.trials);

figure,
oonum = find(cellfun(@(x) strcmp(x.trialType,'oo'), wsArray.trials));
subplot(311), plot(poleUpLength), hold on, plot(oonum,poleUpLength(oonum), 'r.')
subplot(312), plot(poleMovingLength), hold on, plot(oonum,poleMovingLength(oonum), 'r.')
subplot(313), plot(noPoleLength), hold on, plot(oonum,noPoleLength(oonum), 'r.')
%
figure,
poleAxesUpX = [];
poleAxesUpY = [];
poleAxesUpX90 = [];
poleAxesUpY90 = [];
for i = 1 : size(servo_distance_pair,1)*size(servo_distance_pair,2)
% for i = 1
    tn = find(cellfun(@(x) sum([x.angle, x.radialDistance] == servo_distance_pair{i}) == 2, wsArray.trials));
    imshow(wsArray.trials{tn(10)}.binvavg), hold on, 
    for j = 1 : length(tn)
        plot(wsArray.trials{tn(j)}.poleAxesUp{1}(1,:), wsArray.trials{tn(j)}.poleAxesUp{1}(2,:), 'r.'), 
        plot(wsArray.trials{tn(j)}.poleAxesUp{2}(1,:), wsArray.trials{tn(j)}.poleAxesUp{2}(2,:), 'b.')
    end
    if wsArray.trials{tn(j)}.angle == 90
        poleAxesUpX90 = wsArray.trials{tn(1)}.poleAxesUp{1}(1,:);
        poleAxesUpY90 = wsArray.trials{tn(1)}.poleAxesUp{1}(2,:);
    else
        poleAxesUpX = [poleAxesUpX; wsArray.trials{tn(1)}.poleAxesUp{1}(1,:)];
        poleAxesUpY = [poleAxesUpY; wsArray.trials{tn(1)}.poleAxesUp{1}(2,:)];
    end
    waitforbuttonpress
    hold off, imshow(wsArray.trials{tn(10)}.binvavg), hold on, 
    for j = 1 : length(tn)
        plot(wtArray{tn(j)}.w.poleAxesUp{1}(1,:), wtArray{tn(j)}.w.poleAxesUp{1}(2,:), 'r.'), 
        plot(wtArray{tn(j)}.w.poleAxesUp{2}(1,:), wtArray{tn(j)}.w.poleAxesUp{2}(2,:), 'b.')
    end
    waitforbuttonpress
    hold off
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

olnum = 18; % outlier number

figure, 
if ~isempty(wsArray.trials{olnum}.topPix)
    subplot(121)
    plot(wsArray.trials{olnum}.topPix(:,1)), 
end

title([wsArray.trials{olnum}.trackerFileName, '  ', wsArray.trials{olnum}.trialType])
subplot(122), imshow(wsArray.trials{olnum}.binvavg), hold on, 
plot(wsArray.trials{olnum}.poleAxesUp{1}(1,:), wsArray.trials{olnum}.poleAxesUp{1}(2,:), 'r.'), 
plot(wsArray.trials{olnum}.poleAxesUp{2}(1,:), wsArray.trials{olnum}.poleAxesUp{2}(2,:), 'b.')
%%
figure, plot(wsArray.trials{olnum}.topPix(:,2)), title(wsArray.trials{olnum}.trackerFileName)
%%
figure, histogram(abs(diff(wsArray.trials{olnum}.topPix(:,1))))
%%
prctile(wsArray.trials{olnum}.topPix(:,1),10)