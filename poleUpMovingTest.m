close all
mouse = 'JK041';
session = 'S03';
dirBase = 'E:\WhiskerVideo\';
dirName = [dirBase, mouse, session];
cd(dirName)
load('*_touch_hp.mat')
wsArray = Whisker.WhiskerSignalTrialArray_2pad(dirName);
poleUpLength = cellfun(@(x) length(x.poleUpFrames), wsArray.trials);
poleMovingLength = cellfun(@(x) length(x.poleMovingFrames), wsArray.trials);
noPoleLength = cellfun(@(x) x.nof - length(x.poleUpFrames) - length(x.poleMovingFrames), wsArray.trials);

figure,
oonum = find(cellfun(@(x) strcmp(x.trialType,'oo'), wsArray.trials));
subplot(311), plot(poleUpLength), hold on, plot(oonum,poleUpLength(oonum), 'r.')
subplot(312), plot(poleMovingLength), hold on, plot(oonum,poleMovingLength(oonum), 'r.')
subplot(313), plot(noPoleLength), hold on, plot(oonum,noPoleLength(oonum), 'r.')


for i = 1 : length(servo_distance_pair)
    unique(cellfun(@(x) sum([x.angle, x.radialDistance] == , wsArray.trials));


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
    plot(wsArray.trials{olnum}.topPix(:,1)), 
end

title([wsArray.trials{olnum}.trackerFileName, '  ', wsArray.trials{olnum}.trialType])
%%
figure, plot(wsArray.trials{olnum}.topPix(:,2)), title(wsArray.trials{olnum}.trackerFileName)
%%
figure, histogram(abs(diff(wsArray.trials{olnum}.topPix(:,1))))
%%
prctile(wsArray.trials{olnum}.topPix(:,1),10)