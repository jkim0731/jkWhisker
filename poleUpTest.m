function poleUpTest(mouse, session, dirBase)
close all
% sessions = {[4,19,22],[3,16,17],[3,21,22],[1,17,18,91],[7],[2],[1,22:25],[3]};
% mouse = 'JK053';
% session = 'S04';
% dirBase = 'L:\tracked\';
dirName = [dirBase, mouse, session];
if exist(dirName, 'dir')
    cd(dirName)
else
    disp(['No directory named as ', dirName])
    return
end
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
%%
figure,
oonum = find(cellfun(@(x) strcmp(x.trialType,'oo'), wtArray));
subplot(311), plot(poleUpLength), hold on, plot(oonum,poleUpLength(oonum), 'r.')
subplot(312), plot(poleMovingLength), hold on, plot(oonum,poleMovingLength(oonum), 'r.')
subplot(313), plot(noPoleLength), hold on, plot(oonum,noPoleLength(oonum), 'r.')
%%
figure,
poleAxesUpX = [];
poleAxesUpY = [];
poleAxesUpX90 = [];
poleAxesUpY90 = [];
angles = unique(cellfun(@(x) x.angle, wtArray));
rds = unique(cellfun(@(x) x.radialDistance, wtArray));
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
        hold off, imshow(wtArray{tn(10)}.binvavg), hold on, 
        for j = 1 : length(tn)
            plot(wtArray{tn(j)}.poleAxesUp{1}(1,:), wtArray{tn(j)}.poleAxesUp{1}(2,:), 'r.'), 
            plot(wtArray{tn(j)}.poleAxesUp{2}(1,:), wtArray{tn(j)}.poleAxesUp{2}(2,:), 'b.')
        end
        waitforbuttonpress
        hold off
    end
end
figure, hold on
for i = 1 : size(poleAxesUpX,1)
    plot(poleAxesUpX(i,:), poleAxesUpY(i,:), 'k-')
end
plot(poleAxesUpX90, poleAxesUpY90, 'r-')