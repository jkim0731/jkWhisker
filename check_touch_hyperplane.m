clear
mouseName = 'JK056';
sessionName = 'S11';

wDir = 'D:\WhiskerVideo\';
sdir = [wDir,mouseName,sessionName];

cd(sdir)
hpFn = [mouseName, sessionName, '_touch_hp.mat'];
load(hpFn)

wsa = Whisker.WhiskerSignalTrialArray_2pad(sdir);
%%
close all
for i = 1:size(touch_hp,1)*size(touch_hp,2)
    figure('units','normalized','position', [0 0 1 1])
    angle = servo_distance_pair{i}(1);
    dist = servo_distance_pair{i}(2);
    tinds = find(cellfun(@(x) x.angle == angle && x.radialDistance == dist, wsa.trials));
    points = [];
    for j = 1 : length(tinds)
        frames = wsa.trials{tinds(j)}.poleUpFrames;
        points = [points, [wsa.trials{tinds(j)}.whiskerEdgeCoord(frames,1)'; wsa.trials{tinds(j)}.whiskerEdgeCoord(frames,2)'; ones(1,length(frames))*wsa.trials{tinds(j)}.apUpPosition]];
    end
    
    if psi1(i) > 90
        A = viewmtx(psi1(i),-90+psi2(i));
    else
        A = viewmtx(psi1(i),90-psi2(i));
    end
    points_4d = [points; ones(1,size(points,2))];
    points_2d = A*points_4d;
    points_2d = unique(round(points_2d(1:2,:)',2),'rows');
    
    th_4d1 = [touch_hp{i}(1,:) + hp_peaks{i}(1)     + 0;    touch_hp{i}(2:3,:);ones(1,size(touch_hp{i},2))];
    th_2d1 = A*th_4d1;
    th_2d1 = unique(th_2d1(1:2,:)','rows');
    th_4d2 = [touch_hp{i}(1,:) + hp_peaks{i}(2)     + 0;    touch_hp{i}(2:3,:);ones(1,size(touch_hp{i},2))];
    th_2d2 = A*th_4d2;
    th_2d2 = unique(th_2d2(1:2,:)','rows');
    scatter(points_2d(:,1),points_2d(:,2),'k.'), hold on, scatter(th_2d1(:,1), th_2d1(:,2),'r.'), scatter(th_2d2(:,1), th_2d2(:,2),'r.')
    title(['Angle = ', num2str(angle), ', Dist = ', num2str(dist)])
    
    
%     th_4d1 = [touch_hp{i}(1,:) + hp_peaks{i}(1)     + 1;    touch_hp{i}(2:3,:);ones(1,size(touch_hp{i},2))];
%     th_2d1 = A*th_4d1;
%     th_2d1 = unique(th_2d1(1:2,:)','rows');
%     th_4d2 = [touch_hp{i}(1,:) + hp_peaks{i}(2)     + 1;    touch_hp{i}(2:3,:);ones(1,size(touch_hp{i},2))];
%     th_2d2 = A*th_4d2;
%     th_2d2 = unique(th_2d2(1:2,:)','rows');
%     subplot(1,2,2),scatter(points_2d(:,1),points_2d(:,2),'k.'), hold on, scatter(th_2d1(:,1), th_2d1(:,2),'r.'), scatter(th_2d2(:,1), th_2d2(:,2),'r.')
%     title('+1 HP')    
    
end