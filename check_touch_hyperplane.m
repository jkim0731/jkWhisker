clear
mouseName = 'JK030';
sessionName = 'S08';

wDir = 'Y:\Whiskernas\JK\whisker\tracked\';
sdir = [wDir,mouseName,sessionName];

cd(sdir)
hpFn = [mouseName, sessionName, '_touch_hp.mat'];
load(hpFn)

wsa = Whisker.WhiskerSignalTrialArray_2pad(sdir);
%%
close all
% for sdi = 1 : length(servo_distance_pair)
for sdi = 2
    
    figure('units','normalized','position', [0 0 1 1])
    angle = servo_distance_pair{sdi}(1);
    dist = servo_distance_pair{sdi}(2);
    tinds = find(cellfun(@(x) x.angle == angle && x.radialDistance == dist, wsa.trials));
    points = [];
    for j = 1 : length(tinds)
        frames = wsa.trials{tinds(j)}.poleUpFrames;
        points = [points, [wsa.trials{tinds(j)}.whiskerEdgeCoord(frames,1)'; wsa.trials{tinds(j)}.whiskerEdgeCoord(frames,2)'; ones(1,length(frames))*wsa.trials{tinds(j)}.apUpPosition]];
    end
    
    if psi1(sdi) > 90
        A = viewmtx(psi1(sdi),-90+psi2(sdi));
    else
        A = viewmtx(psi1(sdi),90-psi2(sdi));
    end
    points_4d = [points; ones(1,size(points,2))];
    points_2d = A*points_4d;
    points_2d = unique(round(points_2d(1:2,:)',2),'rows');
    
    th_4d1 = [touch_hp{sdi}(1,:) + hp_peaks{sdi}(1)     + 0;    touch_hp{sdi}(2:3,:);ones(1,size(touch_hp{sdi},2))];
    th_2d1 = A*th_4d1;
    th_2d1 = unique(th_2d1(1:2,:)','rows');
    th_4d2 = [touch_hp{sdi}(1,:) + hp_peaks{sdi}(2)     + 0;    touch_hp{sdi}(2:3,:);ones(1,size(touch_hp{sdi},2))];
    th_2d2 = A*th_4d2;
    th_2d2 = unique(th_2d2(1:2,:)','rows');
    scatter(points_2d(:,1),points_2d(:,2),'k.'), hold on, scatter(th_2d1(:,1), th_2d1(:,2),'r.'), scatter(th_2d2(:,1), th_2d2(:,2),'r.')
    title(['Angle = ', num2str(angle), ', Dist = ', num2str(dist)])
    figure, 
    
    plot3(points(1,:), points(2,:), points(3,:), 'k.'), hold on
    plot3(touch_hp{sdi}(1,:) + hp_peaks{sdi}(1), touch_hp{sdi}(2,:), touch_hp{sdi}(3,:), 'r-')
    plot3(touch_hp{sdi}(1,:) + hp_peaks{sdi}(2), touch_hp{sdi}(2,:), touch_hp{sdi}(3,:), 'r-')
end

%%
% for sdi = 1 : length(servo_distance_pair)
for sdi = 1
    inds = find(cellfun(@(x) x.angle == servo_distance_pair{sdi}(1) && x.radialDistance == servo_distance_pair{sdi}(2), wsa.trials));
    points3d = [];
    for ii = 1 : length(inds)
        wec = wsa.trials{inds(ii)}.whiskerEdgeCoord;
        wecind = find(isfinite(sum(wec,2)));
        points3d = [points3d; wec(wecind,:), wsa.trials{inds(ii)}.apUpPosition * ones(length(wecind),1)];
    end
    
    figure, 
    plot3(points3d(:,1), points3d(:,2), points3d(:,3), 'k.'), hold on
    plot3(touch_hp{sdi}(1,:) + hp_peaks{sdi}(1), touch_hp{sdi}(2,:), touch_hp{sdi}(3,:), 'r-')
    plot3(touch_hp{sdi}(1,:) + hp_peaks{sdi}(2), touch_hp{sdi}(2,:), touch_hp{sdi}(3,:), 'r-')
end