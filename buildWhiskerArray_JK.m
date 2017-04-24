%% Setup whisker array builder 
base_angle = 21;

behavior_base_dir = 'Z:\Data\2p\soloData\';
whisker_base_dir = 'Z:\Data\Video\JK\';

sessions = {'AH0648S02','AH0650S02','AH0651S02','AH0652S04','AH0653S03'};
[mouseName, sessionName] = strtok(sessions{1},'S');
% mouse = 'AH0648';
% session = 'S02';

behavior_d = [behavior_base_dir mouseName '\'];
whisker_d = [whisker_base_dir mouseName sessionName '\'];
cd(behavior_d)
load('behavior.mat') % loading b of the mouse (all the sessions)

bb = cellfun(@(x) x.sessionName,b,'UniformOutput',false);
b_ind = find(cellfun(@(x) strcmp(x,sessionName),bb));
b_session = b{b_ind};

cd(whisker_d)

load_fn = [mouseName sessionName '_post.mat'];
load(load_fn); % loading errorlist
%%
filelist=dir([whisker_d '*.measurements']);

dirTrialNums=zeros(1,size(filelist,1));
% trialNums=[];  % enter which trial nums to process 

%%
% Assign the trial numbers to existing .measurements files in the directory
% NOTE : This assumes that the .measurements files have leading numbers
% corresponding to trial number in string positions 1:end-13 of the file
% name. These index numbers may need to be changed to match up to the
% numerical code of the trial number.  (2016/09/05 JK)

for i=1:length(filelist);
    dirTrialNums(i)=str2double(filelist(i).name(1:end-13)); % extract out the trial number from each measurements file present in directory
end
dirTrialNums = setdiff(dirTrialNums,errorlist);
trialNums = sort(dirTrialNums);
trialNums = trialNums(~isnan(trialNums));

includef=cell(size(trialNums,1),1);
for i = 1: length(trialNums)
    includef{i} = num2str(trialNums(i));
end

%% Optional section for cross correlating behavior and video trials, if you didn't pay attention to trial numbers
% vv = nan(max(dirTrialNums),1);
% 
% for i = 1:length(dirTrialNums) 
%     
%     
%     if ~isempty(find(dirTrialNums == i,1)); %create matrix showing trial location
%         fidx = find(dirTrialNums == i); %show index in matrix where trial is 
% 
%         tmp = load([filelist(fidx).name(1:end-13) '.bar']);
%          vv(i) = tmp(1,2);
%     else
%     end
%     
% end
% 
% gngThreshold = nanmean(vv) % for continuous pole postion with equal width go/nogo ranges, mean vv = the transition between go/nogo.
% vv2 = vv >= gngThreshold; % threshold for pole position on 1 side or the other
% vv3 = vv < gngThreshold;
% vvDiff = vv2-vv3;
% figure
% plot([vvDiff],'.')
% hold on
% plot(vv3,'ro')

%%  Build behavior number vector
% bv = zeros(max(b.trialNums),1);
% 
% bv(b.trialNums) = b.trialTypes*2-1; %1=GO and -1=NOGO and 0=NAN
% [c, lags] = xcorr(bv,vvDiff); %correlation between bv (behavioral GO and NOGO) and vv2 (video GO and NOGO trials)
% [mc, mx] = max(c); %find max c value and max x value
% lag_shift = lags(mx) %find lag in between both bv and vv2
% 
% for i = 1: length(filelist)
%     includef{i} = filelist(i).name(1:end-13);
% end
% 
% % self inputted values since tossed first 132 trials and know lag shift
% % lag_shift = -1
% 
% trialNums = dirTrialNums+lag_shift;  % correct the trial numbers
% 
% hold on
% plot(bv,'go')
% 
% % restrict processing to trials...
% startTrial = 1;
% endTrial = 99999999;
% includef = includef(trialNums >= startTrial & trialNums <=endTrial);
% trialNums =  trialNums(trialNums >= startTrial & trialNums <=endTrial);

%% Step 2 - Run with mask! for the first 8 trials and see how it looks
% 
% temp_tn = [50];
% 
% Whisker.makeAllDirectory_WhiskerTrial(d,[0 1],'mask', {[maskx(1,:);masky(1,:)],[maskx(2,:);masky(2,:)]},...
%     'trial_nums',trialNums(temp_tn),'include_files',includef(temp_tn),...
%     'barRadius',15.3,'faceSideInImage', 'bottom', 'framePeriodInSec',.0032,...
%     'imagePixelDimsXY',[vwidth vheight],'pxPerMm',26.23,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward')
% 
% % Whisker.makeAllDirectory_WhiskerSignalTrial(d,'include_files',includef(temp_tn),'polyRoiInPix',[40-20 40+20],'follicleExtrapDistInPix',26);
% Whisker.makeAllDirectory_WhiskerSignalTrial(d,'include_files',includef(temp_tn),'polyRoiInPix',[40-20 40+20]);
% Whisker.makeAllDirectory_WhiskerTrialLiteI(d,'include_files',includef(temp_tn),'r_in_mm',1,'calc_forces',false);
% wl = Whisker.WhiskerTrialLiteArray(d);
% 
% tid = [0 1]; % Set trajectory ID to view
% Whisker.viewdouble_WhiskerTrialLiteArray(wl,tid)

%% check mask - load .WST file first
% tp = [0 4];
% figure;ws.plot_fitted_whisker_time_projection(0,'k',tp), grid on
% hold on; ws.plot_fitted_whisker_time_projection(1,'k',tp)
% hold on; ws.plot_fitted_whisker_ROI_time_projection(0,'r',tp)
% hold on; ws.plot_fitted_whisker_ROI_time_projection(1,'r',tp)
% hold on; ws.plot_mask(0,'g',tp)
% hold on; ws.plot_mask(1,'g',tp)

%% Step 3 - run everything
% select matching files
%tmp = cellfun(@(x)str2num(x(15:end)),includef);
%incf_idx = find(tmp>= 12 & tmp <=85);

%% Make whisker-pole touch space for each type of trial, from 10 randomly selected trials (of each type)
% Currently, only dealing with 4 types of trials: 'rc', 'rf', 'lc', 'lf'
% Should make something different for straight pole touch in S00. 
% 2017/04/11 JK
trial_types = {'rc', 'rf', 'lc', 'lf'};
tt_ind = cell(1,length(trial_types));
wl_array = cell(1,length(trial_types));
touch_points = cell(1,length(trial_types));
touch_hp = cell(1,length(trial_type)); % touch hyperplanes
%%
for trial_type_num = 2 : length(trial_types)    
% trial_type_num = 1
    tt_ind{trial_type_num} = find(cellfun(@(x) strcmp(x.trialType,trial_types{trial_type_num}),b_session.trials));
%     if length(tt_ind{i}) > 10
%         tt_ind{i} = sort(randsample(tt_ind{i},10));
%     end
    temp_files = cell(length(tt_ind{trial_type_num}),1);
    for j = 1 : length(tt_ind{trial_type_num})
        temp_files{j} = num2str(tt_ind{trial_type_num}(j));
    end
%     Whisker.makeAllDirectory_WhiskerSignalTrial(whisker_d,'include_files',temp_files,'polyRoiInPix',[20 80],'pole_available_timepoints',[550:);
    Whisker.makeAllDirectory_WhiskerSignalTrial_2pad(whisker_d,'include_files',temp_files,'polyRoiInPix',[20 80]);
    Whisker.makeAllDirectory_WhiskerTrialLiteI(whisker_d,'include_files',temp_files,'r_in_mm',2,'calc_forces',false,'behavior',b_session);
    wl = Whisker.WhiskerTrialLiteArray(whisker_d,'include_files',temp_files);

%     tid = [0 1]; % Set trajectory ID to view
%     Whisker.viewdouble_WhiskerTrialLiteArray(wl,tid)
    wl_array{trial_type_num} = wl;
end

%%
for trial_type_num = 4
    intersect_3d = [];
    wl = wl_array{trial_type_num};
    figure, hold all
    for tnum = 1 : length(wl.trials)
        top_ind = find(wl.trials{tnum}.intersect_coord(:,1) > 50);
        front_ind = find(wl.trials{tnum}.intersect_coord(:,2) > 50);
        intersect_ind = intersect(wl.trials{tnum}.pole_available_timepoints,intersect(top_ind,front_ind));
        try
            plot3(wl.trials{tnum}.intersect_coord(intersect_ind,1), wl.trials{tnum}.intersect_coord(intersect_ind,2), ones(1,length(intersect_ind))*wl.trials{tnum}.pole_pos, 'k.', 'MarkerSize', 3)
            intersect_3d = [intersect_3d; wl.trials{tnum}.intersect_coord(intersect_ind,1), wl.trials{tnum}.intersect_coord(intersect_ind,2), ones(length(intersect_ind),1)*wl.trials{tnum}.pole_pos];
        catch
            fprintf('Skipping trial #%d because of index problems',tnum);        
        end
    end
    title(wl.trials{1}.trial_type), xlabel('Top-view intersection coord'), ylabel('Front-view intersection coord'), zlabel('Pole position')
end

%% Calculate psi1 % takes ~ 15 sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Manual selection
ind_opt = 1; % optimal peak index. Starting from 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_zeros = 0;
angle_steps = 0:180;
stds_pre = zeros(length(angle_steps),1);
for i = 1 : length(angle_steps)
    A = viewmtx(angle_steps(i),0);
    x4d = [intersect_3d, ones(size(intersect_3d,1),1)]';
    x2d = A*x4d;
    x2d_pix = [floor(x2d(1,:));floor(x2d(2,:)/100)];
    x2d_dim = [max(x2d_pix(2,:)) - min(x2d_pix(2,:)) + 1, max(x2d_pix(1,:)) - min(x2d_pix(1,:)) + 1];
    x2d_proj = zeros(x2d_dim);
%     for i = 1:x2d_dim(1) 
%         temp1d_ind = find(x2d_pix(2,:) == (max(x2d_pix(2,:)) - i + 1));
%         for j = 1:x2d_dim(2)
%             x2d_proj(i,j) = numel(find(x2d_pix(1,temp1d_ind) == (min(x2d_pix(1,:)) + j - 1)));
%         end
%     end
    j_offset = min(x2d_pix(1,:)) - 1;
    i_offset = min(x2d_pix(2,:)) - 1;
    for j = 1 : length(x2d_pix)
        x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) = x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) + 1;
    end
    stds_pre(i) = std(x2d_proj(find(x2d_proj(:))));
%     temp_std = std(x2d_proj(:));
%     temp_std = sum(x2d_proj(:)==0);    
end
[P, I] = findpeaks(smooth(smooth(stds_pre)));
[~, I2] = sort(P, 'descend');
max_psi1_pre = angle_steps(I(I2(ind_opt)));

figure, plot(1:length(stds_pre), smooth(smooth(stds_pre,5))), hold on, plot(I(I2(ind_opt)),stds_pre(I(I2(ind_opt))),'ro')

max_std = 0;
psi1 = 0;
angle_steps = max_psi1_pre-4.99:0.01:max_psi1_pre+4.99;
stds = zeros(length(angle_steps),1);
x2d_final = zeros(size(x2d_proj));
for i = 1:length(angle_steps)
    A = viewmtx(angle_steps(i),0);
    x4d = [intersect_3d, ones(size(intersect_3d,1),1)]';
    x2d = A*x4d;
    x2d_pix = [floor(x2d(1,:));floor(x2d(2,:)/100)];
    x2d_dim = [max(x2d_pix(2,:)) - min(x2d_pix(2,:)) + 1, max(x2d_pix(1,:)) - min(x2d_pix(1,:)) + 1];
    x2d_proj = zeros(x2d_dim);

    j_offset = min(x2d_pix(1,:)) - 1;
    i_offset = min(x2d_pix(2,:)) - 1;
    for j = 1 : length(x2d_pix)
        x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) = x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) + 1;
    end
    temp_std = std(x2d_proj(find(x2d_proj(:))));
%     temp_std = std(x2d_proj(:));
%     temp_std = sum(x2d_proj(:)==0);
    stds(i) = temp_std;
    if temp_std > max_std
        max_std = temp_std;
        psi1 = angle_steps(i);
        x2d_final = x2d_proj;
    end
end
psi1 = psi1-90;
figure, plot(1:length(stds), stds)

A = viewmtx(psi1+90,0);
% A = viewmtx(44,0);
x4d = [intersect_3d, ones(size(intersect_3d,1),1)]';
x2d = A*x4d;
figure, plot(x2d(1,:), x2d(2,:),'k.', 'MarkerSize',3)
%% Draw polygon to select regions for radon transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Manual selection
x2d_flip = flip(x2d_final,1);
BW = roipoly(x2d_flip);
figure, imshow(BW)
x2d_edge = x2d_flip .* BW;
figure, imagesc(x2d_edge)

%% Calculate psi2 % ~ 3 sec

theta = 0:0.01:180;
R = radon(x2d_edge, theta);
[~, max_ind] = max(max(R));
psi2 = (max_ind-1)*0.01;
psi2 = atand(tand(psi2)/100); % psi2 adjusted because it was calculated with pole position divided by 100

figure, imagesc(R), xlabel('Angle = 0:0.01:180'), ylabel('Projected values') 
figure, imagesc(x2d_flip, [0 1000]), hold on, line([54 54+tand((max_ind-1)*0.01)*size(x2d_flip,1)],[1 size(x2d_flip,1)], 'LineWidth', 3, 'Color', [1 1 1])
line([69 69+tand((max_ind-1)*0.01)*size(x2d_flip,1)],[1 size(x2d_flip,1)], 'LineWidth', 3, 'Color', [1 1 1])

%% Calculate touch hyperplanes
xmin = -200; xmax = 200; zmin_data = min(intersect_3d(:,3)); zmax_data = max(intersect_3d(:,3)); zdiff = zmax_data - zmin_data; zmax = zmax_data + zdiff; zmin = zmin_data - zdiff;
z = zmin:zmax;
xyz = zeros((length(z))*(xmax-xmin+1),3);
for i = xmin:xmax
    xyz((i-xmin)*length(z)+1 : (i-xmin+1)*length(z),:) = [ones(length(z),1)*i, zeros(length(z),1), z'];
end

% xyz = xyz';
% xyz(:,xyz(1,:) < 0) = [];
% xyz(:,xyz(1,:) > xmax) = [];
% xyz(:,xyz(2,:) < 0) = [];
% xyz(:,xyz(2,:) > xmax) = [];
% xyz(:,xyz(3,:) < zmin_data) = [];
% xyz(:,xyz(3,:) > zmax_data) = [];
% 
% figure, plot3(xyz(1,:), xyz(2,:), xyz(3,:), 'r.', 'MarkerSize',3), hold on,
% plot3(intersect_3d(:,1),intersect_3d(:,2), intersect_3d(:,3),'k.', 'MarkerSize',3), xlabel('top'), ylabel('front'), zlabel('pos')

[xyz_psi1, ~, ~] = AxelRot(xyz',psi1,[0 0 1], 0); % rotate psi1 degrees counterclockwise around z axis

% xyz_psi1(:,xyz_psi1(1,:) < 0) = [];
% xyz_psi1(:,xyz_psi1(1,:) > xmax) = [];
% xyz_psi1(:,xyz_psi1(2,:) < 0) = [];
% xyz_psi1(:,xyz_psi1(2,:) > xmax) = [];
% xyz_psi1(:,xyz_psi1(3,:) < zmin_data) = [];
% xyz_psi1(:,xyz_psi1(3,:) > zmax_data) = [];
% 
% figure, plot3(xyz_psi1(1,:), xyz_psi1(2,:), xyz_psi1(3,:), 'r.', 'MarkerSize',3), hold on,
% plot3(intersect_3d(:,1),intersect_3d(:,2), intersect_3d(:,3),'k.', 'MarkerSize',3), xlabel('top'), ylabel('front'), zlabel('pos')

zcenter = floor(mean([zmax_data, zmin_data]));
x0 = [0 0 zcenter];
u = [1 tand(psi1) 0];
[xyz_psi2, ~, ~] = AxelRot(xyz_psi1, psi2, u, x0); 
% xyz_psi2(:,xyz_psi2(1,:) <=0) = [];
% xyz_psi2(:,xyz_psi2(1,:) > xmax) = [];
% xyz_psi2(:,xyz_psi2(2,:) <= 0) = [];
% xyz_psi2(:,xyz_psi2(2,:) > xmax) = [];
xyz_psi2(:,xyz_psi2(3,:) < zmin_data) = [];
xyz_psi2(:,xyz_psi2(3,:) > zmax_data) = [];

figure, plot3(intersect_3d(:,1),intersect_3d(:,2), intersect_3d(:,3),'k.', 'MarkerSize',3), xlabel('top'), ylabel('front'), zlabel('pos'), hold on
plot3(xyz_psi2(1,:), xyz_psi2(2,:), xyz_psi2(3,:), 'r.', 'MarkerSize',3)

%% ~ 8 min
%%%%%%%%%%%%%%%%%%%%%% manual selection
steps = 150:300;
%%%%%%%%%%%%%%%%%%%%%% try as short as possible to reduce time next step


if abs(tand(psi1+90)) > 1
    transvec = [1/tand(psi1+90) 1 0]; % translation vector, moving xyz_psi2 along its normal vector after projecting to x-y plane.
else
    transvect = [1 tand(psi1+90) 0];
end

intersect_pix = round(intersect_3d);

num_points_in_hp = zeros(max(steps),1);
for i = steps % this is time consuming...
    hp = round(xyz_psi2);
    hp(1,:) = hp(1,:)+i;
    [~,points,~] = intersect(intersect_pix, hp','rows');
    num_points_in_hp(i) = numel(points);
end

figure, plot(steps,num_points_in_hp(steps), 'k-', 'LineWidth', 3), xlabel('translocation (pix)'), ylabel('# intersection')

%%
touch_hp{trial_type_num} = xyz_psi2; % Don't round them! (at least at this saving process)

%end
%% Identifying touch frames in each type of trial
% trial_types = {'rc', 'rf', 'lc', 'lf'};
% for i = 1 : 4
%     max_trial_num(i) = max(cellfun(@(x) x.trialNum,wl_array{i}.trials));
% end
% max_trial_num = max(max_trial_num);
% touch_frames = cell(1,max_trial_num); % this is based just on the trial num, to make it easier for later comparison with manual inspection of the touches (obvious touches and obvious non-touches)
% for i = 1 : 4
%     for j = 1 : length(wl_array{i}.trials)
%         t = wl_array{i}.trials{j};
%         switch t.trial_type
%             case trial_types{1}
%                 touch_frames{t.trialNum} = intersect(touch_hp{1}
%         end
%     end
% end

%%
trialnum = 70;
temp_ta = round(wl_array{2}.trials{trialnum}.intersect_coord);

temp_touch_hp = touch_hp{2};
temp_touch_hp(1,:) = temp_touch_hp(1,:) + 46;
temp_touch_hp = sortrows(temp_touch_hp',3)';
temp_ind = find(temp_touch_hp(3,:) == wl_array{1}.trials{trialnum}.pole_pos);
temp_touch_hp_ = temp_touch_hp(1:2,temp_ind);
% figure, plot(temp_touch_hp_(1,:), temp_touch_hp_(2,:),'r.'), hold on, plot(temp_ta(550:end-200,1),temp_ta(550:end-200,2),'k.')

temp_touch_hp_ = [temp_touch_hp_(1,:), temp_touch_hp_(1,:) + 1, temp_touch_hp_(1,:) - 1, temp_touch_hp_(1,:) - 2; ...
    temp_touch_hp_(2,:),temp_touch_hp_(2,:),temp_touch_hp_(2,:),temp_touch_hp_(2,:)];
temp_touch_hp_ = temp_touch_hp_';
temp_touch_hp_ = unique(temp_touch_hp_,'rows');
touch_frame = find(ismember(temp_ta,temp_touch_hp_,'rows')==1);

wl_array{2}.trials{trialnum}.trackerFileName

figure,  plot(temp_touch_hp_(:,1),temp_touch_hp_(:,2),'r.'), hold on, plot(temp_ta(550:end-200,1),temp_ta(550:end-200,2),'k.')

%%
tt_ind{1}
%%
% touch_points{1} = cell(1,length(wl.trialNums));
trial_nums = [138,216,344,414,418,498,705,764,787]; % 322
touch_points{1} = cell(1,length(trial_nums));
touch_points{1}{1} = [592, 615, 649, 675, 689, 707, 722, 1121, 1146, 1163];
touch_points{1}{2} = [1130, 1143, 1153];
touch_points{1}{3} = [682, 772, 780, 856, 891, 905, 919, 944, 1014, 1020, 1032, 1146, 1275, 1344,1345];
touch_points{1}{4} = [596,980,1005,1223,2020,2085,2092,2419, 2434,2850,3188,3221,3340,3496,3508,3526];
touch_points{1}{5} = [1110,1130];
touch_points{1}{6} = [534:538];
touch_points{1}{7} = [1439,1559,1609,1625,2853,2858,2872,2899,2954];
touch_points{1}{8} = [612,740, 746,843,873,911,1095,1134,1190,1205];
touch_points{1}{9} = [589,605,664,679,716,725,732,885,917,1292,1536,1574,1592,1677];
% touch_points{1}{10} = [];

%%
temp_files = cell(1,length(trial_nums));
for i = 1 : length(trial_nums)
    temp_files{i} = num2str(trial_nums(i));
end
    Whisker.makeAllDirectory_WhiskerSignalTrial(whisker_d,'include_files',temp_files,'polyRoiInPix',[20 80]);
    Whisker.makeAllDirectory_WhiskerTrialLiteI(whisker_d,'include_files',temp_files,'r_in_mm',2,'calc_forces',false,'behavior',b_session);
    wl = Whisker.WhiskerTrialLiteArray(whisker_d,'include_files',temp_files);
%%
temp_trial = 1;
figure,
subplot(311), plot(1:length(wl.trials{temp_trial}.intersect_coord),smooth(wl.trials{temp_trial}.intersect_coord(:,1)),'k.'), hold on,
plot(1:length(wl.trials{temp_trial}.intersect_coord),smooth(wl.trials{temp_trial}.intersect_coord(:,2)),'b.'),
% plot(touch_points{1}{temp_trial}+1,wl.trials{temp_trial}.intersect_coord(touch_points{1}{temp_trial}+1,1),'r.')
% plot(touch_points{1}{temp_trial}+1,wl.trials{temp_trial}.intersect_coord(touch_points{1}{temp_trial}+1,2),'r.')
subplot(312), plot(1:length(wl.trials{temp_trial}.intersect_coord),smooth(wl.trials{temp_trial}.deltaKappa{1}),'k.'), hold on,
plot(1:length(wl.trials{temp_trial}.intersect_coord),smooth(wl.trials{temp_trial}.deltaKappa{2}),'b.')
subplot(313), plot(1:length(wl.trials{temp_trial}.intersect_coord),smooth(wl.trials{temp_trial}.thetaAtBase{1}),'k.'), hold on,
plot(1:length(wl.trials{temp_trial}.intersect_coord),smooth(wl.trials{temp_trial}.thetaAtBase{2}),'b.')

%%
figure, hold all
for trial_num = 1 : length(trial_nums)
    plot3(wl.trials{trial_num}.intersect_coord(touch_points{1}{trial_num}+1,1)',wl.trials{trial_num}.intersect_coord(touch_points{1}{trial_num}+1,2)',ones(length(touch_points{1}{trial_num}),1)*wl.trials{trial_num}.pole_pos);
end
%%
list = [1,3,4,7,8,9];
edge_fit = cell(size(list));
for i = 1:length(list)
    edge_fit{i} = polyfit(wl.trials{list(i)}.intersect_coord(touch_points{1}{list(i)}+1,1)',wl.trials{list(i)}.intersect_coord(touch_points{1}{list(i)}+1,2)',1);
end
%%
pole_pos = zeros(size(list));
for i = 1 : length(list)
    pole_pos(i) = wl.trials{list(i)}.pole_pos;
end
figure, plot(pole_pos,cellfun(@(x) x(2),edge_fit))

%%
Whisker.makeAllDirectory_WhiskerTrial(whisker_d,[0 1],'mask', {[maskx(1,:);masky(1,:)],[maskx(2,:);masky(2,:)]},...
    'trial_nums',trialNums,'include_files',includef,...
    'barRadius',15.3,'faceSideInImage', 'top', 'framePeriodInSec',.0032,...
    'imagePixelDimsXY',[width height],'pxPerMm',26.23,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','leftward')

Whisker.makeAllDirectory_WhiskerSignalTrial(whisker_d,'include_files',includef,'polyRoiInPix',[20 80]);
Whisker.makeAllDirectory_WhiskerTrialLiteI(whisker_d,'include_files',includef,'r_in_mm',3,'calc_forces',true,'whisker_radius_at_base', 36.5,'whisker_length', 18,'baseline_time_or_kappa_value',0);
wl = Whisker.WhiskerTrialLiteArray(whisker_d);
save([d mouseName sessionName '-WTLIA.mat'],'wl');

tid = 0; % Set trajectory ID to view
Whisker.view_WhiskerTrialLiteArray(wl,tid)
