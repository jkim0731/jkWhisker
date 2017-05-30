

%% Setup whisker array builder 
cd('C:\Users\shires\Documents\MATLAB')
base_angle = 21;

behavior_base_dir = 'Z:\Data\2p\soloData\';
whisker_base_dir = 'Z:\Data\Video\JK\';

mice = {'AH0648','AH0650','AH0651','AH0652','AH0653'};

mouseName = 'AH0648';
sessionName = 'S01';
trial_types = {'rc', 'rf', 'lc', 'lf'};
% trial_types = {'rn', 'ln'};
behavior_d = [behavior_base_dir mouseName '\'];
whisker_d = [whisker_base_dir mouseName sessionName '\'];

if exist([whisker_d, 'touch_hp.mat'],'file')
    error('touch_hp.mat exists.')    
end

load([behavior_d 'behavior.mat']) % loading b of the mouse (all the sessions)

b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
b_session = b{b_ind};

load_fn = [mouseName sessionName '_post.mat'];
load([whisker_d load_fn]); % loading errorlist
% cd('C:\Users\shires\Documents\MATLAB')
%%
if ~isempty(b_ind) % try only the ones with behavior session
    % %%
    filelist=dir([whisker_d '*.measurements']);

    dirTrialNums=zeros(1,size(filelist,1));
    % trialNums=[];  % enter which trial nums to process 

    % %%
    % Assign the trial numbers to existing .measurements files in the directory
    % NOTE : This assumes that the .measurements files have leading numbers
    % corresponding to trial number in string positions 1:end-13 of the file
    % name. These index numbers may need to be changed to match up to the
    % numerical code of the trial number.  (2016/09/05 JK)

    for i=1:length(filelist)
        dirTrialNums(i)=str2double(filelist(i).name(1:end-13)); % extract out the trial number from each measurements file present in directory
    end
    dirTrialNums = setdiff(dirTrialNums,errorlist);
    trialNums = sort(dirTrialNums);
    trialNums = trialNums(~isnan(trialNums));
    trialNums = intersect(trialNums,b{b_ind}.trialNums); % try only the ones with behavior trials

    includef=cell(size(trialNums,1),1);
    for i = 1: length(trialNums)
        includef{i} = num2str(trialNums(i));
    end
end

% %% Make whisker-pole touch space for each type of trial, from 10 randomly selected trials (of each type)
% Currently, only dealing with 4 types of trials: 'rc', 'rf', 'lc', 'lf'
% Should make something different for straight pole touch in S00. 
% 2017/04/11 JK

steps_hp = cell(1,length(trial_types));
num_points_in_hp = cell(1,length(trial_types));
tt_ind = cell(1,length(trial_types));
wl_array = cell(1,length(trial_types));
% touch_points = cell(1,length(trial_types));
touch_hp = cell(1,length(trial_types)); % touch hyperplanes
thp_peak_points = cell(1,length(trial_types)); % touch hyperplane peak points. 2 points for each hyperplane
% %%
% load('wl_array.mat')

for trial_type_num = 1 : length(trial_types)    
% trial_type_num = 1
    tt_ind{trial_type_num} = find(cellfun(@(x) strcmp(x.trialType,trial_types{trial_type_num}),b_session.trials));
    temp_files = cell(length(tt_ind{trial_type_num}),1);
    for j = 1 : length(tt_ind{trial_type_num})
        temp_files{j} = num2str(tt_ind{trial_type_num}(j));
    end
    wl = Whisker.WhiskerTrialLiteArray(whisker_d,'include_files',temp_files);
    wl_array{trial_type_num} = wl;
end
%%
for trial_type_num = 4
    intersect_3d = [];
    wl = wl_array{trial_type_num};
    figure, hold all
    for tnum = 1 : length(wl.trials)
        try        
            top_ind = find(~isnan(wl.trials{tnum}.intersect_coord(:,1)));
            front_ind = find(~isnan(wl.trials{tnum}.intersect_coord(:,2)));
            intersect_ind = intersect(wl.trials{tnum}.pole_available_timepoints,intersect(top_ind,front_ind));
            plot3(wl.trials{tnum}.intersect_coord(intersect_ind,1), wl.trials{tnum}.intersect_coord(intersect_ind,2), ones(1,length(intersect_ind))*wl.trials{tnum}.pole_pos, 'k.', 'MarkerSize', 3)
            intersect_3d = [intersect_3d; wl.trials{tnum}.intersect_coord(intersect_ind,1), wl.trials{tnum}.intersect_coord(intersect_ind,2), ones(length(intersect_ind),1)*wl.trials{tnum}.pole_pos];
        catch
            fprintf('Skipping trial #%d because of index problems \n',tnum);        
        end
    end
    title(wl.trials{1}.trial_type), xlabel('Top-view intersection coord'), ylabel('Front-view intersection coord'), zlabel('Pole position')
end

%% when interested in certain points in the figure
% ttype = 4;
% zvalue = 90050;
% tnum = find(cellfun(@(x) abs(x.pole_pos - zvalue) < 10, wl_array{ttype}.trials))
% wl_array{ttype}.trials{tnum(1)}.trackerFileName
% figure, plot3(wl_array{ttype}.trials{tnum(1)}.intersect_coord(:,1), wl_array{ttype}.trials{tnum(1)}.intersect_coord(:,2), 1:length(wl_array{ttype}.trials{tnum(1)}.intersect_coord(:,1)))
% xlabel('Top-view intersection coord'), ylabel('Front-view intersection coord'), zlabel('Frame #')
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
    j_offset = min(x2d_pix(1,:)) - 1;
    i_offset = min(x2d_pix(2,:)) - 1;
    for j = 1 : length(x2d_pix)
        x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) = x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) + 1;
    end
    stds_pre(i) = std(x2d_proj(find(x2d_proj(:))));
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

bwind = 1;
line1x1 = [];
while(isempty(line1x1))
    line1x1 = find(BW(bwind,:),1,'first'); 
    bwind = bwind + 1;
end
line1y1 = 1; 
line1x2 = line1x1 + tand((max_ind-1)*0.01)*size(x2d_flip,1);
if line1x2 < 1
    line1x2 = 1;
    line1y2 = line1x1 * tand((max_ind-1)*0.01);
elseif line1x2 > size(x2d_flip,2)
    line1x2 = size(x2d_flip,2);
    line1y2 = (size(x2d_flip,2) - line1x1) / tand((max_ind-1)*0.01);
else
    line1y2 = size(x2d_flip,1);
end

bwind = 1;
line2x1 = [];
while(isempty(line2x1))
    line2x1 = find(BW(bwind,:),1,'last'); 
    bwind = bwind + 1;
end
line2y1 = 1;
line2x2 = line2x1 + tand((max_ind-1)*0.01)*size(x2d_flip,1); 
if line2x2 < 1
    line2x2 = 1;
    line2y2 = line2x1 * tand((max_ind-1)*0.01);
elseif line2x2 > size(x2d_flip,2)
    line2x2 = size(x2d_flip,2);
    line2y2 = (size(x2d_flip,2) - line2x1) / tand((max_ind-1)*0.01);
else
    line2y2 = size(x2d_flip,1);
end
figure, imagesc(R), xlabel('Angle = 0:0.01:180'), ylabel('Projected values') 
figure, imagesc(x2d_flip, [0 1000]), hold on, line([line1x1 line1x2],[line1y1 line1y2], 'LineWidth', 3, 'Color', [1 1 1])
line([line2x1 line2x2],[line2y1 line2y2], 'LineWidth', 3, 'Color', [1 1 1])

%% Calculate touch hyperplanes
xmin = -200; xmax = 200; zmin_data = min(intersect_3d(:,3)); zmax_data = max(intersect_3d(:,3)); zdiff = zmax_data - zmin_data; zmax = zmax_data + zdiff; zmin = zmin_data - zdiff;
ymin_data = min(intersect_3d(:,2)); ymax_data = max(intersect_3d(:,2));
z = zmin:zmax;
xyz = zeros((length(z))*(xmax-xmin+1),3);
for i = xmin:xmax
    xyz((i-xmin)*length(z)+1 : (i-xmin+1)*length(z),:) = [ones(length(z),1)*i, zeros(length(z),1), z'];
end

[xyz_psi1, ~, ~] = AxelRot(xyz',psi1,[0 0 1], 0); % rotate psi1 degrees counterclockwise around z axis

zcenter = floor(mean([zmax_data, zmin_data]));
x0 = [0 0 zcenter];
u = [1 tand(psi1) 0];
[xyz_psi2, ~, ~] = AxelRot(xyz_psi1, psi2, u, x0); 
xyz_psi2(:,xyz_psi2(3,:) < zmin_data) = [];
xyz_psi2(:,xyz_psi2(3,:) > zmax_data) = [];
xyz_psi2(:,xyz_psi2(2,:) < ymin_data) = [];
xyz_psi2(:,xyz_psi2(2,:) > ymax_data) = [];

figure, plot3(intersect_3d(:,1),intersect_3d(:,2), intersect_3d(:,3),'k.', 'MarkerSize',3), xlabel('top'), ylabel('front'), zlabel('pos'), hold on
plot3(xyz_psi2(1,:), xyz_psi2(2,:), xyz_psi2(3,:), 'r.', 'MarkerSize',3)

%% ~ 0.5 min (depending on the length of "steps" and the size of xyz_psi2)
%%%%%%%%%%%%%%%%%%%%%% manual selection
steps = 150:200;
%%%%%%%%%%%%%%%%%%%%%% try as short as possible to reduce time next step


% if abs(tand(psi1+90)) > 1
%     transvec = [1/tand(psi1+90) 1 0]; % translation vector, moving xyz_psi2 along its normal vector after projecting to x-y plane.
% else
%     transvec = [1 tand(psi1+90) 0];
% end

intersect_pix = round(intersect_3d);

num_points = zeros(max(steps),1);
parfor i = steps % this is time consuming...
    hp = round(xyz_psi2);
    hp(1,:) = hp(1,:)+i;    
    num_points(i) = sum(ismember(intersect_pix, hp','rows'));
end

figure, plot(steps,num_points(steps), 'k-', 'LineWidth', 3), xlabel('translocation (pix)'), ylabel('# intersection')

steps_hp{trial_type_num} = steps;
num_points_in_hp{trial_type_num} = num_points;
touch_hp{trial_type_num} = xyz_psi2; % Don't round them! (at least at this saving process)
sprintf('trial type #%d processed',trial_type_num)
%% Optional confirmation
% hp_offset = 69;
% figure, plot3(intersect_3d(:,1),intersect_3d(:,2), intersect_3d(:,3),'k.', 'MarkerSize',3), xlabel('top'), ylabel('front'), zlabel('pos'), hold on
% plot3(xyz_psi2(1,:) + hp_offset, xyz_psi2(2,:), xyz_psi2(3,:), 'r.', 'MarkerSize', 3) 

%% Observing peaks
% tt = 1;
% figure, plot(steps_hp{tt}, num_points_in_hp{tt}(steps_hp{tt}))
%%
%% Manual recording of the peaks (must be 2 numbers. If there is only one, type that number twice)
hp_peaks = {[27, 38],[45, 56],[173, 173],[163, 173]};
disp('hp_peaks saved')
%%
save([whisker_d 'touch_hp.mat'],'touch_hp','num_points_in_hp','steps_hp','hp_peaks')
%% Find corners (4 points) of the hyperplane
hp_corners = cell(1,length(touch_hp));
for i = 1 : length(touch_hp)
    zmax = max(touch_hp{i}(3,:));
    zmaxind = find(touch_hp{i}(3,:) == zmax);
    zmaxline = touch_hp{i}(1:2,zmaxind);
    xminind = find(zmaxline(1,:) == min(zmaxline(1,:)));
    zmax_xymin = [zmaxline(:,xminind); ones(1,length(xminind))*zmax];
    xmaxind = find(zmaxline(1,:) == max(zmaxline(1,:)));
    zmax_xymax = [zmaxline(:,xmaxind); ones(1,length(xmaxind))*zmax];
    
    zmin = min(touch_hp{i}(3,:));
    zminind = find(touch_hp{i}(3,:) == zmin);
    zminline = touch_hp{i}(1:2,zminind);
    xminind = find(zminline(1,:) == min(zminline(1,:)));
    zmin_xymin = [zminline(:,xminind); ones(1,length(xminind))*zmin];
    xmaxind = find(zminline(1,:) == max(zminline(1,:)));
    zmin_xymax = [zminline(:,xmaxind); ones(1,length(xmaxind))*zmin];
    
    hp_corners{i} = [zmax_xymin, zmax_xymax, zmin_xymax, zmin_xymin];
end
save('hp_corners_170504.mat','hp_corners','hp_peaks')
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
% trialnum = 70;
% temp_ta = round(wl_array{2}.trials{trialnum}.intersect_coord);
% 
% temp_touch_hp = touch_hp{2};
% temp_touch_hp(1,:) = temp_touch_hp(1,:) + 46;
% temp_touch_hp = sortrows(temp_touch_hp',3)';
% temp_ind = find(temp_touch_hp(3,:) == wl_array{1}.trials{trialnum}.pole_pos);
% temp_touch_hp_ = temp_touch_hp(1:2,temp_ind);
% % figure, plot(temp_touch_hp_(1,:), temp_touch_hp_(2,:),'r.'), hold on, plot(temp_ta(550:end-200,1),temp_ta(550:end-200,2),'k.')
% 
% temp_touch_hp_ = [temp_touch_hp_(1,:), temp_touch_hp_(1,:) + 1, temp_touch_hp_(1,:) - 1, temp_touch_hp_(1,:) - 2; ...
%     temp_touch_hp_(2,:),temp_touch_hp_(2,:),temp_touch_hp_(2,:),temp_touch_hp_(2,:)];
% temp_touch_hp_ = temp_touch_hp_';
% temp_touch_hp_ = unique(temp_touch_hp_,'rows');
% touch_frame = find(ismember(temp_ta,temp_touch_hp_,'rows')==1);
% 
% wl_array{2}.trials{trialnum}.trackerFileName
% 
% figure,  plot(temp_touch_hp_(:,1),temp_touch_hp_(:,2),'r.'), hold on, plot(temp_ta(550:end-200,1),temp_ta(550:end-200,2),'k.')

%%
for i = 1
    wl = wl_array{i};
    trial_type = wl.trials{1}.trial_type;
    eval(['temp_touch_list = touch_list.AH0648S02.',trial_type,';'])
    total_touch_num = 0;
    touch_right_num = 0;
    total_non_touch_num = 0;
    non_touch_wrong_num = 0;
    touch_polygon = [hp_corners{i}+ [hp_peaks{i}(1);0;0], hp_corners{i}+ [hp_peaks{i}(2);0;0]];
    for j = 1 : 10
        touch_frames = [temp_touch_list.touch_protraction{j} + 1; temp_touch_list.touch_retraction{j} + 1];
        non_touch_frames = temp_touch_list.non_touch{j} + 1;        
        
        temp_intersect_coord = [wl.trials{temp_touch_list.trial_num(j)}.intersect_coord,ones(size(wl.trials{temp_touch_list.trial_num(j)}.intersect_coord,1),1)*wl.trials{temp_touch_list.trial_num(j)}.pole_pos];
        touch_intersect_coord = temp_intersect_coord(touch_frames,:);
        non_touch_intersect_coord = temp_intersect_coord(non_touch_frames,:);
        
        total_touch_num = total_touch_num + length(touch_frames);
        touch_right_num = touch_right_num + sum(inhull(touch_intersect_coord,touch_polygon'));
        
        total_non_touch_num = total_non_touch_num + length(non_touch_frames);
        non_touch_wrong_num = non_touch_wrong_num + sum(inhull(non_touch_intersect_coord,touch_polygon'));
    end
end




%%
% pro_touch_right = zeros(9,4); % peak / peak & -1 / peak & -1 & -2 / peak & +1 / peak & +1 & +2 / peak & -1 & +1 / peak & -1 & -2 & +1 / peak & -1 & +1 & +2 / peak & -1 & -2 & +1 & +2
% pro_touch_wrong = zeros(9,4);
% ret_touch_right = zeros(9,4); % peak / peak & +1 / peak & +1 & +2 / peak & -1 / peak & -1 & -2 / peak & +1 & -1 / peak & +1 & +2 & -1 / peak & +1 & -1 & -2 / peak & +1 & +2 & -1 & -2
% ret_touch_wrong = zeros(9,4);
% non_touch_right = zeros(9,4);
% non_touch_wrong = zeros(9,4);
% for i = 1 : 4
%     wl = wl_array{i};
%     trial_type = wl.trials{1}.trial_type;
%     eval(['temp_touch_list = touch_list.AH0648S02.',trial_type,';'])
%     temp_th = touch_hp{i}; % temp_touch_hyperplane
%     temp_trial_ind = temp_touch_list.trial_num;
%     % checking tracker file names
%     for j = 1 : 10
%         if wl.trialNums(temp_trial_ind(j)) ~= str2double(temp_touch_list.tracker_filename(j));
%             error('Tracker file name mismatch in trial index #%d of trial type #%d', j, i);
%         end
%     end
%     for j = 1 : 10
%         temp_pro_frames = temp_touch_list.touch_protraction{j} +1; % make frames start from 1, not 0
%         temp_ret_frames = temp_touch_list.touch_retraction{j} +1; % make frames start from 1, not 0
%         temp_non_frames = temp_touch_list.non_touch{j} +1; % make frames start from 1, not 0
%         temp_trial = wl.trials(temp_trial_ind(j));
%         v = VideoReader([temp_trial.trackerFileName,'.mp4']);
%         if v.NumberOfFrames ~= length(temp_trial.intersect_coord)
%             error(['# of video frames does not match that of whisker-tracker in trial #',temp_trial.trackerFileName])
%         end
%         for k = 1 : 9
%             temp_th_pp_ind = find(temp_th(3,:) == temp_trial.pole_pos); % temp_touch_hyperplane_pole_position_index
%             temp_th_pp = temp_th(1:2,temp_th_pp_ind); % temp_touch_hyperplane_pole_position
%             temp_th_pro_peak = temp_th_pp + [ones(1,length(temp_th_pp)) * hp_peaks{i}(2); zeros(1,length(temp_th_pp))];
%             temp_th_ret_peak = temp_th_pp + [ones(1,length(temp_th_pp)) * hp_peaks{i}(1); zeros(1,length(temp_th_pp))];      
%             temp_intersection = round(temp_trial.intersect_coord);
%             switch k
%                 case 1
%                     temp_th_pro = temp_th_pro_peak';
%                     temp_th_ret = temp_th_ret_peak';                    
%                 case 2
%                     temp_th_pro = [temp_th_pro_peak(1,:), temp_th_pro_peak(1,:) - 1; ...
%                         temp_th_pro_peak(2,:), temp_th_pro_peak(2,:)]';
%                     temp_th_ret = [temp_th_ret_peak(1,:), temp_th_ret_peak(1,:) + 1; ...
%                         temp_th_ret_peak(2,:), temp_th_ret_peak(2,:)]';
%                 case 3
%                     temp_th_pro = [temp_th_pro_peak(1,:), temp_th_pro_peak(1,:) - 1, temp_th_pro_peak(1,:) - 2; ...
%                         temp_th_pro_peak(2,:), temp_th_pro_peak(2,:), temp_th_pro_peak(2,:)]';
%                     temp_th_ret = [temp_th_ret_peak(1,:), temp_th_ret_peak(1,:) + 1, temp_th_ret_peak(1,:) + 2; ...
%                         temp_th_ret_peak(2,:), temp_th_ret_peak(2,:), temp_th_ret_peak(2,:)]';
%                 case 4
%                     temp_th_pro = [temp_th_pro_peak(1,:), temp_th_pro_peak(1,:) + 1; ...
%                         temp_th_pro_peak(2,:), temp_th_pro_peak(2,:)]';
%                     temp_th_ret = [temp_th_ret_peak(1,:), temp_th_ret_peak(1,:) - 1; ...
%                         temp_th_ret_peak(2,:), temp_th_ret_peak(2,:)]';
%                 case 5
%                     temp_th_pro = [temp_th_pro_peak(1,:), temp_th_pro_peak(1,:) + 1, temp_th_pro_peak(1,:) + 2; ...
%                         temp_th_pro_peak(2,:), temp_th_pro_peak(2,:), temp_th_pro_peak(2,:)]';
%                     temp_th_ret = [temp_th_ret_peak(1,:), temp_th_ret_peak(1,:) - 1, temp_th_ret_peak(1,:) - 2; ...
%                         temp_th_ret_peak(2,:), temp_th_ret_peak(2,:), temp_th_ret_peak(2,:)]';
%                 case 6
%                     temp_th_pro = [temp_th_pro_peak(1,:), temp_th_pro_peak(1,:) - 1, temp_th_pro_peak(1,:) + 1; ...
%                         temp_th_pro_peak(2,:), temp_th_pro_peak(2,:), temp_th_pro_peak(2,:)]';
%                     temp_th_ret = [temp_th_ret_peak(1,:), temp_th_ret_peak(1,:) + 1, temp_th_ret_peak(1,:) - 1; ...
%                         temp_th_ret_peak(2,:), temp_th_ret_peak(2,:), temp_th_ret_peak(2,:)]';
%                 case 7
%                      temp_th_pro = [temp_th_pro_peak(1,:), temp_th_pro_peak(1,:) - 1, temp_th_pro_peak(1,:) - 2, temp_th_pro_peak(1,:) + 1; ...
%                         temp_th_pro_peak(2,:), temp_th_pro_peak(2,:), temp_th_pro_peak(2,:), temp_th_pro_peak(2,:)]';
%                     temp_th_ret = [temp_th_ret_peak(1,:), temp_th_ret_peak(1,:) + 1, temp_th_ret_peak(1,:) + 2, temp_th_ret_peak(1,:) - 1; ...
%                         temp_th_ret_peak(2,:), temp_th_ret_peak(2,:), temp_th_ret_peak(2,:), temp_th_ret_peak(2,:)]';        
%                 case 8
%                     temp_th_pro = [temp_th_pro_peak(1,:), temp_th_pro_peak(1,:) - 1, temp_th_pro_peak(1,:) + 1, temp_th_pro_peak(1,:) + 2; ...
%                         temp_th_pro_peak(2,:), temp_th_pro_peak(2,:), temp_th_pro_peak(2,:), temp_th_pro_peak(2,:)]';
%                     temp_th_ret = [temp_th_ret_peak(1,:), temp_th_ret_peak(1,:) + 1, temp_th_ret_peak(1,:) - 1, temp_th_ret_peak(1,:) - 2; ...
%                         temp_th_ret_peak(2,:), temp_th_ret_peak(2,:), temp_th_ret_peak(2,:), temp_th_ret_peak(2,:)]';        
%                 case 9
%                     temp_th_pro = [temp_th_pro_peak(1,:), temp_th_pro_peak(1,:) - 1, temp_th_pro_peak(1,:) - 2, temp_th_pro_peak(1,:) + 1, temp_th_pro_peak(1,:) + 2; ...
%                         temp_th_pro_peak(2,:), temp_th_pro_peak(2,:), temp_th_pro_peak(2,:), temp_th_pro_peak(2,:), temp_th_pro_peak(2,:)]';
%                     temp_th_ret = [temp_th_ret_peak(1,:), temp_th_ret_peak(1,:) + 1, temp_th_ret_peak(1,:) + 2, temp_th_ret_peak(1,:) - 1, temp_th_ret_peak(1,:) - 2; ...
%                         temp_th_ret_peak(2,:), temp_th_ret_peak(2,:), temp_th_ret_peak(2,:), temp_th_ret_peak(2,:), temp_th_ret_peak(2,:)]';        
%             end
%             temp_th_pro = unique(temp_th_pro,'rows');
%             temp_th_ret = unique(temp_th_ret,'rows');
%             
%             temp_th_pro_frames = find(ismember(temp_intersection,temp_th_pro,'rows') == 1);
%             temp_th_ret_frames = find(ismember(temp_intersection,temp_th_ret,'rows') == 1);            
%                        
%             pro_touch_right(k,i) = pro_touch_right(k,i) + sum(ismember(temp_pro_frames,temp_th_pro_frames));
%             pro_touch_wrong(k,i) = pro_touch_wrong(k,i) + length(temp_pro_frames) - sum(ismember(temp_pro_frames,temp_th_pro_frames));
%             ret_touch_right(k,i) = ret_touch_right(k,i) + sum(ismember(temp_ret_frames,temp_th_ret_frames));
%             ret_touch_wrong(k,i) = ret_touch_wrong(k,i) + length(temp_ret_frames) - sum(ismember(temp_ret_frames,temp_th_ret_frames));
%             non_touch_right(k,i) = non_touch_right(k,i) + length(temp_non_frames) - sum(ismember(temp_non_frames,temp_th_pro_frames)) - sum(ismember(temp_non_frames,temp_th_ret_frames));
%             non_touch_wrong(k,i) = non_touch_wrong(k,i) + sum(ismember(temp_non_frames,temp_th_pro_frames)) + sum(ismember(temp_non_frames,temp_th_ret_frames));
%         end
%     end
% end
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

%%
tid = [0 1]; % Set trajectory ID to view
wl = Whisker.WhiskerTrialLiteArray_2pad(whisker_d);
Whisker.viewdouble_WhiskerTrialLiteArray(wl,tid)
