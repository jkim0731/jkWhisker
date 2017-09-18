

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
    'barRadius',15.3,'faceSideInImage', 'bottom', 'framePeriodInSec',.0032,...
    'imagePixelDimsXY',[width height],'pxPerMm',26.23,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward')

Whisker.makeAllDirectory_WhiskerSignalTrial(whisker_d,'include_files',includef,'polyRoiInPix',[20 80]);
Whisker.makeAllDirectory_WhiskerTrialLiteI(whisker_d,'include_files',includef,'r_in_mm',3,'calc_forces',true,'whisker_radius_at_base', 36.5,'whisker_length', 18,'baseline_time_or_kappa_value',0);
wl = Whisker.WhiskerTrialLiteArray(whisker_d);
save([d mouseName sessionName '-WTLIA.mat'],'wl');

%%
tid = [0 1]; % Set trajectory ID to view
wl = Whisker.WhiskerTrialLiteArray_2pad(whisker_d);
Whisker.viewdouble_WhiskerTrialLiteArray(wl,tid)
%%
ttf_ind = [];
for i = 1 : length(wl.trials)
    if ~isempty(wl.trials{i}.th_touch_frames)
        ttf_ind = [ttf_ind; i];
    end
end

%%
ttf_fn = cell(length(ttf_ind),1);
for i = 1 : length(ttf_ind)
    ttf_fn{i} = wl.trials{ttf_ind(i)}.trackerFileName;
end
