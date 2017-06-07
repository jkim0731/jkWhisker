close all

behavior_base_dir = 'Z:\Data\2p\soloData\';
whisker_base_dir = 'Z:\Data\Video\JK\';

mice = {'AH0648','AH0650','AH0651','AH0652','AH0653'};

mouseName = 'AH0650';
sessionName = 'S21';
% trial_types = {'rc', 'rf', 'lc', 'lf'};
trial_types = {'rn', 'ln'};
behavior_d = [behavior_base_dir mouseName '\'];
whisker_d = [whisker_base_dir mouseName sessionName '\'];

load([behavior_d 'behavior.mat']) % loading b of the mouse (all the sessions)

b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
b_session = b{b_ind};

load_fn = [mouseName sessionName '_post.mat'];
load([whisker_d load_fn]); % loading errorlist

if ~isempty(b_ind) % try only the ones with behavior session
    filelist=dir([whisker_d '*.measurements']);
    dirTrialNums=zeros(1,size(filelist,1));

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

% %%
% tid = [0 1]; % Set trajectory ID to view
% wl = Whisker.WhiskerTrialLiteArray_2pad(whisker_d);
% Whisker.viewdouble_WhiskerTrialLiteArray(wl,tid)
% %%
% ttf_ind = [];
% for i = 1 : length(wl.trials)
%     if ~isempty(wl.trials{i}.th_touch_frames)
%         ttf_ind = [ttf_ind; i];
%     end
% end
% 
% %%
% ttf_fn = cell(length(ttf_ind),1);
% for i = 1 : length(ttf_ind)
%     ttf_fn{i} = wl.trials{ttf_ind(i)}.trackerFileName;
% end
% 
% %%
% t_ind = 170;
% confirmed_touch_frame = 551:556;
% figure, 
% % plot(wl.trials{t_ind}.intersect_coord(wl.trials{t_ind}.pole_available_timepoints,1), wl.trials{t_ind}.intersect_coord(wl.trials{t_ind}.pole_available_timepoints,2), 'b.'), hold on,
% plot(wl.trials{t_ind}.th_polygon(:,1), wl.trials{t_ind}.th_polygon(:,2), 'k.'), hold on,
% % plot(wl.trials{t_ind}.intersect_coord(confirmed_touch_frame,1), wl.trials{t_ind}.intersect_coord(confirmed_touch_frame,2), 'b.')
% plot(wl.trials{t_ind}.intersect_coord(:,1), wl.trials{t_ind}.intersect_coord(:,2), 'b.')

%%
tt_ind = cell(1,length(trial_types));
wl_array = cell(1,length(trial_types));

for trial_type_num = 1 : length(trial_types)    
    tt_ind{trial_type_num} = find(cellfun(@(x) strcmp(x.trialType,trial_types{trial_type_num}),b_session.trials));
    temp_files = cell(length(tt_ind{trial_type_num}),1);
    for j = 1 : length(tt_ind{trial_type_num})
        temp_files{j} = num2str(tt_ind{trial_type_num}(j));
    end
    wl = Whisker.WhiskerTrialLiteArray_2pad(whisker_d,'include_files',temp_files);
    wl_array{trial_type_num} = wl;
end
%%
base_angle = 21;
fw1 = []; fw2 = []; fw3 = []; fw4 = []; fw5 = []; fw6 = [];% free-whisking
for wl_array_ind = 1 : length(trial_types)
    figure, hold all
    sp1 = []; sp2 = []; sp3 = []; sp4 = []; sp5 = []; sp6 = []; %subplot 1~6
    for wl_ind = 1 : length(wl_array{wl_array_ind}.trials)
        wt = wl_array{wl_array_ind}.trials{wl_ind};      
        fw_frames = find(wt.time{1} < 0.0032 * 500); % first 500 frames
        sp1 = [sp1; wt.thetaAtBase{1}(wt.th_touch_frames)'+base_angle, wt.thetaAtBase{2}(wt.th_touch_frames)'];
        sp2 = [sp2; wt.deltaKappa{1}(wt.th_touch_frames)', wt.deltaKappa{2}(wt.th_touch_frames)'];
        sp3 = [sp3; wt.thetaAtBase{1}(wt.th_touch_frames)'+base_angle, wt.deltaKappa{1}(wt.th_touch_frames)'];
        sp4 = [sp4; wt.thetaAtBase{2}(wt.th_touch_frames)', wt.deltaKappa{2}(wt.th_touch_frames)'];
        sp5 = [sp5; wt.thetaAtBase{1}(wt.th_touch_frames)'+base_angle, wt.deltaKappa{2}(wt.th_touch_frames)'];
        sp6 = [sp6; wt.thetaAtBase{1}(wt.th_touch_frames)'+base_angle, wt.intersect_coord(wt.th_touch_frames,1)];
        fw1 = [fw1; wt.thetaAtBase{1}(fw_frames)'+base_angle, wt.thetaAtBase{2}(fw_frames)'];
        fw2 = [fw2; wt.deltaKappa{1}(fw_frames)', wt.deltaKappa{2}(fw_frames)'];
        fw3 = [fw3; wt.thetaAtBase{1}(fw_frames)'+base_angle, wt.deltaKappa{1}(fw_frames)'];
        fw4 = [fw4; wt.thetaAtBase{2}(fw_frames)', wt.deltaKappa{2}(fw_frames)'];
        fw5 = [fw5; wt.thetaAtBase{1}(fw_frames)'+base_angle, wt.deltaKappa{2}(fw_frames)'];        
    end
    subplot(3,2,1), % top-view theta vs front-view theta
    plot(sp1(:,1),sp1(:,2),'k.'), title('top-view theta vs front-view theta'), xlim([-30 30]), ylim([-10 30]);
    subplot(3,2,2), % top-view kappa vs front-view kappa
    plot(sp2(:,1),sp2(:,2),'k.'), title('top-view kappa vs front-view kappa'), xlim([-0.03 0.02]), ylim([-0.02 0.02]);
    subplot(3,2,3) % top-view theta vs top-view kappa
    plot(sp3(:,1),sp3(:,2),'k.'), title('top-view theta vs kappa'), xlim([-30 30]), ylim([-0.05 0.05]);
    subplot(3,2,4) % front-view theta vs front-view kappa
    plot(sp4(:,1),sp4(:,2),'k.'), title('front-view theta vs kappa'), xlim([-10 30]), ylim([-0.02 0.02]);
    subplot(3,2,5) % top-view theta vs front-view kappa
    plot(sp5(:,1),sp5(:,2),'k.'), title('top-view theta vs front-view kappa'), xlim([-30 30]), ylim([-0.02 0.02]);
    subplot(3,2,6) % top-view theta vs top-view intersection coordinate
    plot(sp6(:,1),sp6(:,2),'k.'), title('top-view theta vs top-view ISC'), xlim([-30 30]), ylim([100 400]);    
end

% free whisking
figure, 
subplot(3,2,1), % top-view theta vs front-view theta
    plot(fw1(:,1),fw1(:,2),'k.'), title(sprintf('top-view theta vs front-view theta')), xlim([-30 30]), ylim([-10 30]);
    subplot(3,2,2), % top-view kappa vs front-view kappa
    plot(fw2(:,1),fw2(:,2),'k.'), title(sprintf('top-view kappa vs front-view kappa')), xlim([-0.03 0.02]), ylim([-0.02 0.02]);
    subplot(3,2,3) % top-view theta vs top-view kappa
    plot(fw3(:,1),fw3(:,2),'k.'), title(sprintf('top-view theta vs kappa')), xlim([-30 30]), ylim([-0.05 0.05]);
    subplot(3,2,4) % front-view theta vs front-view kappa
    plot(fw4(:,1),fw4(:,2),'k.'), title(sprintf('front-view theta vs kappa')), xlim([-10 30]), ylim([-0.02 0.02]);
    subplot(3,2,5) % top-view theta vs front-view kappa
    plot(fw5(:,1),fw5(:,2),'k.'), title(sprintf('top-view theta vs front-view kappa')), xlim([-30 30]), ylim([-0.02 0.02]);



        
        
        