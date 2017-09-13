close all
% clear
base_angle = 21;
behavior_base_dir = 'Z:\Data\2p\soloData\';
whisker_base_dir = 'Z:\Data\Video\JK\';
trial_types = {'rc', 'rf', 'lc', 'lf'};
mice = {'AH0648','AH0650','AH0651','AH0652','AH0653'};
% sessionNum = {[1:4,6:15],[1,2,4:6,8:10],[2,4:6,8:19],[6,8:13],[2,4:8,13:19]};

sessions{1} = {[4,6,7], [13:15]}; % 0648
sessions{2} = {[4:6], [10]}; % 0650
sessions{3} = {[4:6], [6,8,10]}; % 0651
sessions{4} = {[6,8], [10,11,13]}; % 0652
sessions{5} = {[4:6], [16:18]}; % 0653

mousensession = '323'; % defining which mouse, which session. (1) mouse, (2) before/after (3) session number within before/after
try
    mouseName = mice{str2double(mousensession(1))};
    sessionName = ['S', num2str(sessions{str2double(mousensession(1))}{str2double(mousensession(2))}(str2double(mousensession(3))))];
catch
    error('No matched mouse or session')
end

whisker_d = [whisker_base_dir mouseName sessionName '\'];
if ~exist('b','var') || ~iscell(b) || ~isfield(b{1},'mouseName') || ~strcmp(b{1}.mouseName,mouseName)
    load([behavior_base_dir mouseName filesep 'behavior.mat']) % load b
end

b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
b_session = b{b_ind};
tt_ind = cell(1,length(trial_types));
wl = cell(1,length(trial_types));

for ti = 1 : length(trial_types)    
    tt_ind{ti} = find(cellfun(@(x) strcmp(x.trialType,trial_types{ti}),b_session.trials));
    temp_files = cell(length(tt_ind{ti}),1);
    for j = 1 : length(tt_ind{ti})
        temp_files{j} = num2str(tt_ind{ti}(j));
    end
    wl{ti} = Whisker.WhiskerTrialLiteArray_2pad(whisker_d,'include_files',temp_files);     
end
% assign deltakappa and theta. deltakappa is a true deltakappa, calculated by subtracting the kappa value of one frame before the first touch frame in each chunk
%%
deltakappa_t = cell(1,length(trial_types));
deltakappa_f = cell(1,length(trial_types));
theta_t = cell(1,length(trial_types));
theta_f = cell(1,length(trial_types));
deltatheta_t = cell(1,length(trial_types));
deltatheta_f = cell(1,length(trial_types));
theta_whisking_t_prot = cell(1,length(trial_types)); % Don't need to divide by trial types, but just for the sake of computation...
theta_whisking_t_ret = cell(1,length(trial_types)); % Don't need to divide by trial types, but just for the sake of computation...
theta_whisking_f_prot = cell(1,length(trial_types)); % Don't need to divide by trial types, but just for the sake of computation...
theta_whisking_f_ret = cell(1,length(trial_types)); % Don't need to divide by trial types, but just for the sake of computation...
kappa_whisking_t_prot = cell(1,length(trial_types)); % Don't need to divide by trial types, but just for the sake of computation...
kappa_whisking_t_ret = cell(1,length(trial_types)); % Don't need to divide by trial types, but just for the sake of computation...
kappa_whisking_f_prot = cell(1,length(trial_types)); % Don't need to divide by trial types, but just for the sake of computation...
kappa_whisking_f_ret = cell(1,length(trial_types)); % Don't need to divide by trial types, but just for the sake of computation...


for ti = 1 : length(trial_types)
% for ti = 1     
    for j = 1 : length(wl{ti}.trials)
%     for j = 1 
        wt = wl{ti}.trials{j};
        trial_temp_ind = find(cellfun(@(x) x.trialNum == wt.trialNum,b_session.trials));
        trial_temp = b_session.trials{trial_temp_ind};
        
        if isempty(find(diff(wt.time{1}) > 0.004,1)) && length(wt.time{1}) == wt.videoFrames % Only calculate from trials with NO lost frame
            temp_touch = wt.th_touch_frames;
            temp_whisking = setdiff(wt.pole_available_timepoints,wt.th_touch_frames);
            if ~isempty(temp_touch)
                touch_diff_inds = [0;find(diff(temp_touch) - 1);length(temp_touch)];
                for cind = 1 : length(touch_diff_inds)-1
                    temp_ind = temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1));
                    if length(wt.thetaAtBase{1}(temp_ind)) * length(wt.deltaKappa{1}(temp_ind)) * length(wt.thetaAtBase{2}(temp_ind)) * length(wt.deltaKappa{2}(temp_ind)) > 0 % non of these are empty
                        deltakappa_t{ti}{end+1} = wt.deltaKappa{1}(temp_ind) - wt.deltaKappa{1}(temp_ind(1)-1);
                        deltakappa_f{ti}{end+1} = wt.deltaKappa{2}(temp_ind) - wt.deltaKappa{2}(temp_ind(1)-1);
                        theta_t{ti}{end+1} = wt.thetaAtBase{1}(temp_ind);
                        theta_f{ti}{end+1} = wt.thetaAtBase{2}(temp_ind);
                        deltatheta_t{ti}{end+1} = wt.thetaAtBase{1}(temp_ind) - wt.thetaAtBase{1}(temp_ind(1));
                        deltatheta_f{ti}{end+1} = wt.thetaAtBase{2}(temp_ind) - wt.thetaAtBase{2}(temp_ind(1));
                    end
                end
            end
            if ~isempty(temp_whisking)
                whisking_diff_inds = [0, find(diff(temp_whisking) - 1), length(temp_whisking)];
                for cind = 1 : length(whisking_diff_inds) - 1
                    temp_whisking_bout = temp_whisking(whisking_diff_inds(cind)+1:whisking_diff_inds(cind+1)); % in frames
                    if length(temp_whisking_bout) > 2
                        temp_whisking_diff = [0, diff(wt.thetaAtBase{1}(temp_whisking_bout))];
                        
                        temp_prot = find(temp_whisking_diff > 0); % indices of temp_whisking_bout
                        if length(temp_prot) > 2
                            temp_prot_diff = [0, find(diff(temp_prot)-1), length(temp_prot)]; 
                            for pind = 1 : length(temp_prot_diff) - 1
                                temp_p_ind = temp_whisking_bout(temp_prot(temp_prot_diff(pind)+1):temp_prot(temp_prot_diff(pind+1)));
                                if length(wt.thetaAtBase{1}(temp_p_ind)) * length(wt.deltaKappa{1}(temp_p_ind)) * length(wt.thetaAtBase{2}(temp_p_ind)) * length(wt.deltaKappa{2}(temp_p_ind)) > 0 % non of these are empty
                                    theta_whisking_t_prot{ti}{end+1} = wt.thetaAtBase{1}(temp_p_ind);
                                    theta_whisking_f_prot{ti}{end+1} = wt.thetaAtBase{2}(temp_p_ind);
                                    kappa_whisking_t_prot{ti}{end+1} = wt.deltaKappa{1}(temp_p_ind);
                                    kappa_whisking_f_prot{ti}{end+1} = wt.deltaKappa{2}(temp_p_ind);
                                end
                            end                        
                        end
                        
                        temp_ret = find(temp_whisking_diff < 0); % indices of temp_whisking_bout
                        if length(temp_ret) > 2
                            temp_ret_diff = [0, find(diff(temp_ret)-1), length(temp_ret)]; 
                            for rind = 1 : length(temp_ret_diff) - 1
                                temp_r_ind = temp_whisking_bout(temp_ret(temp_ret_diff(rind)+1):temp_ret(temp_ret_diff(rind+1)));
                                if length(wt.thetaAtBase{1}(temp_r_ind)) * length(wt.deltaKappa{1}(temp_r_ind)) * length(wt.thetaAtBase{2}(temp_r_ind)) * length(wt.deltaKappa{2}(temp_r_ind)) > 0 % non of these are empty
                                    theta_whisking_t_ret{ti}{end+1} = wt.thetaAtBase{1}(temp_r_ind);
                                    theta_whisking_f_ret{ti}{end+1} = wt.thetaAtBase{2}(temp_r_ind);
                                    kappa_whisking_t_ret{ti}{end+1} = wt.deltaKappa{1}(temp_r_ind);
                                    kappa_whisking_f_ret{ti}{end+1} = wt.deltaKappa{2}(temp_r_ind);
                                end
                            end       
                        end
                    end
                end
            end
        end              
    end
end

%%

figure, hold all
% for i = 1 : length(theta_whisking_t_prot)
for i = 1
    for j = 1 : length(theta_whisking_t_prot{i})
        try
            plot3(theta_whisking_t_prot{i}{j}+base_angle,theta_whisking_f_prot{i}{j},kappa_whisking_t_prot{i}{j},'r.');        
        catch
        end
    end
end

% for i = 1 : length(theta_whisking_t_ret)
for i = 1
    for j = 1 : length(theta_whisking_t_ret{i})    
        try
            plot3(theta_whisking_t_ret{i}{j}+base_angle,theta_whisking_f_ret{i}{j},kappa_whisking_t_ret{i}{j},'b.');
        catch
        end
    end
end

title({[mouseName sessionName]; 'Free whisking theta'}), xlabel('Top-view theta'), ylabel('Front-view theta'), zlabel('Top-view kappa')
%%
data_prot = [];
data_ret = [];
for i = 1 : length(theta_whisking_t_prot)
% for i = 1
    for j = 1 : length(theta_whisking_t_prot{i})
        try
            data_prot = [data_prot; theta_whisking_t_prot{i}{j}'+base_angle,theta_whisking_f_prot{i}{j}',kappa_whisking_t_prot{i}{j}'];        
        catch
        end
    end
end

for i = 1 : length(theta_whisking_t_ret)
% for i = 1
    for j = 1 : length(theta_whisking_t_ret{i})    
        try
            data_ret = [data_ret; theta_whisking_t_ret{i}{j}'+base_angle,theta_whisking_f_ret{i}{j}',kappa_whisking_t_ret{i}{j}'];
        catch
        end
    end
end


A = viewmtx(90,0);
x4d_prot = [data_prot, ones(size(data_prot,1),1)]';
x2d_prot = A*x4d_prot;
x4d_ret = [data_ret, ones(size(data_ret,1),1)]';
x2d_ret = A*x4d_ret;

h = figure; plot(x2d_prot(1,:), x2d_prot(2,:),'r.', 'MarkerSize',1); hold on; plot(x2d_ret(1,:), x2d_ret(2,:),'b.', 'MarkerSize',1);

title({[mouseName sessionName]; 'Free whisking theta'}), xlabel('Top-view theta'), ylabel('Front-view theta'), zlabel('Top-view kappa');
h.XLim = [-21 -19];