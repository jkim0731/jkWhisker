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
                        deltakappa_t{ti}{end+1} = wt.deltaKappa{1}(temp_ind) - wt.deltaKappa{1}(temp_ind(1));
                        deltakappa_f{ti}{end+1} = wt.deltaKappa{2}(temp_ind) - wt.deltaKappa{2}(temp_ind(1));
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

%% collapse and average only those from protraction touch
% %% assign deltakappa and theta. deltakappa is a true deltakappa, calculated by subtracting the kappa value of one frame before the first touch frame in each chunk
protraction_deltakappa_t = cell(1,length(trial_types));
protraction_deltakappa_f = cell(1,length(trial_types));
protraction_theta_t = cell(1,length(trial_types));
protraction_theta_f = cell(1,length(trial_types));
protraction_deltatheta_t = cell(1,length(trial_types));
protraction_deltatheta_f = cell(1,length(trial_types));
for ti = 1 : length(trial_types)
    for j = 1 : length(wl{ti}.trials) 
        wt = wl{ti}.trials{j};
        trial_temp_ind = find(cellfun(@(x) x.trialNum == wt.trialNum,b_session.trials));
        trial_temp = b_session.trials{trial_temp_ind};

        temp_touch = wt.th_touch_frames;
        if ~isempty(temp_touch)
            touch_diff_inds = [0;find(diff(temp_touch) - 1);length(temp_touch)];
            for cind = 1 : length(touch_diff_inds)-1
                if length(wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * ...
                        length(wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) ...
                        > 0 ... % non of these are empty                    
                    && wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1)+1) >= wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1)) % and theta_top is increasing when the touch started
                    protraction_deltakappa_t{ti}{end+1} = wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1));
                    protraction_deltakappa_f{ti}{end+1} = wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1));
                    protraction_theta_t{ti}{end+1} = wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                    protraction_theta_f{ti}{end+1} = wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                    protraction_deltatheta_t{ti}{end+1} = wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1));
                    protraction_deltatheta_f{ti}{end+1} = wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1));
                end
            end
        end
    end
end

% %% collapsing deltakappa and theta values into one matrix for easier averaging

protraction_dkt = cell(1,length(trial_types)); %delta kappa _t
protraction_dkf = cell(1,length(trial_types)); % delta kappa _f
protraction_tht = cell(1,length(trial_types)); % theta _t
protraction_thf = cell(1,length(trial_types)); % theta _f
protraction_dtt = cell(1,length(trial_types)); % delta theta _t
protraction_dtf = cell(1,length(trial_types)); % delta theta _f

for ti = 1 : length(trial_types)
    protraction_dkt{ti} = zeros(length(protraction_deltakappa_t{ti}),max(cellfun(@(x) length(x),protraction_deltakappa_t{ti})));
    protraction_dkf{ti} = zeros(size(protraction_dkt{ti}));
    protraction_tht{ti} = zeros(size(protraction_dkt{ti}));
    protraction_thf{ti} = zeros(size(protraction_dkt{ti}));
    protraction_dtt{ti} = zeros(size(protraction_dkt{ti}));
    protraction_dtf{ti} = zeros(size(protraction_dkt{ti}));
end
for tind = 1 : length(trial_types)
    for i = 1 : length(protraction_deltakappa_t{tind})
        protraction_dkt{tind}(i,1:length(protraction_deltakappa_t{tind}{i})) = protraction_deltakappa_t{tind}{i}; if length(protraction_deltakappa_t{tind}{i}) < size(protraction_dkt,2), protraction_dkt{tind}(i,length(protraction_deltakappa_t{tind}{i})+1:end) = NaN; end;
        protraction_dkf{tind}(i,1:length(protraction_deltakappa_f{tind}{i})) = protraction_deltakappa_f{tind}{i}; if length(protraction_deltakappa_f{tind}{i}) < size(protraction_dkt,2), protraction_dkf{tind}(i,length(protraction_deltakappa_f{tind}{i})+1:end) = NaN; end;
        protraction_tht{tind}(i,1:length(protraction_theta_t{tind}{i})) = protraction_theta_t{tind}{i}; if length(protraction_theta_t{tind}{i}) < size(protraction_dkt,2), protraction_tht{tind}(i,length(protraction_theta_t{tind}{i})+1:end) = NaN; end;
        protraction_thf{tind}(i,1:length(protraction_theta_f{tind}{i})) = protraction_theta_f{tind}{i}; if length(protraction_theta_f{tind}{i}) < size(protraction_dkt,2), protraction_thf{tind}(i,length(protraction_theta_f{tind}{i})+1:end) = NaN; end;
        protraction_dtt{tind}(i,1:length(protraction_deltatheta_t{tind}{i})) = protraction_deltatheta_t{tind}{i}; if length(protraction_deltatheta_t{tind}{i}) < size(protraction_dkt,2), protraction_dtt{tind}(i,length(protraction_deltatheta_t{tind}{i})+1:end) = NaN; end;
        protraction_dtf{tind}(i,1:length(protraction_deltatheta_f{tind}{i})) = protraction_deltatheta_f{tind}{i}; if length(protraction_deltatheta_f{tind}{i}) < size(protraction_dkt,2), protraction_dtf{tind}(i,length(protraction_deltatheta_f{tind}{i})+1:end) = NaN; end;        
    end
end

%% averaging deltakappa and theta values altogether 
xlim_max = 15;
figure, hold all
subplot(221), hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(1:xlim_max,nanmean(protraction_dkt{i}(:,1:xlim_max)),nanstd(protraction_dkt{i}(:,1:xlim_max)), 'Color', 'b', 'LineWidth', 3) 
        case 2
            errorbar(1:xlim_max,nanmean(protraction_dkt{i}(:,1:xlim_max)),nanstd(protraction_dkt{i}(:,1:xlim_max)), 'Color', 'c', 'LineWidth', 3) 
        case 3
            errorbar(1:xlim_max,nanmean(protraction_dkt{i}(:,1:xlim_max)),nanstd(protraction_dkt{i}(:,1:xlim_max)), 'Color', 'r', 'LineWidth', 3) 
        case 4
            errorbar(1:xlim_max,nanmean(protraction_dkt{i}(:,1:xlim_max)),nanstd(protraction_dkt{i}(:,1:xlim_max)), 'Color', 'm', 'LineWidth', 3) 
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Frames (/3.2 ms)'), ylabel('\Delta\kappa (/mm)'), title('\Delta\kappa_T_o_p', 'FontSize', 25)

subplot(222), hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(1:xlim_max,nanmean(protraction_dkf{i}(:,1:xlim_max)),nanstd(protraction_dkf{i}(:,1:xlim_max)), 'Color', 'b', 'LineWidth', 3) 
        case 2
            errorbar(1:xlim_max,nanmean(protraction_dkf{i}(:,1:xlim_max)),nanstd(protraction_dkf{i}(:,1:xlim_max)), 'Color', 'c', 'LineWidth', 3) 
        case 3
            errorbar(1:xlim_max,nanmean(protraction_dkf{i}(:,1:xlim_max)),nanstd(protraction_dkf{i}(:,1:xlim_max)), 'Color', 'r', 'LineWidth', 3) 
        case 4
            errorbar(1:xlim_max,nanmean(protraction_dkf{i}(:,1:xlim_max)),nanstd(protraction_dkf{i}(:,1:xlim_max)), 'Color', 'm', 'LineWidth', 3) 
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Frames (/3.2 ms)'), ylabel('\Delta\kappa (/mm)'), title('\Delta\kappa_F_r_o_n_t', 'FontSize', 25)

subplot(223), hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(1:xlim_max,nanmean(protraction_dtt{i}(:,1:xlim_max)),nanstd(protraction_dtt{i}(:,1:xlim_max)), 'Color', 'b', 'LineWidth', 3) 
        case 2
            errorbar(1:xlim_max,nanmean(protraction_dtt{i}(:,1:xlim_max)),nanstd(protraction_dtt{i}(:,1:xlim_max)), 'Color', 'c', 'LineWidth', 3) 
        case 3
            errorbar(1:xlim_max,nanmean(protraction_dtt{i}(:,1:xlim_max)),nanstd(protraction_dtt{i}(:,1:xlim_max)), 'Color', 'r', 'LineWidth', 3) 
        case 4
            errorbar(1:xlim_max,nanmean(protraction_dtt{i}(:,1:xlim_max)),nanstd(protraction_dtt{i}(:,1:xlim_max)), 'Color', 'm', 'LineWidth', 3) 
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Frames (/3.2 ms)'), ylabel('\Delta\theta (\circ)'), title('\Delta\theta_T_o_p', 'FontSize', 25)

subplot(224), hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(1:xlim_max,nanmean(protraction_dtf{i}(:,1:xlim_max)),nanstd(protraction_dtf{i}(:,1:xlim_max)), 'Color', 'b', 'LineWidth', 3) 
        case 2
            errorbar(1:xlim_max,nanmean(protraction_dtf{i}(:,1:xlim_max)),nanstd(protraction_dtf{i}(:,1:xlim_max)), 'Color', 'c', 'LineWidth', 3) 
        case 3
            errorbar(1:xlim_max,nanmean(protraction_dtf{i}(:,1:xlim_max)),nanstd(protraction_dtf{i}(:,1:xlim_max)), 'Color', 'r', 'LineWidth', 3) 
        case 4
            errorbar(1:xlim_max,nanmean(protraction_dtf{i}(:,1:xlim_max)),nanstd(protraction_dtf{i}(:,1:xlim_max)), 'Color', 'm', 'LineWidth', 3) 
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Frames (/3.2 ms)'), ylabel('\Delta\theta (\circ)'), title('\Delta\theta_F_r_o_n_t', 'FontSize', 25)
legend({'Close Down','Far Down','Close Up','Far Up'},'Location', 'eastoutside'), suptitle([mouseName, ' ', sessionName],'FontSize',35)

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

%%
A = viewmtx(0,90);
x4d_prot = [data_prot, ones(size(data_prot,1),1)]';
x2d_prot = A*x4d_prot;
x4d_ret = [data_ret, ones(size(data_ret,1),1)]';
x2d_ret = A*x4d_ret;

h = figure; plot(x2d_prot(1,:), x2d_prot(2,:),'r.', 'MarkerSize',1); hold on; plot(x2d_ret(1,:), x2d_ret(2,:),'b.', 'MarkerSize',1);

title({[mouseName sessionName]; 'Free whisking theta'}), xlabel('Top-view theta'), ylabel('Front-view theta');

%%
A = viewmtx(90,0);
xlim = [-21 -19];
data_prot_res = data_prot( (data_prot(:,1)-mean(xlim)) < floor((max(xlim)-min(xlim))/2), : );
data_ret_res = data_ret( (data_ret(:,1)-mean(xlim)) < floor((max(xlim)-min(xlim))/2), : );
x4d_prot_res = [data_prot_res, ones(size(data_prot_res,1),1)]';
x2d_prot_res = A*x4d_prot_res;
x4d_ret_res = [data_ret_res, ones(size(data_ret_res,1),1)]';
x2d_ret_res = A*x4d_ret_res;

h = figure; plot(x2d_prot_res(1,:), x2d_prot_res(2,:),'r.', 'MarkerSize',1); hold on; plot(x2d_ret_res(1,:), x2d_ret_res(2,:),'b.', 'MarkerSize',1);

title({[mouseName sessionName]; 'Free whisking theta'; ['Top-view theta ' num2str(xlim(1)) ' ~ ' num2str(xlim(2))]}), xlabel('Front-view theta'), ylabel('Top-view kappa');

%%
% %% collapsing deltakappa and theta values into one matrix for easier averaging

protraction_dkt = cell(1,length(trial_types)); %delta kappa _t
protraction_dkf = cell(1,length(trial_types)); % delta kappa _f
protraction_tht = cell(1,length(trial_types)); % theta _t
protraction_thf = cell(1,length(trial_types)); % theta _f
protraction_dtt = cell(1,length(trial_types)); % delta theta _t
protraction_dtf = cell(1,length(trial_types)); % delta theta _f

for ti = 1 : length(trial_types)
    protraction_dkt{ti} = zeros(length(protraction_deltakappa_t{tind}),max(cellfun(@(x) length(x),protraction_deltakappa_t{tind})));
    protraction_dkf{ti} = zeros(size(protraction_dkt{ti}));
    protraction_tht{ti} = zeros(size(protraction_dkt{ti}));
    protraction_thf{ti} = zeros(size(protraction_dkt{ti}));
    protraction_dtt{ti} = zeros(size(protraction_dkt{ti}));
    protraction_dtf{ti} = zeros(size(protraction_dkt{ti}));
end
for tind = 1 : length(trial_types)
    for i = 1 : length(protraction_deltakappa_t{tind})
        protraction_dkt{tind}(i,1:length(protraction_deltakappa_t{tind}{i})) = protraction_deltakappa_t{tind}{i}; if length(protraction_deltakappa_t{tind}{i}) < size(protraction_dkt,2), protraction_dkt{tind}(i,length(protraction_deltakappa_t{tind}{i})+1:end) = NaN; end;
        protraction_dkf{tind}(i,1:length(protraction_deltakappa_f{tind}{i})) = protraction_deltakappa_f{tind}{i}; if length(protraction_deltakappa_f{tind}{i}) < size(protraction_dkt,2), protraction_dkf{tind}(i,length(protraction_deltakappa_f{tind}{i})+1:end) = NaN; end;
        protraction_tht{tind}(i,1:length(protraction_theta_t{tind}{i})) = protraction_theta_t{tind}{i}; if length(protraction_theta_t{tind}{i}) < size(protraction_dkt,2), protraction_tht{tind}(i,length(protraction_theta_t{tind}{i})+1:end) = NaN; end;
        protraction_thf{tind}(i,1:length(protraction_theta_f{tind}{i})) = protraction_theta_f{tind}{i}; if length(protraction_theta_f{tind}{i}) < size(protraction_dkt,2), protraction_thf{tind}(i,length(protraction_theta_f{tind}{i})+1:end) = NaN; end;
        protraction_dtt{tind}(i,1:length(protraction_deltatheta_t{tind}{i})) = protraction_deltatheta_t{tind}{i}; if length(protraction_deltatheta_t{tind}{i}) < size(protraction_dkt,2), protraction_dtt{tind}(i,length(protraction_deltatheta_t{tind}{i})+1:end) = NaN; end;
        protraction_dtf{tind}(i,1:length(protraction_deltatheta_f{tind}{i})) = protraction_deltatheta_f{tind}{i}; if length(protraction_deltatheta_f{tind}{i}) < size(protraction_dkt,2), protraction_dtf{tind}(i,length(protraction_deltatheta_f{tind}{i})+1:end) = NaN; end;        
    end
end

% %% averaging deltakappa and theta values altogether 
for i = 1 : length(trial_types)
    figure, subplot(221), errorbar(1:size(protraction_dkt{i},2),nanmean(protraction_dkt{i}),nanstd(protraction_dkt{i})), title('delta kappa top'), xlim([0 20])
    subplot(222), errorbar(1:size(protraction_dkf{i},2),nanmean(protraction_dkf{i}),nanstd(protraction_dkf{i})), title('delta kappa front'), xlim([0 20])
    subplot(223), errorbar(1:size(protraction_dtt{i},2),nanmean(protraction_dtt{i}),nanstd(protraction_dtt{i})), title('delta theta top'), xlim([0 20])
    subplot(224), errorbar(1:size(protraction_dtf{i},2),nanmean(protraction_dtf{i}),nanstd(protraction_dtf{i})), title('delta theta front'), xlim([0 20])
end


%%
figure,
subplot(221), title('delta kappa top'), hold all
for i = 1 : 4
    for j = 1 : length(deltakappa_t{j})
        switch i
            case 1 % close down
                plot(1:length(deltakappa_t{i}{j}), deltakappa_t{i}{j}, 'b.')
            case 2 % far down
                plot(1:length(deltakappa_t{i}{j}), deltakappa_t{i}{j}, 'c.')
            case 3 % close up    
                plot(1:length(deltakappa_t{i}{j}), deltakappa_t{i}{j}, 'r.')
            case 4 % far up    
                plot(1:length(deltakappa_t{i}{j}), deltakappa_t{i}{j}, 'm.')
        end    
    end
end


%%


