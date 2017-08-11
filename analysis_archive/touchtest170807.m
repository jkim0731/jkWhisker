% close all
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
    sessionName = sprintf('S%02d', sessions{str2double(mousensession(1))}{str2double(mousensession(2))}(str2double(mousensession(3))));
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
%             temp_whisking = setdiff(wt.pole_available_timepoints,wt.th_touch_frames);
%             temp_whisking = 1:wt.pole_available_timepoints(1)-1;
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
%                 whisking_diff_inds = [0, find(diff(temp_whisking) - 1), length(temp_whisking)];
%                 for cind = 1 : length(whisking_diff_inds) - 1
%                     temp_whisking_bout = temp_whisking(whisking_diff_inds(cind)+1:whisking_diff_inds(cind+1)); % in frames
%                     if length(temp_whisking_bout) > 2
                        temp_whisking_bout = 1:wt.pole_available_timepoints(1)-1; % free whisking before pole_up 2017/08/10
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
%                     end
%                 end
            end
        end              
    end
end

%% Touch angle (top) histogram

max_theta = max([max(cellfun(@(x) max(x),theta_t{1})), max(cellfun(@(x) max(x),theta_t{2})), max(cellfun(@(x) max(x),theta_t{3})), max(cellfun(@(x) max(x),theta_t{4}))]); % 4 for # of trial types
min_theta = min([min(cellfun(@(x) min(x),theta_t{1})), min(cellfun(@(x) min(x),theta_t{2})), min(cellfun(@(x) min(x),theta_t{3})), min(cellfun(@(x) min(x),theta_t{4}))]);
min_edge = floor(min_theta/5)*5; % 5 for angle division
max_edge = ceil(max_theta/5)*5;
hist_edge = min_edge:5:max_edge;
hist_values = zeros(length(trial_types), length(hist_edge)-1);
for i = 1 : length(trial_types)
    hist_values(i,:) = histcounts(cell2mat(theta_t{i}),hist_edge);
end
plot_xval = hist_edge(1:end-1) + 5/2 + base_angle;
figure, hold on;
for i = 1 : length(trial_types)
    switch i
        case 1
            plot(plot_xval, hist_values(i,:),'b-', 'LineWidth', 3)
        case 2
            plot(plot_xval, hist_values(i,:),'c-', 'LineWidth', 3)
        case 3
            plot(plot_xval, hist_values(i,:),'r-', 'LineWidth', 3)
        case 4
            plot(plot_xval, hist_values(i,:),'m-', 'LineWidth', 3)
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15;
legend({'Close Down','Far Down','Close Up','Far Up'}), xlabel('\theta_T_o_p (\circ)'), ylabel('# of frames'), title({[mouseName, ' ', sessionName]; 'Touch angle'},'FontSize',15)

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

%% Calculate average linear fit of free whisking that spanned over 10 degrees
figure, hold on,
for nn = 5:5:20
    slopes = [];
    for i = 1 : length(theta_whisking_t_prot)
        for j = 1 : length(theta_whisking_t_prot{i})
            if theta_whisking_t_prot{i}{j}(end) - theta_whisking_t_prot{i}{j}(1) > nn
                slopes = [slopes;polyfit(theta_whisking_t_prot{i}{j},theta_whisking_f_prot{i}{j},1)];
            end
        end
    end
    histogram(slopes(:,1))
end
ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';
legend({'>5 frames', '>10 frames', '>15 frames', '>20 frames'}), xlabel('Linear fit slope from \theta_T_o_p VS \theta_F_r_o_n_t'), ylabel('Counts'),
title({[mouseName, ' ', sessionName]; 'Free whisking theta linear fit histogram'}, 'FontSize',15)

%% Linear fitting curve and R2 value between (1) theta_top VS theta_front, (2) theta_top VS kappa_top, and
% quadratic fitting curve and R2 value between (3) theta_top VS kappa_front
% all during free whisking
theta_t_temp = [];
theta_f_temp = [];
kappa_t_temp = [];
kappa_f_temp = [];
for i = 1 : length(theta_whisking_t_prot)    
    theta_t_temp = [theta_t_temp, cell2mat(cellfun(@(x) x,theta_whisking_t_prot{i},'UniformOutput', false))];    
    theta_f_temp = [theta_f_temp, cell2mat(cellfun(@(x) x,theta_whisking_f_prot{i},'UniformOutput', false))];
    kappa_t_temp = [kappa_t_temp, cell2mat(cellfun(@(x) x,kappa_whisking_t_prot{i},'UniformOutput', false))];    
    kappa_f_temp = [kappa_f_temp, cell2mat(cellfun(@(x) x,kappa_whisking_f_prot{i},'UniformOutput', false))];    
end
for i = 1 : length(theta_whisking_t_ret)    
    theta_t_temp = [theta_t_temp, cell2mat(cellfun(@(x) x,theta_whisking_t_ret{i},'UniformOutput', false))];    
    theta_f_temp = [theta_f_temp, cell2mat(cellfun(@(x) x,theta_whisking_f_ret{i},'UniformOutput', false))];
    kappa_t_temp = [kappa_t_temp, cell2mat(cellfun(@(x) x,kappa_whisking_t_ret{i},'UniformOutput', false))];    
    kappa_f_temp = [kappa_f_temp, cell2mat(cellfun(@(x) x,kappa_whisking_f_ret{i},'UniformOutput', false))];    
end

nanind = isnan(theta_t_temp) + isnan(theta_f_temp) + isnan(kappa_t_temp) + isnan(kappa_f_temp);
theta_t_total = theta_t_temp(nanind == 0);
theta_f_total = theta_f_temp(nanind == 0);
kappa_t_total = kappa_t_temp(nanind == 0);
kappa_f_total = kappa_f_temp(nanind == 0);

% (1) theta_top VS theta_front
p1 = polyfit(theta_t_total,theta_f_total,1);
yfit = polyval(p1,theta_t_total);
yres = theta_f_total - yfit;
SSres = sum(yres.^2);
SStotal = (length(theta_f_total)-1) * var(theta_f_total);
rsq1 = 1 - SSres/SStotal;
% (2) theta_top VS kappa_top
p2 = polyfit(theta_t_total,kappa_t_total,1);
yfit = polyval(p2,theta_t_total);
yres = kappa_t_total - yfit;
SSres = sum(yres.^2);
SStotal = (length(kappa_t_total)-1) * var(kappa_t_total);
rsq2 = 1 - SSres/SStotal;
% (3) theta_top VS kappa_front
p3 = polyfit(theta_t_total,kappa_f_total,2);
yfit = polyval(p3,theta_t_total);
yres = kappa_f_total - yfit;
SSres = sum(yres.^2);
SStotal = (length(kappa_f_total)-1) * var(kappa_f_total);
rsq3 = 1 - SSres/SStotal;

%% Calculate average linear fit of free whisking that spanned over 10 degrees, and calculate average R2
slopes = [];
for i = 1 : length(theta_whisking_t_prot)
    for j = 1 : length(theta_whisking_t_prot{i})
        if theta_whisking_t_prot{i}{j}(end) - theta_whisking_t_prot{i}{j}(1) > 10
            slopes = [slopes;polyfit(theta_whisking_t_prot{i}{j},theta_whisking_f_prot{i}{j},1)];
        end
    end
end
avg_slope = mean(slopes);

r2 = [];
for i = 1 : length(theta_whisking_t_prot)
    for j = 1 : length(theta_whisking_t_prot{i})
        if theta_whisking_t_prot{i}{j}(end) - theta_whisking_t_prot{i}{j}(1) > 10
            
            r2 = [r2;polyfit(theta_whisking_t_prot{i}{j},theta_whisking_f_prot{i}{j},1)];
        end
    end
end



%% 3d plot of every points - dead
% 
% figure, hold all
% % for i = 1 : length(theta_whisking_t_prot)
% for i = 1
%     for j = 1 : length(theta_whisking_t_prot{i})
%         try
%             plot3(theta_whisking_t_prot{i}{j}+base_angle,theta_whisking_f_prot{i}{j},kappa_whisking_t_prot{i}{j},'r.');        
%         catch
%         end
%     end
% end
% 
% % for i = 1 : length(theta_whisking_t_ret)
% for i = 1
%     for j = 1 : length(theta_whisking_t_ret{i})    
%         try
%             plot3(theta_whisking_t_ret{i}{j}+base_angle,theta_whisking_f_ret{i}{j},kappa_whisking_t_ret{i}{j},'b.');
%         catch
%         end
%     end
% end
% 
% title({[mouseName sessionName]; 'Free whisking theta'}), xlabel('Top-view theta'), ylabel('Front-view theta'), zlabel('Top-view kappa')
%% Data reorganization for fast 3D view
data_prot = [];
data_ret = [];
for i = 1 : length(theta_whisking_t_prot)
% for i = 1
    for j = 1 : length(theta_whisking_t_prot{i})
        try
            data_prot = [data_prot; theta_whisking_t_prot{i}{j}'+base_angle,theta_whisking_f_prot{i}{j}',kappa_whisking_t_prot{i}{j}',kappa_whisking_f_prot{i}{j}'];        
        catch
        end
    end
end

for i = 1 : length(theta_whisking_t_ret)
% for i = 1
    for j = 1 : length(theta_whisking_t_ret{i})    
        try
            data_ret = [data_ret; theta_whisking_t_ret{i}{j}'+base_angle,theta_whisking_f_ret{i}{j}',kappa_whisking_t_ret{i}{j}',kappa_whisking_f_ret{i}{j}'];
        catch
        end
    end
end

%% fast 3D view
figure,
subplot(231)
A = viewmtx(0,90);
x4d_prot = [data_prot(:,1:3), ones(size(data_prot,1),1)]'; x2d_prot = A*x4d_prot;
x4d_ret = [data_ret(:,1:3), ones(size(data_ret,1),1)]'; x2d_ret = A*x4d_ret;
plot(x2d_prot(1,:), x2d_prot(2,:),'r.', 'MarkerSize',1); hold on; plot(x2d_ret(1,:), x2d_ret(2,:),'b.', 'MarkerSize',1);
xlabel('\theta_T_o_p'), ylabel('\theta_F_r_o_n_t'); title('\theta_T_o_p VS \theta_F_r_o_n_t')
ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';

subplot(232)
A = viewmtx(0,90);
x4d_prot = [data_prot(:,[1,3,4]), ones(size(data_prot,1),1)]'; x2d_prot = A*x4d_prot;
x4d_ret = [data_ret(:,[1,3,4]), ones(size(data_ret,1),1)]'; x2d_ret = A*x4d_ret;
plot(x2d_prot(1,:), x2d_prot(2,:),'r.', 'MarkerSize',1); hold on; plot(x2d_ret(1,:), x2d_ret(2,:),'b.', 'MarkerSize',1);
xlabel('\theta_T_o_p'), ylabel('\kappa_T_o_p'); title('\theta_T_o_p VS \kappa_T_o_p')
ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';

subplot(233)
A = viewmtx(0,90);
x4d_prot = [data_prot(:,[1,4,3]), ones(size(data_prot,1),1)]'; x2d_prot = A*x4d_prot;
x4d_ret = [data_ret(:,[1,4,3]), ones(size(data_ret,1),1)]'; x2d_ret = A*x4d_ret;
plot(x2d_prot(1,:), x2d_prot(2,:),'r.', 'MarkerSize',1); hold on; plot(x2d_ret(1,:), x2d_ret(2,:),'b.', 'MarkerSize',1);
xlabel('\theta_T_o_p'), ylabel('\kappa_F_r_o_n_t'); title('\theta_T_o_p VS \kappa_F_r_o_n_t')
ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';

subplot(234)
A = viewmtx(0,90);
x4d_prot = [data_prot(:,[2,4,3]), ones(size(data_prot,1),1)]'; x2d_prot = A*x4d_prot;
x4d_ret = [data_ret(:,[2,4,3]), ones(size(data_ret,1),1)]'; x2d_ret = A*x4d_ret;
plot(x2d_prot(1,:), x2d_prot(2,:),'r.', 'MarkerSize',1); hold on; plot(x2d_ret(1,:), x2d_ret(2,:),'b.', 'MarkerSize',1);
xlabel('\theta_F_r_o_n_t'), ylabel('\kappa_F_r_o_n_t'); title('\theta_F_r_o_n_t VS \kappa_F_r_o_n_t')
ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';

subplot(235)
A = viewmtx(0,90);
x4d_prot = [data_prot(:,[2,3,4]), ones(size(data_prot,1),1)]'; x2d_prot = A*x4d_prot;
x4d_ret = [data_ret(:,[2,3,4]), ones(size(data_ret,1),1)]'; x2d_ret = A*x4d_ret;
plot(x2d_prot(1,:), x2d_prot(2,:),'r.', 'MarkerSize',1); hold on; plot(x2d_ret(1,:), x2d_ret(2,:),'b.', 'MarkerSize',1);
xlabel('\theta_F_r_o_n_t'), ylabel('\kappa_T_o_p'); title('\theta_F_r_o_n_t VS \kappa_T_o_p')
ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';

subplot(236)
A = viewmtx(0,90);
x4d_prot = [data_prot(:,[3,4,1]), ones(size(data_prot,1),1)]'; x2d_prot = A*x4d_prot;
x4d_ret = [data_ret(:,[3,4,1]), ones(size(data_ret,1),1)]'; x2d_ret = A*x4d_ret;
plot(x2d_prot(1,:), x2d_prot(2,:),'r.', 'MarkerSize',1); hold on; plot(x2d_ret(1,:), x2d_ret(2,:),'b.', 'MarkerSize',1);
ylabel('\kappa_F_r_o_n_t'), xlabel('\kappa_T_o_p'); title('\kappa_T_o_p VS \kappa_F_r_o_n_t')
ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';

legend({'Protraction','Retraction'},'Location', 'eastoutside', 'FontSize', 15)
suptitle([mouseName, ' ', sessionName, ' Free whisking'], 'FontSize', 30)
%
%%
% A = viewmtx(90,0);
% xlim = [-21 -19];
% data_prot_res = data_prot( (data_prot(:,1)-mean(xlim)) < floor((max(xlim)-min(xlim))/2), : );
% data_ret_res = data_ret( (data_ret(:,1)-mean(xlim)) < floor((max(xlim)-min(xlim))/2), : );
% x4d_prot_res = [data_prot_res, ones(size(data_prot_res,1),1)]';
% x2d_prot_res = A*x4d_prot_res;
% x4d_ret_res = [data_ret_res, ones(size(data_ret_res,1),1)]';
% x2d_ret_res = A*x4d_ret_res;
% 
% h = figure; plot(x2d_prot_res(1,:), x2d_prot_res(2,:),'r.', 'MarkerSize',1); hold on; plot(x2d_ret_res(1,:), x2d_ret_res(2,:),'b.', 'MarkerSize',1);
% 
% title({[mouseName sessionName]; 'Free whisking theta'; ['Top-view theta ' num2str(xlim(1)) ' ~ ' num2str(xlim(2))]}), xlabel('Front-view theta'), ylabel('Top-view kappa');



