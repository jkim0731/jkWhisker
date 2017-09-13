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
ic_t = cell(1,length(trial_types));
ic_f = cell(1,length(trial_types));

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
            temp_whisking = 1:wt.pole_available_timepoints(1)-1;
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
                        ic_t{ti}{end+1} = wt.intersect_coord(temp_ind,1);
                        ic_f{ti}{end+1} = wt.intersect_coord(temp_ind,2);                        
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

%% Touch angle (front) histogram

max_theta = max([max(cellfun(@(x) max(x),theta_f{1})), max(cellfun(@(x) max(x),theta_f{2})), max(cellfun(@(x) max(x),theta_f{3})), max(cellfun(@(x) max(x),theta_f{4}))]); % 4 for # of trial types
min_theta = min([min(cellfun(@(x) min(x),theta_f{1})), min(cellfun(@(x) min(x),theta_f{2})), min(cellfun(@(x) min(x),theta_f{3})), min(cellfun(@(x) min(x),theta_f{4}))]);
min_edge = floor(min_theta/5)*5; % 5 for angle division
max_edge = ceil(max_theta/5)*5;
hist_edge = min_edge:5:max_edge;
hist_values = zeros(length(trial_types), length(hist_edge)-1);
for i = 1 : length(trial_types)
    hist_values(i,:) = histcounts(cell2mat(theta_f{i}),hist_edge);
end
plot_xval = hist_edge(1:end-1) + 5/2;
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
legend({'Close Down','Far Down','Close Up','Far Up'},'Location','west'), xlabel('\theta_F_r_o_n_t (\circ)'), ylabel('# of frames'), title({[mouseName, ' ', sessionName]; 'Touch angle'},'FontSize',15)

%% collapse and average only those from protraction touch
% %% assign deltakappa and theta. deltakappa is a true deltakappa, calculated by subtracting the kappa value of one frame before the first touch frame in each chunk
protraction_deltakappa_t = cell(1,length(trial_types));
protraction_deltakappa_f = cell(1,length(trial_types));
protraction_theta_t = cell(1,length(trial_types));
protraction_theta_f = cell(1,length(trial_types));
protraction_deltatheta_t = cell(1,length(trial_types));
protraction_deltatheta_f = cell(1,length(trial_types));
ic_t_prot = cell(1,length(trial_types));
ic_f_prot = cell(1,length(trial_types));
ic_t_ret = cell(1,length(trial_types));
ic_f_ret = cell(1,length(trial_types));

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
                    && wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1)+1) > wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1)) ... % and theta_top is increasing when the touch started
                    && wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1)) > wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1)-1) % and theta_top is increasing when the touch started 
                    
                    protraction_deltakappa_t{ti}{end+1} = wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1));
                    protraction_deltakappa_f{ti}{end+1} = wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1));
                    protraction_theta_t{ti}{end+1} = wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                    protraction_theta_f{ti}{end+1} = wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                    protraction_deltatheta_t{ti}{end+1} = wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1));
                    protraction_deltatheta_f{ti}{end+1} = wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1));
                    ic_t_prot{ti}{end+1} = wt.intersect_coord(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)),1);
                    ic_f_prot{ti}{end+1} = wt.intersect_coord(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)),2);
                end
            end
        end
    end
end

% %% collapsing deltakappa and theta values into one matrix for easier averaging

protraction_dkt = cell(1,length(trial_types)); % delta kappa _t
protraction_dkf = cell(1,length(trial_types)); % delta kappa _f
protraction_tht = cell(1,length(trial_types)); % theta _t
protraction_thf = cell(1,length(trial_types)); % theta _f
protraction_dtt = cell(1,length(trial_types)); % delta theta _t
protraction_dtf = cell(1,length(trial_types)); % delta theta _f
dic_t_prot = cell(1,length(trial_types)); % delta intersection coordinate _t
dic_f_prot = cell(1,length(trial_types));

for ti = 1 : length(trial_types)
    protraction_dkt{ti} = zeros(length(protraction_deltakappa_t{ti}),max(cellfun(@(x) length(x),protraction_deltakappa_t{ti})));
    protraction_dkf{ti} = zeros(size(protraction_dkt{ti}));
    protraction_tht{ti} = zeros(size(protraction_dkt{ti}));
    protraction_thf{ti} = zeros(size(protraction_dkt{ti}));
    protraction_dtt{ti} = zeros(size(protraction_dkt{ti}));
    protraction_dtf{ti} = zeros(size(protraction_dkt{ti}));
    dic_t_prot{ti} = zeros(size(protraction_dkt{ti}));
    dic_f_prot{ti} = zeros(size(protraction_dkt{ti}));
end
for tind = 1 : length(trial_types)
    for i = 1 : length(protraction_deltakappa_t{tind})
%         protraction_dkt{tind}(i,1:length(protraction_deltakappa_t{tind}{i})) = protraction_deltakappa_t{tind}{i}; if length(protraction_deltakappa_t{tind}{i}) < size(protraction_dkt{tind},2), protraction_dkt{tind}(i,length(protraction_deltakappa_t{tind}{i})+1:end) = NaN; end;
%         protraction_dkf{tind}(i,1:length(protraction_deltakappa_f{tind}{i})) = protraction_deltakappa_f{tind}{i}; if length(protraction_deltakappa_f{tind}{i}) < size(protraction_dkt{tind},2), protraction_dkf{tind}(i,length(protraction_deltakappa_f{tind}{i})+1:end) = NaN; end;
%         protraction_tht{tind}(i,1:length(protraction_theta_t{tind}{i})) = protraction_theta_t{tind}{i}; if length(protraction_theta_t{tind}{i}) < size(protraction_dkt{tind},2), protraction_tht{tind}(i,length(protraction_theta_t{tind}{i})+1:end) = NaN; end;
%         protraction_thf{tind}(i,1:length(protraction_theta_f{tind}{i})) = protraction_theta_f{tind}{i}; if length(protraction_theta_f{tind}{i}) < size(protraction_dkt{tind},2), protraction_thf{tind}(i,length(protraction_theta_f{tind}{i})+1:end) = NaN; end;
%         protraction_dtt{tind}(i,1:length(protraction_deltatheta_t{tind}{i})) = protraction_deltatheta_t{tind}{i}; if length(protraction_deltatheta_t{tind}{i}) < size(protraction_dkt{tind},2), protraction_dtt{tind}(i,length(protraction_deltatheta_t{tind}{i})+1:end) = NaN; end;
%         protraction_dtf{tind}(i,1:length(protraction_deltatheta_f{tind}{i})) = protraction_deltatheta_f{tind}{i}; if length(protraction_deltatheta_f{tind}{i}) < size(protraction_dkt{tind},2), protraction_dtf{tind}(i,length(protraction_deltatheta_f{tind}{i})+1:end) = NaN; end;        
%         dic_t_prot{tind}(i,1:length(ic_t_prot{tind}{i})) = ic_t_prot{tind}{i}; if length(ic_t_prot{tind}{i}) < size(protraction_dkt{tind},2), dic_t_prot{tind}(i,length(ic_t_prot{tind}{i})+1:end) = NaN; end;
%         dic_f_prot{tind}(i,1:length(ic_f_prot{tind}{i})) = ic_f_prot{tind}{i}; if length(ic_f_prot{tind}{i}) < size(protraction_dkt{tind},2), dic_f_prot{tind}(i,length(ic_f_prot{tind}{i})+1:end) = NaN; end;
        protraction_dkt{tind}(i,1:length(protraction_deltakappa_t{tind}{i})) = protraction_deltakappa_t{tind}{i}; 
        protraction_dkf{tind}(i,1:length(protraction_deltakappa_f{tind}{i})) = protraction_deltakappa_f{tind}{i}; 
        protraction_tht{tind}(i,1:length(protraction_theta_t{tind}{i})) = protraction_theta_t{tind}{i}; 
        protraction_thf{tind}(i,1:length(protraction_theta_f{tind}{i})) = protraction_theta_f{tind}{i}; 
        protraction_dtt{tind}(i,1:length(protraction_deltatheta_t{tind}{i})) = protraction_deltatheta_t{tind}{i}; 
        protraction_dtf{tind}(i,1:length(protraction_deltatheta_f{tind}{i})) = protraction_deltatheta_f{tind}{i}; 
        dic_t_prot{tind}(i,1:length(ic_t_prot{tind}{i})) = ic_t_prot{tind}{i} - ic_t_prot{tind}{i}(1); 
        dic_f_prot{tind}(i,1:length(ic_f_prot{tind}{i})) = ic_f_prot{tind}{i} - ic_f_prot{tind}{i}(1); 
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

%% ************************************************************************

%% Calculate average linear fit of free whisking that spanned over 10 degrees

%% ************************************************************************
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
avg_slope = nanmean(slopes);

r2 = [];
f2 = [];
c2_f = {};
c2_t = {};
x2 = [];
v2 = [];
l2 = [];
for i = 1 : length(theta_whisking_t_prot)
    for j = 1 : length(theta_whisking_t_prot{i})
        if theta_whisking_t_prot{i}{j}(end) - theta_whisking_t_prot{i}{j}(1) > 10
            fun = @(x) sum((theta_whisking_f_prot{i}{j}-theta_whisking_t_prot{i}{j}*avg_slope(1)-x).^2);
            [x,fval] = fminsearch(fun,0);
            if ~isnan(fval)
                f2 = [f2;fval];
                rsq = 1 - fval/((length(theta_whisking_f_prot{i}{j})-1) * var(theta_whisking_f_prot{i}{j}));
                v2 = [v2;var(theta_whisking_f_prot{i}{j})];
                l2 = [l2;length(theta_whisking_f_prot{i}{j})];
                r2 = [r2;rsq];
                c2_f{end+1} = theta_whisking_f_prot{i}{j};
                c2_t{end+1} = theta_whisking_t_prot{i}{j};
                x2 = [x2;x];
            end
        end
    end
end
%% make a figure
figure, histogram(r2,[0:0.05:1]), xlabel('R-squared'), ylabel('Counts'), ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';
title({[mouseName, ' ', sessionName]; 'Free whisking theta linear fit histogram'}, 'FontSize',15)

%% test
iii = 909;
length(c2_t{iii})
figure, plot(c2_t{iii},c2_f{iii},'k.','MarkerSize',20), hold on, plot(c2_t{iii}, polyval([avg_slope(1),x2(iii)],c2_t{iii}),'r-','LineWidth',3)
ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';
xlabel('\theta_T_o_p'), ylabel('\theta_F_r_o_n_t'),
title({[mouseName, ' ', sessionName]; ['Free whisking protraction e.g. #' num2str(iii)]}, 'FontSize',15)
%% Make all negative values to zero, normalize, and plot the histogram again
r2(r2<0) = 0;
figure, histogram(r2,[0:0.05:1],'Normalization','probability'), xlabel('R-squared'), ylabel('Proportion'), ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';
title({[mouseName, ' ', sessionName]; 'Free whisking theta linear fit histogram'}, 'FontSize',15)

%% re-calculating normalized theta_front based on the linear fitting to theta_top
prot_dtheta_frontfit = cell(1,length(trial_types));
for ti = 1 : length(trial_types)    
    prot_dtheta_frontfit{ti} = zeros(size(protraction_dkt{ti}));
end

for i = 1 : length(trial_types)
    for j = 1 : length(protraction_theta_t{i})
        prot_dtheta_frontfit{i}(j,1:length(protraction_theta_f{i}{j})) = protraction_theta_f{i}{j} - protraction_theta_t{i}{j}*avg_slope(1);
        prot_dtheta_frontfit{i}(j,1:length(protraction_theta_f{i}{j})) = prot_dtheta_frontfit{i}(j,1:length(protraction_theta_f{i}{j})) - prot_dtheta_frontfit{i}(j,1);        
    end
end

%% plot it
xlim_max = 15;
figure, hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(0:xlim_max,nanmean(prot_dtheta_frontfit{i}(:,1:xlim_max+1)),nanstd(prot_dtheta_frontfit{i}(:,1:xlim_max+1)), 'Color', 'b', 'LineWidth', 3) 
        case 2
            errorbar(0:xlim_max,nanmean(prot_dtheta_frontfit{i}(:,1:xlim_max+1)),nanstd(prot_dtheta_frontfit{i}(:,1:xlim_max+1)), 'Color', 'c', 'LineWidth', 3) 
        case 3
            errorbar(0:xlim_max,nanmean(prot_dtheta_frontfit{i}(:,1:xlim_max+1)),nanstd(prot_dtheta_frontfit{i}(:,1:xlim_max+1)), 'Color', 'r', 'LineWidth', 3) 
        case 4
            errorbar(0:xlim_max,nanmean(prot_dtheta_frontfit{i}(:,1:xlim_max+1)),nanstd(prot_dtheta_frontfit{i}(:,1:xlim_max+1)), 'Color', 'm', 'LineWidth', 3) 
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Frames (/3.2 ms)'), ylabel('\Delta\theta (\circ)'), title({[mouseName, ' ', sessionName];'\Delta\theta_F_r_o_n_t'}, 'FontSize', 15)
legend({'Close Down','Far Down','Close Up','Far Up'})
plot(0:xlim_max,zeros(1,xlim_max+1),'k--', 'LineWidth', 2)

%% re-calculating normalized kappa_front based on the quadratic fitting to theta_top
curves = [];
for i = 1 : length(theta_whisking_t_prot)
    for j = 1 : length(theta_whisking_t_prot{i})
        if theta_whisking_t_prot{i}{j}(end) - theta_whisking_t_prot{i}{j}(1) > 10 && length(theta_whisking_t_prot{i}{j}) > 2
            curves = [curves;polyfit(theta_whisking_t_prot{i}{j},kappa_whisking_f_prot{i}{j},2)];
        end
    end
end
avg_curve = nanmean(curves);

r2 = [];
f2 = [];
c2_f = {};
c2_t = {};
x2 = [];
v2 = [];
l2 = [];
for i = 1 : length(theta_whisking_t_prot)
    for j = 1 : length(theta_whisking_t_prot{i})
        if theta_whisking_t_prot{i}{j}(end) - theta_whisking_t_prot{i}{j}(1) > 10
            fun = @(x) sum((kappa_whisking_f_prot{i}{j}- avg_curve(1)*(theta_whisking_t_prot{i}{j}.^2) - avg_curve(2)*theta_whisking_t_prot{i}{j} - x).^2);
            [x,fval] = fminsearch(fun,0);
            if ~isnan(fval)
                f2 = [f2;fval];
                rsq = 1 - fval/((length(kappa_whisking_f_prot{i}{j})-1) * var(kappa_whisking_f_prot{i}{j}));
                v2 = [v2;var(kappa_whisking_f_prot{i}{j})];
                l2 = [l2;length(kappa_whisking_f_prot{i}{j})];
                r2 = [r2;rsq];
                c2_f{end+1} = kappa_whisking_f_prot{i}{j};
                c2_t{end+1} = kappa_whisking_t_prot{i}{j};
                x2 = [x2;x];
            end
        end
    end
end
%% make a figure
figure, histogram(r2), xlabel('R-squared'), ylabel('Counts'), ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';
title({[mouseName, ' ', sessionName]; 'Free whisking \kappa_F_r_o_n_t quadratic fit histogram'}, 'FontSize',15)

%%
r2(r2<0) = 0;
figure, histogram(r2,[0:0.05:1],'Normalization','probability'), xlabel('R-squared'), ylabel('Proportion'), ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';
title({[mouseName, ' ', sessionName]; 'Free whisking \kappa_F_r_o_n_t quadratic fit histogram'}, 'FontSize',15)

%% 
rsq_quad = [];
for i = 1 : length(theta_whisking_t_prot)
    for j = 1 : length(theta_whisking_t_prot{i})
        yfit = polyval(p3,theta_whisking_t_prot{i}{j});
        yres = kappa_whisking_f_prot{i}{j} - yfit;
        SSres = sum(yres.^2);
        SStotal = (length(kappa_whisking_f_prot{i}{j})-1) * var(kappa_whisking_f_prot{i}{j});
        rsq_quad = [rsq_quad; 1 - SSres/SStotal];
    end
end

% rsq_quad(rsq_quad<0) = -0.05;
%%
figure, histogram(rsq_quad,-100:1,'Normalization','probability'), xlabel('R-squared'), ylabel('Proportion'), ax = gca; ax.FontSize = 15; ax.FontWeight = 'bold'; ax.LineWidth = 3; ax.Box = 'off';
title({[mouseName, ' ', sessionName]; 'Free whisking \kappa_F_r_o_n_t quadratic fit histogram';'Fixed constant'}, 'FontSize',15)

%% re-calculating normalized kappa_front based on the quadratic fitting to theta_top
prot_dkappa_frontfit = cell(1,length(trial_types));
for ti = 1 : length(trial_types)    
    prot_dkappa_frontfit{ti} = zeros(size(protraction_dkt{ti}));
end

for i = 1 : length(trial_types)
    for j = 1 : length(protraction_theta_t{i})
        prot_dkappa_frontfit{i}(j,1:length(protraction_deltakappa_f{i}{j})) = protraction_deltakappa_f{i}{j} - avg_curve(1)*(protraction_theta_t{i}{j}.^2) - avg_curve(2)*protraction_theta_t{i}{j} - x;
        prot_dkappa_frontfit{i}(j,1:length(protraction_deltakappa_f{i}{j})) = prot_dkappa_frontfit{i}(j,1:length(protraction_deltakappa_f{i}{j})) - prot_dkappa_frontfit{i}(j,1);
    end
end

%% plot it
xlim_max = 15;
figure, hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(0:xlim_max,nanmean(prot_dkappa_frontfit{i}(:,1:xlim_max+1)),nanstd(prot_dkappa_frontfit{i}(:,1:xlim_max+1)), 'Color', 'b', 'LineWidth', 3) 
        case 2
            errorbar(0:xlim_max,nanmean(prot_dkappa_frontfit{i}(:,1:xlim_max+1)),nanstd(prot_dkappa_frontfit{i}(:,1:xlim_max+1)), 'Color', 'c', 'LineWidth', 3) 
        case 3
            errorbar(0:xlim_max,nanmean(prot_dkappa_frontfit{i}(:,1:xlim_max+1)),nanstd(prot_dkappa_frontfit{i}(:,1:xlim_max+1)), 'Color', 'r', 'LineWidth', 3) 
        case 4
            errorbar(0:xlim_max,nanmean(prot_dkappa_frontfit{i}(:,1:xlim_max+1)),nanstd(prot_dkappa_frontfit{i}(:,1:xlim_max+1)), 'Color', 'm', 'LineWidth', 3) 
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Frames (/3.2 ms)'), ylabel('\Delta\kappa (/mm)'), title({[mouseName, ' ', sessionName];'\Delta\kappa_F_r_o_n_t'}, 'FontSize', 15)
legend({'Close Down','Far Down','Close Up','Far Up'})
plot(0:xlim_max,zeros(1,xlim_max+1),'k--', 'LineWidth', 2)

%% ************************************************************************

%% change of intersection coordinate for sliding detection

%% ************************************************************************
% first, for top view intersection coordinate, histogram (each protraction
% touch, max intersection coord change)

figure, 
subplot(121), hold all
for i = 1 : length(trial_types)
    max_dic_t_prot = -min(dic_t_prot{i},[],2);
    histogram(max_dic_t_prot,0:0.5:20);
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Maximum top-view protraction (AU)'), ylabel('Counts'), title({[mouseName, ' ', sessionName];'Maximum intersection coordinate protraction'; '(top-view)'}, 'FontSize', 15)
% legend({'Close Down','Far Down','Close Up','Far Up'})

subplot(122), hold all
for i = 1 : length(trial_types)
    if i < 3 % Down pole
        max_dic_f_prot = [-min(dic_f_prot{i},[],2);-max(dic_f_prot{i},[],2)-0.25];
        histogram(max_dic_f_prot,-20:0.5:20);
    else % Up pole
        max_dic_f_prot = [max(dic_f_prot{i},[],2);min(dic_f_prot{i},[],2)-0.25];
        histogram(max_dic_f_prot,-20:0.5:20);
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Maximum front-view change (AU)'), ylabel('Counts'), title({[mouseName, ' ', sessionName];'Maximum intersection coordinate change'; '(front-view)'}, 'FontSize', 15)
legend({'Close Down','Far Down','Close Up','Far Up'})

%% change of intersection coordinate for sliding detection
% first, for top view intersection coordinate, histogram (each protraction
% touch, max intersection coord change) - Normalized

figure, 
subplot(131), hold all
max_dic_t_prot = [];
for i = 1 : length(trial_types)
    max_dic_t_prot = [max_dic_t_prot;-min(dic_t_prot{i}(:,1:5),[],2)];    
%     max_dic_t_prot = [max_dic_t_prot;-min(dic_t_prot{i},[],2)];    
end
% histogram(max_dic_t_prot,0:0.5:20,'Normalization','probability');
histogram(max_dic_t_prot,0:0.5:20,'Normalization','cdf'), ylim([0 1]);
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Maximum top-view protraction (AU)'), ylabel('Cumulative proportion'), title({[mouseName, ' ', sessionName];'Maximum intersection coordinate protraction'; '(top-view)'}, 'FontSize', 15)

subplot(132), hold all
max_dic_f_prot = [];
for i = 1 : length(trial_types)
    if i < 3
        max_dic_f_prot = [max_dic_f_prot;-min(dic_f_prot{i}(:,1:5),[],2)];
%         max_dic_f_prot = [max_dic_f_prot;-min(dic_f_prot{i},[],2)];
    else
        max_dic_f_prot = [max_dic_f_prot;max(dic_f_prot{i}(:,1:5),[],2)];
%         max_dic_f_prot = [max_dic_f_prot;max(dic_f_prot{i},[],2)];
    end    
end
% histogram(max_dic_f_prot,0:0.5:20,'Normalization','probability');
histogram(max_dic_f_prot,0:0.5:20,'Normalization','cdf');
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Maximum front-view change (AU)'), ylabel('Cumulative proportion'), title({[mouseName, ' ', sessionName];'Maximum intersection coordinate change'; '(front-view)'}, 'FontSize', 15)

subplot(133), hold all
min_dic_f_prot = [];
for i = 1 : length(trial_types)
    if i < 3
        min_dic_f_prot = [min_dic_f_prot;-max(dic_f_prot{i}(:,1:5),[],2)];
%         min_dic_f_prot = [min_dic_f_prot;-max(dic_f_prot{i},[],2)];
    else
        min_dic_f_prot = [min_dic_f_prot;min(dic_f_prot{i}(:,1:5),[],2)];
%         min_dic_f_prot = [min_dic_f_prot;min(dic_f_prot{i},[],2)];
    end    
end
% histogram(min_dic_f_prot,-20:0.5:0,'Normalization','probability');
histogram(min_dic_f_prot,-20:0.5:0,'Normalization','cdf');
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Maximum front-view reverse change (AU)'), ylabel('Cumulative proportion'), title({[mouseName, ' ', sessionName];'Maximum intersection coordinate reverse change'; '(front-view)'}, 'FontSize', 15)

%% Correlation between top-view and front-view change
% absolute values

figure, subplot(121), hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            for j = 1 : length(ic_t_prot{i})
                plot(ic_t_prot{i}{j}, ic_f_prot{i}{j}, 'b.');
            end
        case 2
            for j = 1 : length(ic_t_prot{i})
                plot(ic_t_prot{i}{j}, ic_f_prot{i}{j}, 'c.');
            end
        case 3
            for j = 1 : length(ic_t_prot{i})
                plot(ic_t_prot{i}{j}, ic_f_prot{i}{j}, 'r.');
            end
        case 4
            for j = 1 : length(ic_t_prot{i})
                plot(ic_t_prot{i}{j}, ic_f_prot{i}{j}, 'm.');
            end
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Top-view protraction (AU)'), ylabel('Front-view protraction (AU)'), title({[mouseName, ' ', sessionName];'Change in touch intersection coordinates (Absolute)'}, 'FontSize', 15)

% delta values
subplot(122), hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            plot(-dic_t_prot{i}(:), dic_f_prot{i}(:), 'b.');
        case 2
            plot(-dic_t_prot{i}(:), dic_f_prot{i}(:), 'c.');
        case 3
            plot(-dic_t_prot{i}(:), dic_f_prot{i}(:), 'r.');
        case 4
            plot(-dic_t_prot{i}(:), dic_f_prot{i}(:), 'm.');
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Top-view protraction (AU)'), ylabel('Front-view protraction (AU)'), title({[mouseName, ' ', sessionName];'Change in touch intersection coordinates (Relative)'}, 'FontSize', 15)
legend({'Close Down','Far Down','Close Up','Far Up'})

%% change of intersection coordinate for sliding detection
% Histogram from different parts in the data

max_dic_t_prot = [];
for i = 1 : length(trial_types)
    max_dic_t_prot = [max_dic_t_prot;-min(dic_t_prot{i}(:,1:5),[],2)];    
%     max_dic_t_prot = -min(dic_t_prot{i},[],2);    
end

max_dic_f_prot = [];
for i = 1 : length(trial_types)
    if i < 3
        max_dic_f_prot = [max_dic_f_prot;-min(dic_f_prot{i}(:,1:5),[],2)];
%         max_dic_f_prot = [max_dic_f_prot;-min(dic_f_prot{i},[],2)];
    else
        max_dic_f_prot = [max_dic_f_prot;max(dic_f_prot{i}(:,1:5),[],2)];
%         max_dic_f_prot = [max_dic_f_prot;max(dic_f_prot{i},[],2)];
    end    
end

min_dic_f_prot = [];
for i = 1 : length(trial_types)
    if i < 3
        min_dic_f_prot = [min_dic_f_prot;-max(dic_f_prot{i}(:,1:5),[],2)];
%         min_dic_f_prot = [min_dic_f_prot;-max(dic_f_prot{i},[],2)];
    else
        min_dic_f_prot = [min_dic_f_prot;min(dic_f_prot{i}(:,1:5),[],2)];
%         min_dic_f_prot = [min_dic_f_prot;min(dic_f_prot{i},[],2)];
    end    
end

%%
hist_ind1 = find(min_dic_f_prot<-0.5);
hist_ind2 = find(min_dic_f_prot>=-0.5);
max_dic_t_prot_hist1 = max_dic_t_prot(hist_ind1);
max_dic_f_prot_hist1 = max_dic_f_prot(hist_ind1);
min_dic_f_prot_hist1 = min_dic_f_prot(hist_ind1);
max_dic_t_prot_hist2 = max_dic_t_prot(hist_ind2);
max_dic_f_prot_hist2 = max_dic_f_prot(hist_ind2);
min_dic_f_prot_hist2 = min_dic_f_prot(hist_ind2);

figure, 
subplot(131), hold all
histogram(max_dic_t_prot_hist1,0:0.5:20);
histogram(max_dic_t_prot_hist2,0:0.5:20);
% histogram(max_dic_t_prot_hist,0:0.5:20,'Normalization','cdf');
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Maximum top-view protraction (AU)'), ylabel('Count'), title({[mouseName, ' ', sessionName];'Maximum intersection coordinate protraction'; '(top-view)'}, 'FontSize', 15)
% legend({'Close Down','Far Down','Close Up','Far Up'})

subplot(132), hold all
histogram(max_dic_f_prot_hist1,0:0.5:20);
histogram(max_dic_f_prot_hist2,0:0.5:20);
% histogram(max_dic_f_prot_hist,0:0.5:20,'Normalization','cdf');
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Maximum front-view change (AU)'), ylabel('Count'), title({[mouseName, ' ', sessionName];'Maximum intersection coordinate change'; '(front-view)'}, 'FontSize', 15)

subplot(133), hold all
histogram(min_dic_f_prot_hist1,-20:0.5:0);
histogram(min_dic_f_prot_hist2,-20:0.5:0);
% histogram(min_dic_f_prot_hist,-20:0.5:0,'Normalization','cdf');
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Maximum front-view reverse change (AU)'), ylabel('Count'), title({[mouseName, ' ', sessionName];'Maximum intersection coordinate reverse change'; '(front-view)'}, 'FontSize', 15)

%% ************************************************************************

%% time-series data from selected touches

%% ************************************************************************
max_dic_t_prot = cell(1,length(trial_types));
max_dic_f_prot = cell(1,length(trial_types));
% min_dic_f_prot = cell(1,length(trial_types));
ts_ind = cell(1,length(trial_types));
for i = 1 : length(trial_types)
    max_dic_t_prot{i} = -min(dic_t_prot{i}(:,1:5),[],2);
    if i < 3
        max_dic_f_prot{i} = -min(dic_f_prot{i}(:,1:5),[],2);
    else
        max_dic_f_prot{i} = max(dic_f_prot{i}(:,1:5),[],2);
    end      
%     ts_ind{i} = find(max_dic_f_prot{i} > 1.5);
%     if i < 3
%         min_dic_f_prot{i} = -max(dic_f_prot{i}(:,1:5),[],2);
%     else
%         min_dic_f_prot{i} = min(dic_f_prot{i}(:,1:5),[],2);
%     end  
%     ts_ind{i} = find(min_dic_f_prot{i} > -0.5);
%     ts_ind{i} = 1:length(max_dic_t_prot{i});
    ts_ind{i} = intersect(find(max_dic_t_prot{i} > 2), find(max_dic_f_prot{i} > 1.5));
end

prot_theta_f_selected = cell(1,length(trial_types));
prot_theta_t_selected = cell(1,length(trial_types));
prot_kappa_f_selected = cell(1,length(trial_types));
prot_kappa_t_selected = cell(1,length(trial_types));

for ti = 1 : length(trial_types)    
    prot_theta_f_selected{ti} = zeros(size(protraction_dkt{ti}));
    prot_theta_t_selected{ti} = zeros(size(protraction_dkt{ti}));
    prot_kappa_f_selected{ti} = zeros(size(protraction_dkt{ti}));
    prot_kappa_t_selected{ti} = zeros(size(protraction_dkt{ti}));
end

for i = 1 : length(trial_types)
    for j = 1 : length(ts_ind{i})
        tempind = ts_ind{i}(j);
        prot_theta_f_selected{i}(tempind,1:length(protraction_theta_f{i}{tempind})) = protraction_theta_f{i}{tempind} - protraction_theta_t{i}{tempind}*avg_slope(1);
        prot_theta_f_selected{i}(tempind,1:length(protraction_theta_f{i}{tempind})) = prot_theta_f_selected{i}(tempind,1:length(protraction_theta_f{i}{tempind})) - prot_theta_f_selected{i}(tempind,1);        
        prot_theta_t_selected{i}(tempind,1:length(protraction_theta_t{i}{tempind})) = protraction_theta_t{i}{tempind} - protraction_theta_t{i}{tempind}(1);
        prot_kappa_f_selected{i}(tempind,1:length(protraction_deltakappa_f{i}{tempind})) = protraction_deltakappa_f{i}{tempind} - avg_curve(1)*(protraction_theta_t{i}{tempind}.^2) - avg_curve(2)*protraction_theta_t{i}{tempind}; 
        prot_kappa_f_selected{i}(tempind,1:length(protraction_deltakappa_f{i}{tempind})) = prot_kappa_f_selected{i}(tempind,1:length(protraction_deltakappa_f{i}{tempind})) - prot_kappa_f_selected{i}(tempind,1);
        prot_kappa_t_selected{i}(tempind,1:length(protraction_deltakappa_t{i}{tempind})) = protraction_deltakappa_t{i}{tempind} - protraction_deltakappa_t{i}{tempind}(1);
    end
end
%% errorbar plot
xlim_max = 15;
figure, suptitle({[mouseName, ' ', sessionName, ' Sliding touch'];' '}, 'FontSize', 30)
subplot(221),hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(0:xlim_max,nanmean(prot_kappa_t_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_kappa_t_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'b', 'LineWidth', 3) 
        case 2
            errorbar(0:xlim_max,nanmean(prot_kappa_t_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_kappa_t_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'c', 'LineWidth', 3) 
        case 3
            errorbar(0:xlim_max,nanmean(prot_kappa_t_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_kappa_t_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'r', 'LineWidth', 3) 
        case 4
            errorbar(0:xlim_max,nanmean(prot_kappa_t_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_kappa_t_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'm', 'LineWidth', 3) 
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Frames (/3.2 ms)'), ylabel('\Delta\kappa (/mm)'), title('\Delta\kappa_T_o_p', 'FontSize', 15), ylim([-0.002 0.001])
plot(0:xlim_max,zeros(1,xlim_max+1),'k--', 'LineWidth', 2)

subplot(222),hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(0:xlim_max,nanmean(prot_kappa_f_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_kappa_f_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'b', 'LineWidth', 3) 
        case 2
            errorbar(0:xlim_max,nanmean(prot_kappa_f_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_kappa_f_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'c', 'LineWidth', 3) 
        case 3
            errorbar(0:xlim_max,nanmean(prot_kappa_f_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_kappa_f_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'r', 'LineWidth', 3) 
        case 4
            errorbar(0:xlim_max,nanmean(prot_kappa_f_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_kappa_f_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'm', 'LineWidth', 3) 
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Frames (/3.2 ms)'), ylabel('\Delta\kappa (/mm)'), title('\Delta\kappa_F_r_o_n_t', 'FontSize', 15), ylim([-0.0005 0.0005])
plot(0:xlim_max,zeros(1,xlim_max+1),'k--', 'LineWidth', 2)

subplot(223),hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(0:xlim_max,nanmean(prot_theta_t_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_theta_t_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'b', 'LineWidth', 3) 
        case 2
            errorbar(0:xlim_max,nanmean(prot_theta_t_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_theta_t_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'c', 'LineWidth', 3) 
        case 3
            errorbar(0:xlim_max,nanmean(prot_theta_t_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_theta_t_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'r', 'LineWidth', 3) 
        case 4
            errorbar(0:xlim_max,nanmean(prot_theta_t_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_theta_t_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'm', 'LineWidth', 3) 
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Frames (/3.2 ms)'), ylabel('\Delta\theta (\circ)'), title('\Delta\theta_T_o_p', 'FontSize', 15), ylim([-4 8])
plot(0:xlim_max,zeros(1,xlim_max+1),'k--', 'LineWidth', 2)

subplot(224),hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(0:xlim_max,nanmean(prot_theta_f_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_theta_f_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'b', 'LineWidth', 3) 
        case 2
            errorbar(0:xlim_max,nanmean(prot_theta_f_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_theta_f_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'c', 'LineWidth', 3) 
        case 3
            errorbar(0:xlim_max,nanmean(prot_theta_f_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_theta_f_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'r', 'LineWidth', 3) 
        case 4
            errorbar(0:xlim_max,nanmean(prot_theta_f_selected{i}(ts_ind{i},1:xlim_max+1)),nanstd(prot_theta_f_selected{i}(ts_ind{i},1:xlim_max+1)), 'Color', 'm', 'LineWidth', 3) 
    end
end
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; xlabel('Frames (/3.2 ms)'), ylabel('\Delta\theta (\circ)'), title('\Delta\theta_F_r_o_n_t', 'FontSize', 15), ylim([-3 3])
plot(0:xlim_max,zeros(1,xlim_max+1),'k--', 'LineWidth', 2)

legend({'Close Down','Far Down','Close Up','Far Up'})
%% individual plot

i = 1;

figure, 
switch i
    case 1        
        suptitle({[mouseName, ' ', sessionName, ' Close Down'];' ';' '}, 'FontSize', 20)
    case 2
        suptitle({[mouseName, ' ', sessionName, ' Far Down'];' ';' '}, 'FontSize', 20)
    case 3
        suptitle({[mouseName, ' ', sessionName, ' Close Up'];' ';' '}, 'FontSize', 20)
    case 4
        suptitle({[mouseName, ' ', sessionName, ' Far Up'];' ';' '}, 'FontSize', 20)
end

subplot(211),hold all
for j = 1 : length(ts_ind{i})
    if ~isnan(sum(prot_kappa_f_selected{i}(j,1:xlim_max+1)))
        plot(0:xlim_max, prot_kappa_f_selected{i}(ts_ind{i}(j),1:xlim_max+1),'Color',[0.7 0.7 0.7])
    end
end
plot(0:xlim_max, nanmean(prot_kappa_f_selected{i}(ts_ind{i},1:xlim_max+1)),'r-','LineWidth',3)
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; title('\Delta\kappa_F_r_o_n_t'), , ylabel('\Delta\kappa (/mm)'), ylim([-0.002 0.001])
plot(0:xlim_max,zeros(1,xlim_max+1),'k--', 'LineWidth', 3)

subplot(212),hold all
for j = 1 : length(ts_ind{i})
    if ~isnan(sum(prot_theta_f_selected{i}(j,1:xlim_max+1)))
        plot(0:xlim_max, prot_theta_f_selected{i}(ts_ind{i}(j),1:xlim_max+1),'Color',[0.7 0.7 0.7])
    end
end
plot(0:xlim_max, nanmean(prot_theta_f_selected{i}(ts_ind{i},1:xlim_max+1)),'r-','LineWidth',3)
ax = gca; ax.LineWidth = 3; ax.FontWeight = 'bold'; ax.FontSize = 15; title('\Delta\theta_F_r_o_n_t'), xlabel('Frames (/3.2 ms)'), ylabel('\Delta\theta (\circ)'), ylim([-4 6])
plot(0:xlim_max,zeros(1,xlim_max+1),'k--', 'LineWidth', 3)

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



