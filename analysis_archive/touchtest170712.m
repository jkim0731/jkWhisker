% % %% 170712 finding threshold for delta Kappa 
% % close all
% % clear
% % base_angle = 21;
% % behavior_base_dir = 'Z:\Data\2p\soloData\';
% % whisker_base_dir = 'Z:\Data\Video\JK\';
% % trial_types = {'rc', 'rf', 'lc', 'lf'};
% % mice = {'AH0648','AH0650','AH0651','AH0652','AH0653'};
% % % sessionNum = {[1:4,6:15],[1,2,4:6,8:10],[2,4:6,8:19],[6,8:13],[2,4:8,13:19]};
% % 
% % sessions{1} = {[4,6,7], [13:15]}; % 0648
% % sessions{2} = {[4:6], [10]}; % 0650
% % sessions{3} = {[4:6], [6,8,10]}; % 0651
% % sessions{4} = {[6,8], [10,11,13]}; % 0652
% % sessions{5} = {[4:6], [16:18]}; % 0653
% % 
% % mouseName = mice{1};
% % sessionName = 'S13';
% % 
% % whisker_d = [whisker_base_dir mouseName sessionName '\'];
% % if ~exist('b','var') || ~iscell(b) || ~isfield(b{1},'mouseName') || ~strcmp(b{1}.mouseName,mouseName)
% %     load([behavior_base_dir mouseName filesep 'behavior.mat']) % load b
% % end
% % 
% % b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
% % b_session = b{b_ind};
% % tt_ind = cell(1,length(trial_types));
% % wl = cell(1,length(trial_types));
% % 
% % for ti = 1 : length(trial_types)    
% %     tt_ind{ti} = find(cellfun(@(x) strcmp(x.trialType,trial_types{ti}),b_session.trials));
% %     temp_files = cell(length(tt_ind{ti}),1);
% %     for j = 1 : length(tt_ind{ti})
% %         temp_files{j} = num2str(tt_ind{ti}(j));
% %     end
% %     wl{ti} = Whisker.WhiskerTrialLiteArray_2pad(whisker_d,'include_files',temp_files);     
% % end
% %% assign deltakappa and theta. deltakappa is a true deltakappa, calculated by subtracting the kappa value of one frame before the first touch frame in each chunk
% 
% deltakappa_t = cell(1,length(trial_types));
% deltakappa_f = cell(1,length(trial_types));
% theta_t = cell(1,length(trial_types));
% theta_f = cell(1,length(trial_types));
% deltatheta_t = cell(1,length(trial_types));
% deltatheta_f = cell(1,length(trial_types));
% for ti = 1 : length(trial_types)
%     for j = 1 : length(wl{ti}.trials) 
%         wt = wl{ti}.trials{j};
%         trial_temp_ind = find(cellfun(@(x) x.trialNum == wt.trialNum,b_session.trials));
%         trial_temp = b_session.trials{trial_temp_ind};
% 
%         temp_touch = wt.th_touch_frames;
%         if ~isempty(temp_touch)
%             touch_diff_inds = [0;find(diff(temp_touch) - 1);length(temp_touch)];
%             for cind = 1 : length(touch_diff_inds)-1
%                 if length(wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * ...
%                         length(wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) ...
%                         > 0 % non of these are empty
%                     deltakappa_t{ti}{end+1} = wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1)-1);
%                     deltakappa_f{ti}{end+1} = wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1)-1);
%                     theta_t{ti}{end+1} = wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
%                     theta_f{ti}{end+1} = wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
%                     deltatheta_t{ti}{end+1} = wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1));
%                     deltatheta_f{ti}{end+1} = wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1));
%                 end
%             end
%         end
%     end
% end
% 
% %% collapsing deltakappa and theta values into one matrix for easier averaging
% 
% dkt = cell(1,length(trial_types)); %delta kappa _t
% dkf = cell(1,length(trial_types)); % delta kappa _f
% tht = cell(1,length(trial_types)); % theta _t
% thf = cell(1,length(trial_types)); % theta _f
% dtt = cell(1,length(trial_types)); % delta theta _t
% dtf = cell(1,length(trial_types)); % delta theta _f
% 
% for ti = 1 : length(trial_types)
%     dkt{ti} = zeros(length(deltakappa_t{tind}),max(cellfun(@(x) length(x),deltakappa_t{tind})));
%     dkf{ti} = zeros(size(dkt{ti}));
%     tht{ti} = zeros(size(dkt{ti}));
%     thf{ti} = zeros(size(dkt{ti}));
%     dtt{ti} = zeros(size(dkt{ti}));
%     dtf{ti} = zeros(size(dkt{ti}));
% end
% for tind = 1 : length(trial_types)
%     for i = 1 : length(deltakappa_t{tind})
%         dkt{tind}(i,1:length(deltakappa_t{tind}{i})) = deltakappa_t{tind}{i}; if length(deltakappa_t{tind}{i}) < size(dkt,2), dkt{tind}(i,length(deltakappa_t{tind}{i})+1:end) = NaN; end;
%         dkf{tind}(i,1:length(deltakappa_f{tind}{i})) = deltakappa_f{tind}{i}; if length(deltakappa_f{tind}{i}) < size(dkt,2), dkf{tind}(i,length(deltakappa_f{tind}{i})+1:end) = NaN; end;
%         tht{tind}(i,1:length(theta_t{tind}{i})) = theta_t{tind}{i}; if length(theta_t{tind}{i}) < size(dkt,2), tht{tind}(i,length(theta_t{tind}{i})+1:end) = NaN; end;
%         thf{tind}(i,1:length(theta_f{tind}{i})) = theta_f{tind}{i}; if length(theta_f{tind}{i}) < size(dkt,2), thf{tind}(i,length(theta_f{tind}{i})+1:end) = NaN; end;
%         dtt{tind}(i,1:length(deltatheta_t{tind}{i})) = deltatheta_t{tind}{i}; if length(deltatheta_t{tind}{i}) < size(dkt,2), dtt{tind}(i,length(deltatheta_t{tind}{i})+1:end) = NaN; end;
%         dtf{tind}(i,1:length(deltatheta_f{tind}{i})) = deltatheta_f{tind}{i}; if length(deltatheta_f{tind}{i}) < size(dkt,2), dtf{tind}(i,length(deltatheta_f{tind}{i})+1:end) = NaN; end;        
%     end
% end
% 
% %% averaging deltakappa and theta values altogether 
% for i = 1 : length(trial_types)
%     figure, subplot(221), errorbar(1:size(dkt{i},2),nanmean(dkt{i}),nanstd(dkt{i})), title('delta kappa top'), xlim([0 20])
%     subplot(222), errorbar(1:size(dkf{i},2),nanmean(dkf{i}),nanstd(dkf{i})), title('delta kappa front'), xlim([0 20])
%     subplot(223), errorbar(1:size(dtt{i},2),nanmean(dtt{i}),nanstd(dtt{i})), title('delta theta top'), xlim([0 20])
%     subplot(224), errorbar(1:size(dtf{i},2),nanmean(dtf{i}),nanstd(dtf{i})), title('delta theta front'), xlim([0 20])
% end
% 
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
                    protraction_deltakappa_t{ti}{end+1} = wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1)-1);
                    protraction_deltakappa_f{ti}{end+1} = wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) - wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1)-1);
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
figure, hold all
subplot(221), title('delta kappa top'), hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(1:size(protraction_dkt{i},2),nanmean(protraction_dkt{i}),nanstd(protraction_dkt{i}), 'Color', 'b') 
        case 2
            errorbar(1:size(protraction_dkt{i},2),nanmean(protraction_dkt{i}),nanstd(protraction_dkt{i}), 'Color', 'c') 
        case 3
            errorbar(1:size(protraction_dkt{i},2),nanmean(protraction_dkt{i}),nanstd(protraction_dkt{i}), 'Color', 'r') 
        case 4
            errorbar(1:size(protraction_dkt{i},2),nanmean(protraction_dkt{i}),nanstd(protraction_dkt{i}), 'Color', 'm') 
    end
end
subplot(222), title('delta kappa front'), xlim([0 20]), hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(1:size(protraction_dkf{i},2),nanmean(protraction_dkf{i}),nanstd(protraction_dkf{i}), 'Color', 'b') 
        case 2
            errorbar(1:size(protraction_dkf{i},2),nanmean(protraction_dkf{i}),nanstd(protraction_dkf{i}), 'Color', 'c') 
        case 3
            errorbar(1:size(protraction_dkf{i},2),nanmean(protraction_dkf{i}),nanstd(protraction_dkf{i}), 'Color', 'r') 
        case 4
            errorbar(1:size(protraction_dkf{i},2),nanmean(protraction_dkf{i}),nanstd(protraction_dkf{i}), 'Color', 'm') 
    end
end
subplot(223), title('delta theta top'), xlim([0 20]), hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(1:size(protraction_dtt{i},2),nanmean(protraction_dtt{i}),nanstd(protraction_dtt{i}), 'Color', 'b') 
        case 2
            errorbar(1:size(protraction_dtt{i},2),nanmean(protraction_dtt{i}),nanstd(protraction_dtt{i}), 'Color', 'c') 
        case 3
            errorbar(1:size(protraction_dtt{i},2),nanmean(protraction_dtt{i}),nanstd(protraction_dtt{i}), 'Color', 'r') 
        case 4
            errorbar(1:size(protraction_dtt{i},2),nanmean(protraction_dtt{i}),nanstd(protraction_dtt{i}), 'Color', 'm') 
    end
end
subplot(224), title('delta theta front'), xlim([0 20]), hold all
for i = 1 : length(trial_types)
    switch i
        case 1
            errorbar(1:size(protraction_dtf{i},2),nanmean(protraction_dtf{i}),nanstd(protraction_dtf{i}), 'Color', 'b') 
        case 2
            errorbar(1:size(protraction_dtf{i},2),nanmean(protraction_dtf{i}),nanstd(protraction_dtf{i}), 'Color', 'c') 
        case 3
            errorbar(1:size(protraction_dtf{i},2),nanmean(protraction_dtf{i}),nanstd(protraction_dtf{i}), 'Color', 'r') 
        case 4
            errorbar(1:size(protraction_dtf{i},2),nanmean(protraction_dtf{i}),nanstd(protraction_dtf{i}), 'Color', 'm') 
    end
end



% %% collapse and average only those from retraction touch
% 
%     
% %%
% % Whisker.viewdouble_WhiskerTrialLiteArray(wl{1},[0 1])
% 
