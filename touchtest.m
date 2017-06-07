% For 2pad (2-port angle distance), WT, WST, and WL is generated during prebuild_JK.m

%% Setup whisker array builder 
behavior_base_dir = 'Z:\Data\2p\soloData\';
whisker_base_dir = 'Z:\Data\Video\JK\';

mice = {'AH0648','AH0650','AH0651','AH0652','AH0653'};

%%%%%%%%%%%%%%%%%%%%%% manual selection
steps = {[10:70],[20:80],[140:200],[140:200]};
%%%%%%%%%%%%%%%%%%%%%% try as short as possible to reduce time next step
mouseName = 'AH0653';
sessionNum = [1,2,4:8,11,14:18];
% sessionNum =[19:21];
trial_types = {'rc', 'rf', 'lc', 'lf'};
% trial_types = {'rn', 'ln'};
%%
for sessionInd = 1 : length(sessionNum)
% for sessionInd = 1:3    
    sessionName = sprintf('S%02d',sessionNum(sessionInd));

    behavior_d = [behavior_base_dir mouseName '\'];
    whisker_d = [whisker_base_dir mouseName sessionName '\'];

    % if exist([whisker_d, 'touch_hp.mat'],'file')
    %     error('touch_hp.mat exists.')    
    % end

    if exist('b','var')
        if strcmp(b{1}.mouseName, mouseName)
            disp('using the same behavior file')
        else
            disp('loading a new behavior file')
            load([behavior_d 'behavior.mat']) % loading b of the mouse (all the sessions)
        end
    else
        load([behavior_d 'behavior.mat']) % loading b of the mouse (all the sessions)
    end
    b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
    b_session = b{b_ind};

    load_fn = [mouseName sessionName '_post.mat'];
    load([whisker_d load_fn]); % loading errorlist

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
    touch_hp = cell(1,length(trial_types)); % touch hyperplanes
    hp_peaks = cell(1,length(trial_types)); % touch hyperplane peak points. 2 points for each hyperplane

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
    figure,
    for trial_type_num = 1 : length(trial_types)
        intersect_3d_total = [];
        wl = wl_array{trial_type_num};        
        for tnum = 1 : length(wl.trials)
            try        
                top_ind = find(~isnan(wl.trials{tnum}.intersect_coord(:,1)));
                front_ind = find(~isnan(wl.trials{tnum}.intersect_coord(:,2)));
                intersect_ind = intersect(wl.trials{tnum}.pole_available_timepoints,intersect(top_ind,front_ind));
                intersect_3d_total = [intersect_3d_total; wl.trials{tnum}.intersect_coord(intersect_ind,1), wl.trials{tnum}.intersect_coord(intersect_ind,2), ones(length(intersect_ind),1)*wl.trials{tnum}.pole_pos];
            catch
                fprintf('Skipping trial #%d because of index problems \n',tnum);        
            end
        end
        intersect_3d_total = unique(round(intersect_3d_total,2),'rows');
        subplot(2,2,trial_type_num), plot3(intersect_3d(:,1), intersect_3d(:,2), intersect_3d(:,3), 'k.', 'MarkerSize', 0.1)
        title(wl.trials{1}.trial_type), xlabel('Top-view intersection coord'), ylabel('Front-view intersection coord'), zlabel('Pole position')
        drawnow;
        waitforbuttonpress;
    end
end