
%% basic information
% mice = {'AH0648', 'AH0650', 'AH0651', 'AH0652', 'AH0653'};
mice = {'AH0650'};
videoloc = 'JK';
d = (['Z:\Data\Video\' videoloc filesep]);

ppm = 10.56526073;
            % 'pxPerMm': 17.81002608 for telecentric lens
            % 'pxPerMm': 10.56526073 for microVideo lens
% comment out when doing for all of the sessions in the mouse directory
% sessions = {{'S02'}, {'S02'}, {'S02'}, {'S04'}, {'S03'}}; 
sessions = {{'S19','S20','S21'}}; 
all_session = 1; % 1 if using all sessions, 0 if using selected sessions
%% Define follicle points and masks
% saves follicle_n_mask.mat file consists of variables 'maskx','masky','width', 'height', and 'follicle_first'
if all_session == 1
    %% use this code when doing for all of the sessions in the mouse directory 
    for i = 1 : size(mice,2)
        cd(d)
        sn = dir([mice{i},'*']);
        for j = 1 : length(sn)
    %         if sn(j).isdir
            if sn(j).isdir && ~strcmp(sn(j).name,strcat(mice{i},sessions_exc{i}))
                [mouseName, sessionName] = strtok(sn(j).name,'S');            
                follicle_n_mask(mouseName,sessionName,videoloc)
            end
            close all
        end
    end
else
    %% use this code when doing for selected sessions in each mouse directory
    for i = 1 : size(mice,2)
        cd(d)
        sn = sessions{i};
        for j = 1 : length(sn)
            mouseName = mice{i};
            sessionName = sn{j};
            follicle_n_mask(mouseName,sessionName,videoloc)
            close all
        end
    end
end

%%
%%
%%
%%
%% TODO: MAKE MASK CHECKING PROCEDURE
%%
%%
%%
%%

%% Re-measure .measurements file before building whisker arrays 
% saves _post.mat file consists of variables 'maskx','masky', 'includef', 'errorlist', 'width', and 'height'
% postmeasurements uses jk_measurements_dir()

if all_session == 1
    %% use this code when doing for all of the sessions in the mouse directory 
    for i = 1 : size(mice,2)
        cd(d)
        sn = dir([mice{i},'*']);
        for j = 1 : length(sn)
    %         if sn(j).isdir
            if sn(j).isdir && ~strcmp(sn(j).name,strcat(mice{i},sessions_exc{i}))
                [mouseName, sessionName] = strtok(sn(j).name,'S');
                postmeasurements(mouseName,sessionName,videoloc)
            end
        end
    end
else
    %% use this code when doing for selected sessions in each mouse directory
    for i = 1 : size(mice,2)
        sn = sessions{i};
        for j = 1 : length(sn)
            mouseName = mice{i};
            sessionName = sn{j};
            postmeasurements(mouseName,sessionName,videoloc)
        end
    end
end

%% build WT, WST, and WL
behavior_base_dir = 'Z:\Data\2p\soloData\';
if all_session == 1
    for mi = 1 : size(mice,2) % mouse index
        whisker_d = [d, mice{mi}];
        sn = dir([mice{mi},'*']);
        for si = 1 : length(sn)
            if sn(si).isdir && ~strcmp(sn(si).name,strcat(mice{mi},sessions_exc{mi}))
                [mouseName, sessionName] = strtok(sn(si).name,'S');      
                behavior_d = [behavior_base_dir mouseName '\'];
                whisker_d = [d, mice{mi}, sn{si}];
                load([whisker_d, filesep, mice{mi}, sn{si}, '_post.mat'])
            
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
            
                Whisker.makeAllDirectory_WhiskerTrial(whisker_d,[0 1],'mask', {[maskx(1,:);masky(1,:)],[maskx(2,:);masky(2,:)]},...
                    'trial_nums',trialNums,'include_files',includef,...
                    'barRadius',15.3,'faceSideInImage', 'bottom', 'framePeriodInSec',.0032,...
                    'imagePixelDimsXY',[width height],'pxPerMm',ppm,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward');
                Whisker.makeAllDirectory_WhiskerSignalTrial_2pad(whisker_d,'include_files',includef,'polyRoiInPix',[20 80]);
                Whisker.makeAllDirectory_WhiskerTrialLiteI(whisker_d,'include_files',includef,'r_in_mm',3,'calc_forces',false,'behavior',b_session);
            end
        end
    end
else
    for mi = 1 : size(mice,2) % mouse index
        sn = sessions{mi};
        mouseName = mice{mi};
        for si = 1 : length(sn) % session index
            behavior_d = [behavior_base_dir mouseName '\'];
            whisker_d = [d, mice{mi}, sn{si}];
            load([whisker_d, filesep, mice{mi}, sn{si}, '_post.mat'])

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
            
            Whisker.makeAllDirectory_WhiskerTrial(whisker_d,[0 1],'mask', {[maskx(1,:);masky(1,:)],[maskx(2,:);masky(2,:)]},...
                'trial_nums',trialNums,'include_files',includef,...
                'barRadius',15.3,'faceSideInImage', 'bottom', 'framePeriodInSec',0.003225806451613,...
                'imagePixelDimsXY',[width height],'pxPerMm',ppm,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward');
            Whisker.makeAllDirectory_WhiskerSignalTrial_2pad(whisker_d,'include_files',includef,'polyRoiInPix',[20 80]);
            Whisker.makeAllDirectory_WhiskerTrialLiteI(whisker_d,'include_files',includef,'r_in_mm',3,'calc_forces',false,'behavior',b_session);
        end
    end
end        