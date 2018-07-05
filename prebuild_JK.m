%% basic information
% mice = {'JK025','JK027','JK030','JK036','JK037','JK038','JK039','JK041'};
mice = {'JK025'};

videoloc = 'J:\WhiskerVideo\';
if strcmp(videoloc(end),filesep)
    whisker_d = videoloc;
else
    whisker_d = ([videoloc filesep]);
end
behavior_base_dir = 'J:\SoloData\';

ppm = 17.81/2;
            % 'pxPerMm': 17.81002608 for telecentric lens
% comment out when doing for all of the sessions in the mouse directory

rInMm = 3; % mm from the mask along the whisker to calculate delta kappa
%%
%%
%% re-do these
% sessions = {[14,15]};  % for JK039. Did not check JK041 S06~S31 
%%
%%
%%
% sessions = {[4,19,22],[3,16,17],[3,21,22],[1,17,18,91],[7],[2],[22:25],[3]};  
sessions = {[4]};
sessions_pre = {[],[],[],[],[],[],[],[]};
sessions_piezo = {[],[],[],[],[],[],[],[]};
sessions_spont = {[],[],[],[],[],[],[],[]};

all_session = 0; % 1 if using all sessions, 0 if using selected sessions

DoFollicle = 0;
DoRemeasure = 0;
DoWTandWST = 0;
DoWL = 1;

%% Define follicle points and masks
% saves follicle_n_mask.mat file consists of variables 'maskx','masky','width', 'height', and 'follicle_first'

if DoFollicle
    if all_session == 1
        %% use this code when doing for all of the sessions in the mouse directory 
        for i = 1 : size(mice,2)
            cd(whisker_d)
            sn = dir([mice{i},'S*']);
            sn_pre = dir([mice{i},'pre*']);
            sn_piezo = dir([mice{i},'piezo*']);
            sn_spont = dir([mice{i},'spont*']);
            if ~isempty(sn)
                for j = 1 : length(sn)
                    if sn(j).isdir
                        [mouseName, sessionName] = strtok(sn(j).name,'S');
                        if ~isempty(sessionName)
                            follicle_n_mask(mouseName,sessionName,videoloc,'skip')
                        end
                    end
                    close all
                end
            end
            if ~isempty(sn_pre)
                for j = 1 : length(sn_pre)
                    if sn_pre(j).isdir
                        [mouseName, sessionName] = strtok(sn_pre(j).name,'pre');            
                        follicle_n_mask(mouseName,sessionName,videoloc,'skip')                                        
                    end
                    close all
                end
            end
            if ~isempty(sn_piezo)
                for j = 1 : length(sn_piezo)
                    if sn_piezo(j).isdir
                        [mouseName, sessionName] = strtok(sn_piezo(j).name,'piezo');            
                        follicle_n_mask(mouseName,sessionName,videoloc,'skip')                                        
                    end
                    close all
                end
            end  
            if ~isempty(sn_spont)
                for j = 1 : length(sn_spont)
                    if sn_spont(j).isdir
                        [mouseName, sessionName] = strtok(sn_spont(j).name,'spont');            
                        follicle_n_mask(mouseName,sessionName,videoloc,'skip')                                        
                    end
                    close all
                end
            end  
        end
    else
        %% use this code when doing for selected sessions in each mouse directory
        for i = 1 : size(mice,2)
            cd(whisker_d)
            if ~isempty(sessions{i})
                for j = 1 : length(sessions{i})
                    mouseName = mice{i};
                    sessionName = sprintf('S%02d',sessions{i}(j));
                    follicle_n_mask(mouseName,sessionName,videoloc,'skip')
                    close all
                end
            end
            if ~isempty(sessions_pre{i})
                for j = 1 : length(sessions_pre{i})
                    mouseName = mice{i};
                    sessionName = sprintf('pre%d',sessions_pre{i}(j));
                    follicle_n_mask(mouseName,sessionName,videoloc,'skip')
                    close all
                end
            end
            if ~isempty(sessions_piezo{i})
                for j = 1 : length(sessions_piezo{i})
                    mouseName = mice{i};
                    sessionName = sprintf('pre%d',sessions_piezo{i}(j));
                    follicle_n_mask(mouseName,sessionName,videoloc,'skip')
                    close all
                end
            end
            if ~isempty(sessions_spont{i})
                for j = 1 : length(sessions_spont{i})
                    mouseName = mice{i};
                    sessionName = sprintf('pre%d',sessions_spont{i}(j));
                    follicle_n_mask(mouseName,sessionName,videoloc,'skip')
                    close all
                end
            end            
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

if DoRemeasure
    if all_session == 1
        %% use this code when doing for all of the sessions in the mouse directory 
        for i = 1 : size(mice,2)
            cd(whisker_d)
            sn = dir([mice{i},'S*']);
            sn_pre = dir([mice{i}, 'pre*']);
            sn_piezo = dir([mice{i},'piezo*']);
            sn_spont = dir([mice{i},'spont*']);
            if ~isempty(sn)
                for j = 1 : length(sn)
                    if sn(j).isdir
                        [mouseName, sessionName] = strtok(sn(j).name,'S');        
                        if ~isempty(sessionName) 
                            postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
                        end
                    end
                end
            end
            if ~isempty(sn_pre)
                for j = 1 : length(sn_pre)
                    if sn_pre(j).isdir
                        [mouseName, sessionName] = strtok(sn_pre(j).name,'pre');
                        postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
                    end
                end
            end
            if ~isempty(sn_piezo)
                for j = 1 : length(sn_piezo)
                    if sn_piezo(j).isdir
                        [mouseName, sessionName] = strtok(sn_piezo(j).name,'pre');
                        postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
                    end
                end
            end
            if ~isempty(sn_spont)
                for j = 1 : length(sn_spont)
                    if sn_spont(j).isdir
                        [mouseName, sessionName] = strtok(sn_spont(j).name,'pre');
                        postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
                    end
                end
            end            
        end
    else
        %% use this code when doing for selected sessions in each mouse directory
        for i = 1 : size(mice,2)
            cd(whisker_d)
            if ~isempty(sessions{i})
                for j = 1 : length(sessions{i})
                    mouseName = mice{i};
                    sessionName = sprintf('S%02d',sessions{i}(j));
                    postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
                end
            end
            if ~isempty(sessions_pre{i})
                for j = 1 : length(sessions_pre{i})
                    mouseName = mice{i};
                    sessionName = sprintf('pre%d',sessions_pre{i}(j));
                    postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
                end
            end
            if ~isempty(sessions_piezo{i})
                for j = 1 : length(sessions_piezo{i})
                    mouseName = mice{i};
                    sessionName = sprintf('pre%d',sessions_piezo{i}(j));
                    postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
                end
            end
            if ~isempty(sessions_spont{i})
                for j = 1 : length(sessions_spont{i})
                    mouseName = mice{i};
                    sessionName = sprintf('pre%d',sessions_spont{i}(j));
                    postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
                end
            end
        end
    end
end
%% build WT_2pad, WST_2pad, and WL_2pad
% build WL_2pad after touch plane

if DoWTandWST
    cd(whisker_d)
    if all_session == 1
        for mi = 1 : size(mice,2) % mouse index
            cd(whisker_d)
            sn = dir([whisker_d, mice{mi},'S*']);
            for si = 1 : length(sn)
                if sn(si).isdir
                    [mouseName, sessionName] = strtok(sn(si).name,'S');
                    behavior_d = [behavior_base_dir mouseName '\'];
                    if ~isempty(sessionName)
                        if exist('b','var')
                            if strcmp(b{1}.mouseName, mouseName)
                                disp('using the same behavior file')
                            else
                                disp('loading a new behavior file')
                                load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                            end
                        else
                            load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                        end
                        b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        b_session = b{b_ind};

                        buildWTandWST(mouseName, sessionName, whisker_d, b_session, ppm)
                    end
                end
            end

            cd(whisker_d)
            sn_pre = dir([mice{mi},'pre*']);
            for si = 1 : length(sn_pre)
                cd(whisker_d)
                if sn_pre(si).isdir
                    [mouseName, sessionName] = strtok(sn_pre(si).name,'pre');
                    behavior_d = [behavior_base_dir mouseName '\'];
                    if exist('b','var')
                        if strcmp(b{1}.mouseName, mouseName)
                            disp('using the same behavior file')
                        else
                            disp('loading a new behavior file')
                            load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                        end
                    else
                        load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                    end
                    b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                    b_session = b{b_ind};

                    buildWTandWST(mouseName, sessionName, whisker_d, b_session, ppm)
                end
            end
            
            cd(whisker_d)
            sn_piezo = dir([mice{i},'piezo*']);
            for si = 1 : length(sn_piezo)
                cd(whisker_d)
                if sn_piezo(si).isdir
                    [mouseName, sessionName] = strtok(sn_piezo(si).name,'piezo');
                    buildWTandWST(mouseName, sessionName, whisker_d, [], ppm)
                end
            end
            
            cd(whisker_d)
            sn_spont = dir([mice{i},'spont*']);
            for si = 1 : length(sn_spont)
                cd(whisker_d)
                if sn_spont(si).isdir
                    [mouseName, sessionName] = strtok(sn_spont(si).name,'spont');
                    buildWTandWST(mouseName, sessionName, whisker_d, [], ppm)
                end
            end

        end
    else
        for mi = 1 : size(mice,2) % mouse index            
            mouseName = mice{mi};
            if ~isempty(sessions{mi}) 
                for j = 1 : length(sessions{mi})  
                    cd(whisker_d)
                    sessionName = sprintf('S%02d',sessions{mi}(j));
                    if exist([mouseName, sessionName],'dir')
                        behavior_d = [behavior_base_dir mouseName '\'];

                        if exist('b','var')
                            if strcmp(b{1}.mouseName, mouseName)
                                disp('using the same behavior file')
                            else
                                disp('loading a new behavior file')
                                load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                            end
                        else
                            load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                        end
                        if strcmp(sessionName, 'S91')
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S01'), b));
                        else
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        end
                        b_session = b{b_ind};
                        buildWTandWST(mouseName, sessionName, whisker_d, b_session, ppm)
                    end
                end
            end

            if ~isempty(sessions_pre{mi})
                for j = 1 : length(sessions_pre{mi})
                    sessionName = sprintf('pre%d',sessions_pre{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        behavior_d = [behavior_base_dir mouseName '\'];

                        if exist('b','var')
                            if strcmp(b{1}.mouseName, mouseName)
                                disp('using the same behavior file')
                            else
                                disp('loading a new behavior file')
                                load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                            end
                        else
                            load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                        end
                        b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        b_session = b{b_ind};

                        buildWTandWST(mouseName, sessionName, whisker_d, b_session, ppm)
                    end
                end
            end
            
            if ~isempty(sessions_piezo{mi})
                for j = 1 : length(sessions_piezo{mi})
                    sessionName = sprintf('piezo%d',sessions_piezo{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        buildWTandWST(mouseName, sessionName, whisker_d, [], ppm)
                    end
                end
            end
            
            if ~isempty(sessions_spont{mi})
                for j = 1 : length(sessions_spont{mi})
                    sessionName = sprintf('spont%d',sessions_spont{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        buildWTandWST(mouseName, sessionName, whisker_d, [], ppm)
                    end
                end
            end

        end
    end        
end
%% Perfrom touch_hyperplane 
% it includes frame-by-frame estimation of corresponding motor position, based on _WST files
% touch_hyperplane

%% Build WL (Finally)
% it includes touch frame calculation

if DoWL
    cd(whisker_d)
    if all_session == 1
        for mi = 1 : size(mice,2) % mouse index
            cd(whisker_d)
            sn = dir([whisker_d, mice{mi},'S*']);
            for si = 1 : length(sn)
                if sn(si).isdir
                    [mouseName, sessionName] = strtok(sn(si).name,'S');
                    behavior_d = [behavior_base_dir mouseName '\'];
                    if ~isempty(sessionName)
                        if exist('b','var')
                            if strcmp(b{1}.mouseName, mouseName)
                                disp('using the same behavior file')
                            else
                                disp('loading a new behavior file')
                                load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                            end
                        else
                            load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                        end
                        if strcmp(sessionName, 'S91')
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S01'), b));
                        else
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        end
                        b_session = b{b_ind};
                        wd = [whisker_d, mouseName, sessionName];
                        buildWL(wd, b_session, rInMm)
                    end
                end
            end

            cd(whisker_d)
            sn_pre = dir([mice{mi},'pre*']);
            for si = 1 : length(sn_pre)
                cd(whisker_d)
                if sn_pre(si).isdir
                    [mouseName, sessionName] = strtok(sn_pre(si).name,'pre');
                    behavior_d = [behavior_base_dir mouseName '\'];
                    if exist('b','var')
                        if strcmp(b{1}.mouseName, mouseName)
                            disp('using the same behavior file')
                        else
                            disp('loading a new behavior file')
                            load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                        end
                    else
                        load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                    end
                    b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                    b_session = b{b_ind};
                    wd = [whisker_d, mouseName, sessionName];
                    buildWL(wd, b_session, rInMm)
                end
            end

            cd(whisker_d)
            sn_piezo = dir([mice{mi},'piezo*']);
            for si = 1 : length(sn_piezo)
                cd(whisker_d)
                if sn_piezo(si).isdir
                    [mouseName, sessionName] = strtok(sn_piezo(si).name,'piezo');
                    wd = [whisker_d, mouseName, sessionName];
                    buildWL(wd, [], rInMm)
                end
            end
            
            cd(whisker_d)
            sn_spont= dir([mice{mi},'spont*']);
            for si = 1 : length(sn_spont)
                cd(whisker_d)
                if sn_spont(si).isdir
                    [mouseName, sessionName] = strtok(sn_spont(si).name,'piezo');
                    wd = [whisker_d, mouseName, sessionName];
                    buildWL(wd, [], rInMm)
                end
            end
        end
    else
        for mi = 1 : size(mice,2) % mouse index            
            mouseName = mice{mi};
            if ~isempty(sessions{mi}) 
                for j = 1 : length(sessions{mi})  
                    cd(whisker_d)
                    sessionName = sprintf('S%02d',sessions{mi}(j));
                    if exist([mouseName, sessionName],'dir')
                        behavior_d = [behavior_base_dir mouseName '\'];

                        if exist('b','var')
                            if strcmp(b{1}.mouseName, mouseName)
                                disp('using the same behavior file')
                            else
                                disp('loading a new behavior file')
                                load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                            end
                        else
                            load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                        end
                        if strcmp(sessionName, 'S91')
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S01'), b));
                        else
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        end
                        b_session = b{b_ind};
                        wd = [whisker_d, mouseName, sessionName];
                        buildWL(wd, b_session, rInMm)
                    end
                end
            end

            if ~isempty(sessions_pre{mi})
                for j = 1 : length(sessions_pre{mi})
                    sessionName = sprintf('pre%d',sessions_pre{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        behavior_d = [behavior_base_dir mouseName '\'];

                        if exist('b','var')
                            if strcmp(b{1}.mouseName, mouseName)
                                disp('using the same behavior file')
                            else
                                disp('loading a new behavior file')
                                load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                            end
                        else
                            load([behavior_d 'behavior_', mouseName,'.mat']) % loading b of the mouse (all the sessions)
                        end
                        b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        b_session = b{b_ind};
                        wd = [whisker_d, mouseName, sessionName];
                        buildWL(wd, b_session, rInMm)
                    end
                end
            end
            
            if ~isempty(sessions_piezo{mi})
                for j = 1 : length(sessions_piezo{mi})
                    sessionName = sprintf('piezo%d',sessions_piezo{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        wd = [whisker_d, mouseName, sessionName];
                        buildWL(wd, [], rInMm)
                    end
                end
            end
            
            if ~isempty(sessions_spont{mi})
                for j = 1 : length(sessions_spont{mi})
                    sessionName = sprintf('spont%d',sessions_spont{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        wd = [whisker_d, mouseName, sessionName];
                        buildWL(wd, [], rInMm)
                    end
                end
            end            
        end
    end
end





