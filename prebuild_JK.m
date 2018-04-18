
%% basic information
mice = {'JK025'};
videoloc = 'Y:\Whiskernas\JK_temp\whisker\tracked\';
if strcmp(videoloc(end),filesep)
    whisker_d = videoloc;
else
    whisker_d = ([videoloc filesep]);
end
behavior_base_dir = 'Y:\Whiskernas\JK_temp\SoloData\';

ppm = 17.81/2;
            % 'pxPerMm': 17.81002608 for telecentric lens
ppm = ppm / 2; % for binning 2
            % 'pxPerMm': 10.56526073 for microVideo -------------------------------------------------------------------------------------------------------------------------------------lens
% comment out when doing for all of the sessions in the mouse directory
sessions = {[4,19]};  
sessions_pre = {[]};

all_session = 0; % 1 if using all sessions, 0 if using selected sessions
networkfailtime = [];
%% Define follicle points and masks
% saves follicle_n_mask.mat file consists of variables 'maskx','masky','width', 'height', and 'follicle_first'
% if all_session == 1
%     %% use this code when doing for all of the sessions in the mouse directory 
%     for i = 1 : size(mice,2)
%         cd(d)
%         sn = dir([mice{i},'S*']);
%         sn_pre = dir([mice{i},'pre*']); 
%         if ~isempty(sn)
%             for j = 1 : length(sn)
%                 if sn(j).isdir
%                     [mouseName, sessionName] = strtok(sn(j).name,'S');
%                     if ~isempty(sessionName)
%                         follicle_n_mask(mouseName,sessionName,videoloc,'skip')
%                     end
%                 end
%                 close all
%             end
%         end
%         if ~isempty(sn_pre)
%             for j = 1 : length(sn_pre)
%                 if sn_pre(j).isdir
%                     [mouseName, sessionName] = strtok(sn_pre(j).name,'pre');            
%                     follicle_n_mask(mouseName,sessionName,videoloc,'skip')                                        
%                 end
%                 close all
%             end
%         end        
%     end
% else
%     %% use this code when doing for selected sessions in each mouse directory
%     for i = 1 : size(mice,2)
%         cd(d)
%         if ~isempty(sessions{i})
%             for j = 1 : length(sessions{i})
%                 mouseName = mice{i};
%                 sessionName = sprintf('S%02d',sessions{i}(j));
%                 follicle_n_mask(mouseName,sessionName,videoloc,'skip')
%                 close all
%             end
%         end
%         if ~isempty(sessions_pre{i})
%             for j = 1 : length(sessions_pre{i})
%                 mouseName = mice{i};
%                 sessionName = sprintf('pre%d',sessions_pre{i}(j));
%                 follicle_n_mask(mouseName,sessionName,videoloc,'skip')
%                 close all
%             end
%         end
%     end
% end

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

% if all_session == 1
%     %% use this code when doing for all of the sessions in the mouse directory 
%     for i = 1 : size(mice,2)
%         cd(d)
%         sn = dir([mice{i},'S*']);
%         sn_pre = dir([mice{i}, 'pre*']); 
%         if ~isempty(sn)
%             for j = 1 : length(sn)
%                 if sn(j).isdir
%                     [mouseName, sessionName] = strtok(sn(j).name,'S');        
%                     if ~isempty(sessionName) 
%                         postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
%                     end
%                 end
%             end
%         end
%         if ~isempty(sn_pre)
%             for j = 1 : length(sn_pre)
%                 if sn_pre(j).isdir
%                     [mouseName, sessionName] = strtok(sn_pre(j).name,'pre');
%     %                 postmeasurements(mouseName,sessionName,videoloc,ppm)
%                     postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
%                 end
%             end
%         end
%     end
% else
%     %% use this code when doing for selected sessions in each mouse directory
%     for i = 1 : size(mice,2)
%         cd(d)
%         if ~isempty(sessions{i})
%             for j = 1 : length(sessions{i})
%                 mouseName = mice{i};
%                 sessionName = sprintf('S%02d',sessions{i}(j));
% %                 postmeasurements(mouseName,sessionName,videoloc,ppm)
%                 postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
%             end
%         end
%         if ~isempty(sessions_pre{i})
%             for j = 1 : length(sessions_pre{i})
%                 mouseName = mice{i};
%                 sessionName = sprintf('pre%d',sessions_pre{i}(j));
% %                 postmeasurements(mouseName,sessionName,videoloc,ppm)
%                 postmeasurements(mouseName,sessionName,videoloc,ppm,'skip')
%             end
%         end
%     end
% end

%% build WT_2pad, WST_2pad, and WL_2pad
% build WL_2pad after touch plane
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
                    b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
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
    end
end        

%% Perfrom touch_hyperplane 
% it includes frame-by-frame estimation of corresponding motor position, based on _WST files
% touch_hyperplane

%% Build WL (Finally)
% it includes touch frame calculation
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
                    
                    buildWL(whisker_d, b_session)
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
                
                buildWL(whisker_d, b_session)
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
                    b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                    b_session = b{b_ind};

                    buildWL(whisker_d, b_session)
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

                    buildWL(whisker_d, b_session)
                end
            end
        end
    end
end




