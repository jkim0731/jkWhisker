
%% basic information
mice = {'JK025','JK027','JK030','JK036','JK037','JK038','JK039','JK041'};
% mice = {'JK027'};
videoloc = 'Y:\Whiskernas\JK_temp\whisker\tracked\';
if strcmp(videoloc(end),filesep)
    whisker_d = ([videoloc filesep]);
else
    whisker_d = videoloc;
end
behavior_base_dir = 'Y:\Whiskernas\JK_temp\SoloData\';

ppm = 17.81002608;
            % 'pxPerMm': 17.81002608 for telecentric lens
ppm = ppm / 2; % for binning 2
            % 'pxPerMm': 10.56526073 for microVideo -------------------------------------------------------------------------------------------------------------------------------------lens
% comment out when doing for all of the sessions in the mouse directory
sessions = {[8:19,22],[1:22],[1:22],[1:18,21],[1:24],[1:22,24:31],[1:28],[1:19,21:30]};  
sessions_pre = {[1],[1],[1,2],[1],[1,2],[1],[1],[1]};

all_session = 0; % 1 if using all sessions, 0 if using selected sessions

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
cd(videoloc)
if all_session == 1
    for mi = 1 : size(mice,2) % mouse index
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
        
        sn_pre = dir([mice{mi},'pre*']);
        for si = 1 : length(sn_pre)
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
                sessionName = sprintf('S%02d',sessions{mi}(j));
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
        
        if ~isempty(sessions_pre{mi})
            for j = 1 : length(sessions_pre{mi})
                sessionName = sprintf('pre%d',sessions_pre{mi}(j));
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