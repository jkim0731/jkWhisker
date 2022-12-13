%% basic information

%%
%% Error sessions!!

% 2022/12/12 Running all JK025-056 again for WL and WF.
% One thing to change at buildWL_2pad.m
% Might not be a big problem (touch hyperplane adjustment.

% 2018/12/22 Run JK052pre1 W3.

% Run JK027S17 WL_2pad again. error at 2nd running. check with wlarray first.
% Error sessions for length(protraction frames) ~= length(prtraction distance) -> Run WL on these first
%%
%%
%% re-do these
% sessions = {[14,15]};  % for JK039. Did not check JK041 S06~S31 
%%
%%
%%
%%
%%
%%

% mice = {'JK025','JK027','JK030','JK036','JK037','JK038','JK039','JK041','JK052', 'JK053','JK054','JK056', 'JK070', 'JK074', 'JK075', 'JK076'};
mice = {'JK025','JK027','JK030','JK036','JK037','JK038','JK039','JK041','JK052', 'JK053','JK054','JK056'};
% mice = {'JK036','JK037','JK038','JK039','JK041','JK052', 'JK053','JK054','JK056'};

% 

% mice = {'JK027'};
% mice = {'JK052','JK056'};
% mice = {'JK056'};
% mice = {'JK054','JK056', 'JK070', 'JK074', 'JK075', 'JK076'};
% mice = {'JK025','JK027','JK030','JK036','JK037','JK038','JK039','JK041'};
% mice = {'JK038','JK039','JK041'};

% videoloc = 'D:\TPM\JK\tracked\';
videoloc = 'D:\WhiskerVideo\';
if strcmp(videoloc(end),filesep)
    whisker_d = videoloc;
else
    whisker_d = ([videoloc filesep]);
end
% behavior_base_dir = 'D:\TPM\JK\soloData\';
behavior_base_dir = 'D:\SoloData\';

ppm = 17.81;
            % 'pxPerMm': 17.81002608 for telecentric lens
            % /2 for mice <= JK041, because of binning.
% comment out when doing for all of the sessions in the mouse directory
maskmm = 1; % mm from the face to draw the mask
facePosition = 'bottom';
rInMm = 3; % mm from the mask along the whisker to calculate delta kappa
follicleSkip = 'skip'; % 'skip' or 'noskip'
remeasureSkip = 'skip'; % 'skip' or 'noskip'
touchHyperplaneSkip = 'noskip';
videoFreq = 311.24; % frequency of whisker video imaging. If 0, then use timestamp file (calculated from .seq file)
barRadius = 0.3; % in mm

% parameters for refining touch frames
whiskingAmpThreshold = 2.5; % in degrees
stdHistogramThreshold = 1;
distanceHistogramBinInMm = 0.02; %
distanceHistogramBin = round(ppm*distanceHistogramBinInMm*100)/100; % up to 2 significant numbers
touchBoundaryThickness = 0.3; % in mm
touchBoundaryBuffer = 0.1; % in mm
maxPointsNearHyperplane = videoFreq * 15 / 1000; % mean touch duration ~ 15 ms. 
touchKappaSTDthreshold = 2;

% sessions = {[4,5,18,19,22],[2,3,7,10,14],[3,4,20,21,22],[1,2,16,17,18,91],[7],[2],[1,2,21,22:25],[3]};
% sessions = {[3,4,21:23,25,26],[3],[3,4,25,26],[3:5]};
% 
% sessions = {[1:22],[1:22,99],[1:22],[1:20,91],[1:24],[1:31],[1:25],[1:30]};
% sessions_pre = {[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3]};
% 

% sessions = {    [1:19,22],  [1:22,99], [1:7,9:22], [1:18,91],  [1:10,12:24],   [1:22,24:31],   [1:25], [1:19,21:30],   [1:29,94,95],     [1:3,5:21],     [1:26], [1:13]};
% sessions_pre = {[1],        [1],    [1:2],  [1],        [1:2],          [1],            [1],    [1],            [1],        [1],            [1],    [1]};


% sessions_pre = {[1],        [1],    [1:2],  [1],        [1:2],          [1],            [1],    [1]};

% sessions = {[4]};

% sessions_spont = {[],[],[],[],[],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5]};

% 
% sessions = {[],[],[],[],[],[],[],[],[],[]};
% sessions_pre = {[],[],[],[],[],[],[],[],[],[],[],[]};
% sessions_piezo = {[],[],[],[],[],[],[],[],[],[],[],[]};
% sessions_spont = {[],[],[],[],[],[],[],[],[],[],[],[]};

% sessionsDone = {[4,5,18,19,22],[3,4,15,16,17],[3,4,20,21,22],[1,2,16,17,18],[7],[2],[1,2,21,22:25],[3], [3:5,24:26], [3],[3],[3]};
% sessionsTorun = {   [1:19,22],  [1:22,99], [1:7,9:22], [1:18,91],  [1:10,12:24],   [1:22,24:31],   [1:25], [1:19,21:30]};
% sessionsDone = {[3:5]};
% sessionsTorun = {[8]};
% sessions = cell(1,length(sessionsTorun));
% for i = 1 : length(sessions)
%     sessions{i} = setdiff(sessionsTorun{i}, sessionsDone{i});
% end
% sessions = {[1:29,94,95],     [1:3,5:21],     [1:26], [1:13]};
% sessions = {[91],[91,92],[91],[91],[91],[91],[91],[91],[91]};
% sessions = {[91]};
all_session = 1; % 1 if using all sessions, 0 if using selected sessions

DoFollicle = 0;
DoRemeasure = 0;
doWT = 0;
testPoleUp = 0;
doWST = 0;
makeTouchHyperplane = 0;
doWL = 1;
do3D = 0;
doWF = 1;

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
                            follicle_n_mask(mouseName,sessionName,videoloc, ppm, maskmm, facePosition, follicleSkip)                            
                        end
                    end
                    close all
                end
            end
            if ~isempty(sn_pre)
                for j = 1 : length(sn_pre)
                    if sn_pre(j).isdir
                        [mouseName, sessionName] = strtok(sn_pre(j).name,'pre');            
                        follicle_n_mask(mouseName,sessionName,videoloc, ppm, maskmm, facePosition, follicleSkip)                                        
                    end
                    close all
                end
            end
            if ~isempty(sn_piezo)
                for j = 1 : length(sn_piezo)
                    if sn_piezo(j).isdir
                        [mouseName, sessionName] = strtok(sn_piezo(j).name,'piezo');            
                        follicle_n_mask(mouseName,sessionName,videoloc, ppm, maskmm, facePosition, follicleSkip)                                        
                    end
                    close all
                end
            end  
            if ~isempty(sn_spont)
                for j = 1 : length(sn_spont)
                    if sn_spont(j).isdir
                        [mouseName, sessionName] = strtok(sn_spont(j).name,'spont');            
                        follicle_n_mask(mouseName,sessionName,videoloc, ppm, maskmm, facePosition, follicleSkip)                                        
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
                    follicle_n_mask(mouseName,sessionName,videoloc, ppm, maskmm, facePosition, follicleSkip)
                    close all
                end
            end
            if ~isempty(sessions_pre{i})
                for j = 1 : length(sessions_pre{i})
                    mouseName = mice{i};
                    sessionName = sprintf('pre%d',sessions_pre{i}(j));
                    follicle_n_mask(mouseName,sessionName,videoloc, ppm, maskmm, facePosition, follicleSkip)
                    close all
                end
            end
            if ~isempty(sessions_piezo{i})
                for j = 1 : length(sessions_piezo{i})
                    mouseName = mice{i};
                    sessionName = sprintf('piezo%d',sessions_piezo{i}(j));
                    follicle_n_mask(mouseName,sessionName,videoloc, ppm, maskmm, facePosition, follicleSkip)
                    close all
                end
            end
            if ~isempty(sessions_spont{i})
                for j = 1 : length(sessions_spont{i})
                    mouseName = mice{i};
                    sessionName = sprintf('spont%d',sessions_spont{i}(j));
                    follicle_n_mask(mouseName,sessionName,videoloc, ppm, maskmm, facePosition, follicleSkip)
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
                            postmeasurements(mouseName,sessionName,videoloc,ppm,remeasureSkip)
                        end
                    end
                end
            end
            if ~isempty(sn_pre)
                for j = 1 : length(sn_pre)
                    if sn_pre(j).isdir
                        [mouseName, sessionName] = strtok(sn_pre(j).name,'pre');
                        postmeasurements(mouseName,sessionName,videoloc,ppm,remeasureSkip)
                    end
                end
            end
            if ~isempty(sn_piezo)
                for j = 1 : length(sn_piezo)
                    if sn_piezo(j).isdir
                        [mouseName, sessionName] = strtok(sn_piezo(j).name,'pizeo');
                        postmeasurements(mouseName,sessionName,videoloc,ppm,remeasureSkip)
                    end
                end
            end
            if ~isempty(sn_spont)
                for j = 1 : length(sn_spont)
                    if sn_spont(j).isdir
                        [mouseName, sessionName] = strtok(sn_spont(j).name,'spont');
                        postmeasurements(mouseName,sessionName,videoloc,ppm,remeasureSkip)
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
                    postmeasurements(mouseName,sessionName,videoloc,ppm,remeasureSkip)
                end
            end
            if ~isempty(sessions_pre{i})
                for j = 1 : length(sessions_pre{i})
                    mouseName = mice{i};
                    sessionName = sprintf('pre%d',sessions_pre{i}(j));
                    postmeasurements(mouseName,sessionName,videoloc,ppm,remeasureSkip)
                end
            end
            if ~isempty(sessions_piezo{i})
                for j = 1 : length(sessions_piezo{i})
                    mouseName = mice{i};
                    sessionName = sprintf('piezo%d',sessions_piezo{i}(j));
                    postmeasurements(mouseName,sessionName,videoloc,ppm,remeasureSkip)
                end
            end
            if ~isempty(sessions_spont{i})
                for j = 1 : length(sessions_spont{i})
                    mouseName = mice{i};
                    sessionName = sprintf('spont%d',sessions_spont{i}(j));
                    postmeasurements(mouseName,sessionName,videoloc,ppm,remeasureSkip)
                end
            end
        end
    end
end
%% build WT_2pad, WST_2pad, and WL_2pad
% build WL_2pad after touch plane

if doWT
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
                        elseif strcmp(sessionName, 'S94') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S04'), b));
                        elseif strcmp(sessionName, 'S95') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S05'), b));
                        elseif strcmp(sessionName, 'S99') % JK027
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S17'), b)); 
                        else
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        end                        
                        b_session = b{b_ind};

                        buildWT_2pad(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
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

                    buildWT_2pad(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
                end
            end
            
            cd(whisker_d)
            sn_piezo = dir([mice{mi},'piezo*']);
            for si = 1 : length(sn_piezo)
                cd(whisker_d)
                if sn_piezo(si).isdir
                    [mouseName, sessionName] = strtok(sn_piezo(si).name,'piezo');
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
                    if ~isempty(b_ind)
                        b_session = b{b_ind};
                        buildWT_2pad(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
                    else
                        buildWT_2pad(mouseName, sessionName, whisker_d, [], videoFreq, ppm, barRadius)
                    end
                end
            end
            
            cd(whisker_d)
            sn_spont = dir([mice{mi},'spont*']);
            for si = 1 : length(sn_spont)
                cd(whisker_d)
                if sn_spont(si).isdir
                    [mouseName, sessionName] = strtok(sn_spont(si).name,'spont');
                    buildWT_2pad(mouseName, sessionName, whisker_d, [], videoFreq, ppm, barRadius)
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
                        elseif strcmp(sessionName, 'S94') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S04'), b));
                        elseif strcmp(sessionName, 'S95') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S05'), b));
                        elseif strcmp(sessionName, 'S99') % JK027
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S17'), b)); 
                        else
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        end
                        b_session = b{b_ind};
                        buildWT_2pad(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
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

                        buildWT_2pad(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
                    end
                end
            end
            
            if ~isempty(sessions_piezo{mi})
                for j = 1 : length(sessions_piezo{mi})
                    sessionName = sprintf('piezo%d',sessions_piezo{mi}(j));
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
                        if isempty(b_ind)
                            buildWT_2pad(mouseName, sessionName, whisker_d, [], videoFreq, ppm, barRadius)
                        else
                            b_session = b{b_ind};
                            buildWT_2pad(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
                        end
                    end
                end
            end
            
            if ~isempty(sessions_spont{mi})
                for j = 1 : length(sessions_spont{mi})
                    sessionName = sprintf('spont%d',sessions_spont{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        buildWT_2pad(mouseName, sessionName, whisker_d, [], videoFreq, ppm, barRadius)
                    end
                end
            end

        end
    end        
end

%% Check pole up and moving frames, and also pole edges
if testPoleUp
    if all_session == 1
        %% use this code when doing for all of the sessions in the mouse directory 
        for i = 1 : size(mice,2)
            cd(whisker_d)
            sn = dir([mice{i},'S*']);
            sn_pre = dir([mice{i}, 'pre*']);
            if ~isempty(sn)
                for j = 1 : length(sn)
                    if sn(j).isdir
                        [mouseName, sessionName] = strtok(sn(j).name,'S');        
                        if ~isempty(sessionName) 
                            poleUpTest(mouseName,sessionName,videoloc)
                        end
                    end
                end
            end
            if ~isempty(sn_pre)
                for j = 1 : length(sn_pre)
                    if sn_pre(j).isdir
                        [mouseName, sessionName] = strtok(sn_pre(j).name,'pre');
                        poleUpTest(mouseName,sessionName,videoloc)
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
                    poleUpTest(mouseName,sessionName,videoloc)
                end
            end
            if ~isempty(sessions_pre{i})
                for j = 1 : length(sessions_pre{i})
                    mouseName = mice{i};
                    sessionName = sprintf('pre%d',sessions_pre{i}(j));
                    poleUpTest(mouseName,sessionName,videoloc)
                end
            end
        end
    end
end

%% Build WST
if doWST
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
                        elseif strcmp(sessionName, 'S94') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S04'), b));
                        elseif strcmp(sessionName, 'S95') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S05'), b));
                        elseif strcmp(sessionName, 'S99') % JK027
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S17'), b)); 
                        else
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        end
                        b_session = b{b_ind};

                        buildWST_2pad(mouseName, sessionName, whisker_d, b_session, ppm)
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

                    buildWST_2pad(mouseName, sessionName, whisker_d, b_session, ppm)
                end
            end
            
            cd(whisker_d)
            sn_piezo = dir([mice{mi},'piezo*']);
            for si = 1 : length(sn_piezo)
                cd(whisker_d)
                if sn_piezo(si).isdir
                    [mouseName, sessionName] = strtok(sn_piezo(si).name,'piezo');
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
                    if ~isempty(b_ind)
                        b_session = b{b_ind};
                        buildWST_2pad(mouseName, sessionName, whisker_d, b_session, ppm)
                    else
                        buildWST_2pad(mouseName, sessionName, whisker_d, [], ppm)
                    end
                end
            end
            
            cd(whisker_d)
            sn_spont = dir([mice{mi},'spont*']);
            for si = 1 : length(sn_spont)
                cd(whisker_d)
                if sn_spont(si).isdir
                    [mouseName, sessionName] = strtok(sn_spont(si).name,'spont');
                    buildWST_2pad(mouseName, sessionName, whisker_d, [], ppm)
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
                        elseif strcmp(sessionName, 'S94') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S04'), b));
                        elseif strcmp(sessionName, 'S95') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S05'), b));
                        elseif strcmp(sessionName, 'S99') % JK027
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S17'), b)); 
                        else
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        end
                        b_session = b{b_ind};
                        buildWST_2pad(mouseName, sessionName, whisker_d, b_session, ppm)
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

                        buildWST_2pad(mouseName, sessionName, whisker_d, b_session, ppm)
                    end
                end
            end
            
            if ~isempty(sessions_piezo{mi})
                for j = 1 : length(sessions_piezo{mi})
                    sessionName = sprintf('piezo%d',sessions_piezo{mi}(j));
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
                    end
                end
            end
            
            if ~isempty(sessions_spont{mi})
                for j = 1 : length(sessions_spont{mi})
                    sessionName = sprintf('spont%d',sessions_spont{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        buildWST_2pad(mouseName, sessionName, whisker_d, [], ppm)
                    end
                end
            end

        end
    end        
end

%% Perfrom touch_hyperplane
% it includes frame-by-frame estimation of corresponding motor position, based on _WST files
% touch_hyperplane

if makeTouchHyperplane
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
                        elseif strcmp(sessionName, 'S94') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S04'), b));
                        elseif strcmp(sessionName, 'S95') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S05'), b));
                        elseif strcmp(sessionName, 'S99') % JK027
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S17'), b));
                        else
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        end
                        b_session = b{b_ind};
                        run_touch_hyperplane(mouseName, sessionName, b_session, whisker_d, ppm, touchHyperplaneSkip)
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
                    run_touch_hyperplane(mouseName, sessionName, b_session, whisker_d, ppm, touchHyperplaneSkip)
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
                        elseif strcmp(sessionName, 'S94') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S04'), b));
                        elseif strcmp(sessionName, 'S95') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S05'), b));
                        elseif strcmp(sessionName, 'S99') % JK027
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S17'), b));
                        else
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        end
                        b_session = b{b_ind};
                        run_touch_hyperplane(mouseName, sessionName, b_session, whisker_d, ppm, touchHyperplaneSkip)
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
                        run_touch_hyperplane(mouseName, sessionName, b_session, whisker_d, ppm, touchHyperplaneSkip)
                    end
                end
            end
        end
    end    
end
%% Build WL (Finally)
% it includes touch frame calculation

if doWL
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
                        elseif strcmp(sessionName, 'S94') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S04'), b));
                        elseif strcmp(sessionName, 'S95') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S05'), b));
                        elseif strcmp(sessionName, 'S99') % JK027
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S17'), b));
                        else
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        end
                        b_session = b{b_ind};
                        wd = [whisker_d, mouseName, sessionName];                        
                        buildWL_2pad(wd, rInMm, 'b_session', b_session, 'whiskingAmpThreshold', whiskingAmpThreshold, 'stdHistogramThreshold', stdHistogramThreshold, 'distanceHistogramBin', distanceHistogramBin, 'touchBoundaryThickness', touchBoundaryThickness, 'touchBoundaryBuffer', touchBoundaryBuffer, 'maxPointsNearHyperplane', maxPointsNearHyperplane, 'touchKappaSTDthreshold', touchKappaSTDthreshold)
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
                    buildWL_2pad(wd, rInMm, 'b_session', b_session, 'whiskingAmpThreshold', whiskingAmpThreshold, 'stdHistogramThreshold', stdHistogramThreshold, 'distanceHistogramBin', distanceHistogramBin, 'touchBoundaryThickness', touchBoundaryThickness, 'touchBoundaryBuffer', touchBoundaryBuffer, 'maxPointsNearHyperplane', maxPointsNearHyperplane, 'touchKappaSTDthreshold', touchKappaSTDthreshold)
                end
            end
            
%             cd(whisker_d)
%             sn_spont= dir([mice{mi},'spont*']);
%             for si = 1 : length(sn_spont)
%                 cd(whisker_d)
%                 if sn_spont(si).isdir
%                     [mouseName, sessionName] = strtok(sn_spont(si).name,'piezo');
%                     wd = [whisker_d, mouseName, sessionName];
%                     buildWL_2pad(wd, rInMm)
%                 end
%             end
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
                        elseif strcmp(sessionName, 'S94') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S04'), b));
                        elseif strcmp(sessionName, 'S95') % JK052
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S05'), b));
                        elseif strcmp(sessionName, 'S99') % JK027
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,'S17'), b));
                        else
                            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
                        end
                        b_session = b{b_ind};
                        wd = [whisker_d, mouseName, sessionName];
                        buildWL_2pad(wd, rInMm, 'b_session', b_session, 'whiskingAmpThreshold', whiskingAmpThreshold, 'stdHistogramThreshold', stdHistogramThreshold, 'distanceHistogramBin', distanceHistogramBin, 'touchBoundaryThickness', touchBoundaryThickness, 'touchBoundaryBuffer', touchBoundaryBuffer, 'maxPointsNearHyperplane', maxPointsNearHyperplane, 'touchKappaSTDthreshold', touchKappaSTDthreshold)
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
                        buildWL_2pad(wd, rInMm, 'b_session', b_session, 'whiskingAmpThreshold', whiskingAmpThreshold, 'stdHistogramThreshold', stdHistogramThreshold, 'distanceHistogramBin', distanceHistogramBin, 'touchBoundaryThickness', touchBoundaryThickness, 'touchBoundaryBuffer', touchBoundaryBuffer, 'maxPointsNearHyperplane', maxPointsNearHyperplane, 'touchKappaSTDthreshold', touchKappaSTDthreshold)
                    end
                end
            end
            
            if ~isempty(sessions_piezo{mi})
                for j = 1 : length(sessions_piezo{mi})
                    sessionName = sprintf('piezo%d',sessions_piezo{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        wd = [whisker_d, mouseName, sessionName];
                        buildWL_2pad(wd, rInMm)
                    end
                end
            end
            
            if ~isempty(sessions_spont{mi})
                for j = 1 : length(sessions_spont{mi})
                    sessionName = sprintf('spont%d',sessions_spont{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        wd = [whisker_d, mouseName, sessionName];
                        buildWL_2pad(wd, rInMm)
                    end
                end
            end            
        end
    end
end

%% Build 3D reconstruction

errorsession = {};
errorfn = {};

if do3D
    cd(whisker_d)
    if all_session == 1
        for mi = 1 : size(mice,2) % mouse index
            cd(whisker_d)
            sn = dir([whisker_d, mice{mi},'S*']);
            for si = 1 : length(sn)
                if sn(si).isdir
                    [mouseName, sessionName] = strtok(sn(si).name,'S');
                    wd = [whisker_d, mouseName, sessionName];
                    Whisker.makeAllDirectory_Whisker3D_2pad(wd, 'rInMm', rInMm);
                end
            end

            cd(whisker_d)
            sn_pre = dir([mice{mi},'pre*']);
            for si = 1 : length(sn_pre)
                cd(whisker_d)
                if sn_pre(si).isdir
                    [mouseName, sessionName] = strtok(sn_pre(si).name,'pre');
                    wd = [whisker_d, mouseName, sessionName];
                    Whisker.makeAllDirectory_Whisker3D_2pad(wd, 'rInMm', rInMm);
                end
            end

            cd(whisker_d)
            sn_spont= dir([mice{mi},'spont*']);
            for si = 1 : length(sn_spont)
                cd(whisker_d)
                if sn_spont(si).isdir
                    [mouseName, sessionName] = strtok(sn_spont(si).name,'piezo');
                    wd = [whisker_d, mouseName, sessionName];
                    Whisker.makeAllDirectory_Whisker3D_2pad(wd, 'rInMm', rInMm);
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
                        wd = [whisker_d, mouseName, sessionName];
                        errors = Whisker.makeAllDirectory_Whisker3D_2pad(wd, 'rInMm', rInMm);
                        if ~isempty(errors)
                            errorsession{end+1} = [mouseName,sessionName];
                            errorfn{end+1} = errors;
                        end
                    end
                end
            end

            if ~isempty(sessions_pre{mi})
                for j = 1 : length(sessions_pre{mi})
                    sessionName = sprintf('pre%d',sessions_pre{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        wd = [whisker_d, mouseName, sessionName];
                        errors = Whisker.makeAllDirectory_Whisker3D_2pad(wd, 'rInMm', rInMm);
                        if ~isempty(errors)
                            errorsession{end+1} = [mouseName,sessionName];
                            errorfn{end+1} = errors;
                        end
                    end
                end
            end
            
%             if ~isempty(sessions_spont{mi})
%                 for j = 1 : length(sessions_spont{mi})
%                     sessionName = sprintf('spont%d',sessions_spont{mi}(j));
%                     cd(whisker_d)
%                     if exist([mouseName, sessionName],'dir')
%                         wd = [whisker_d, mouseName, sessionName];
%                         errors = Whisker.makeAllDirectory_Whisker3D_2pad(wd, 'rInMm', rInMm);
%                         if ~isempty(errors)
%                             errorsession{end+1} = [mouseName,sessionName];
%                             errorfn{end+1} = errors;
%                         end
%                     end
%                 end
%             end
        end
    end
end
    
%% Build whisker FINAL
if doWF
    cd(whisker_d)
    if all_session == 1
        for mi = 1 : size(mice,2) % mouse index
            cd(whisker_d)
            sn = dir([whisker_d, mice{mi},'S*']);
            for si = 1 : length(sn)                
                if sn(si).isdir && ~contains(sn(si).name, 'spont')
                    [mouseName, sessionName] = strtok(sn(si).name,'S');
                    wd = [whisker_d, mouseName, sessionName];
                    Whisker.makeAllDirectory_WhiskerFinal_2pad(wd);
%                     cd(sn(si).name)                    
%                     wfdir = dir('*_WF_2pad.mat');
%                     wldir = dir('*_WL_2pad.mat');
%                     if length(wldir) > 2 && length(wfdir) < length(wldir)
%                         Whisker.makeAllDirectory_WhiskerFinal_2pad(pwd);
%                     end
                end
%                 cd(whisker_d)
            end

            cd(whisker_d)
            sn_pre = dir([mice{mi},'pre*']);
            for si = 1 : length(sn_pre)
                cd(whisker_d)
                if sn_pre(si).isdir
                    [mouseName, sessionName] = strtok(sn_pre(si).name,'pre');
                    wd = [whisker_d, mouseName, sessionName];
                    Whisker.makeAllDirectory_WhiskerFinal_2pad(wd);
%                     cd(sn(si).name)                    
%                     wfdir = dir('*_WF_2pad.mat');
%                     wldir = dir('*_WL_2pad.mat');
%                     if length(wldir) > 2 && length(wfdir) < length(wldir)
%                         Whisker.makeAllDirectory_WhiskerFinal_2pad(pwd);
%                     end
%                     cd(whisker_d)
                end
            end

%             cd(whisker_d)
%             sn_spont= dir([mice{mi},'spont*']);
%             for si = 1 : length(sn_spont)
%                 cd(whisker_d)
%                 if sn_spont(si).isdir
%                     [mouseName, sessionName] = strtok(sn_spont(si).name,'piezo');
%                     wd = [whisker_d, mouseName, sessionName];
%                     Whisker.makeAllDirectory_Whisker3D_2pad(wd, 'rInMm', rInMm);
%                 end
%             end
        end
    else
        for mi = 1 : size(mice,2) % mouse index
            mouseName = mice{mi};
            if ~isempty(sessions{mi}) 
                for j = 1 : length(sessions{mi})  
                    cd(whisker_d)
                    sessionName = sprintf('S%02d',sessions{mi}(j));
                    if exist([mouseName, sessionName],'dir')
                        wd = [whisker_d, mouseName, sessionName];
                        Whisker.makeAllDirectory_WhiskerFinal_2pad(wd);
                    end
                end
            end

            if ~isempty(sessions_pre{mi})
                for j = 1 : length(sessions_pre{mi})
                    sessionName = sprintf('pre%d',sessions_pre{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        wd = [whisker_d, mouseName, sessionName];
                        Whisker.makeAllDirectory_WhiskerFinal_2pad(wd);                        
                    end
                end
            end
            
%             if ~isempty(sessions_spont{mi})
%                 for j = 1 : length(sessions_spont{mi})
%                     sessionName = sprintf('spont%d',sessions_spont{mi}(j));
%                     cd(whisker_d)
%                     if exist([mouseName, sessionName],'dir')
%                         wd = [whisker_d, mouseName, sessionName];
%                         errors = Whisker.makeAllDirectory_Whisker3D_2pad(wd, 'rInMm', rInMm);
%                         if ~isempty(errors)
%                             errorsession{end+1} = [mouseName,sessionName];
%                             errorfn{end+1} = errors;
%                         end
%                     end
%                 end
%             end
        end
    end
end