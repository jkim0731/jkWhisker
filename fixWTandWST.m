%% basic information

%%
%% Error sessions!!
% 025-22
% Need to resume from 027-03. for now, run pre sessions first
% Run JK027S17 WL_2pad again. error at 2nd running. check with wlarray first.
% Run all W3 again. 
%%
%%
%%

% mice = {'JK025','JK027','JK030','JK036','JK037','JK038','JK039','JK041'};

% Error sessions for length(protraction frames) ~= length(prtraction distance) -> Run WL on these first
% {'JK030S10','JK030S13','JK030pre1','JK030pre2','JK036S09','JK036pre1','JK037S03','JK037pre1','JK039S24','JK039S25','JK052S02','JK052S03','JK052S08','JK052S14','JK053S07','JK054S26'}
% and then run 2nd processing of WL. 
% Then check proportion of trials without touch frames again, except for these sessions.



mice = {'JK052', 'JK053','JK054', 'JK056', 'JK070'};




videoloc = 'D:\WhiskerVideo\';
if strcmp(videoloc(end),filesep)
    whisker_d = videoloc;
else
    whisker_d = ([videoloc filesep]);
end
behavior_base_dir = 'D:\SoloData\';

ppm = 17.81;
            % 'pxPerMm': 17.81002608 for telecentric lens
            % /2 for mice <= JK041, because of binning.
% comment out when doing for all of the sessions in the mouse directory
maskmm = 1; % mm from the face to draw the mask
facePosition = 'bottom';
rInMm = 2; % mm from the mask along the whisker to calculate delta kappa
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
%%
%%
%% re-do these
% sessions = {[14,15]};  % for JK039. Did not check JK041 S06~S31 
%%
%%
%%
% sessions = {[4,19,22],[3,16,17],[3,21,22],[1,17,18,91],[7],[2],[1,22:25],[3]};
% sessions = {[3,4,21:23,25,26],[3],[3,4,25,26],[3:5]};

sessions = {[],[],[],[],[],[6:31],[1:25],[1:30]};
sessions_pre = {[],[],[],[],[],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3]};
sessions_spont = {[],[],[],[],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5]};

all_session = 1;
doWT = 1;
doWST = 1;
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

                        buildWT_2pad_error(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
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

                    buildWT_2pad_error(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
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
                        buildWT_2pad_error(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
                    else
                        buildWT_2pad_error(mouseName, sessionName, whisker_d, [], videoFreq, ppm, barRadius)
                    end
                end
            end
            
            cd(whisker_d)
            sn_spont = dir([mice{mi},'spont*']);
            for si = 1 : length(sn_spont)
                cd(whisker_d)
                if sn_spont(si).isdir
                    [mouseName, sessionName] = strtok(sn_spont(si).name,'spont');
                    buildWT_2pad_error(mouseName, sessionName, whisker_d, [], videoFreq, ppm, barRadius)
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
                        buildWT_2pad_error(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
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

                        buildWT_2pad_error(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
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
                            buildWT_2pad_error(mouseName, sessionName, whisker_d, [], videoFreq, ppm, barRadius)
                        else
                            b_session = b{b_ind};
                            buildWT_2pad_error(mouseName, sessionName, whisker_d, b_session, videoFreq, ppm, barRadius)
                        end
                    end
                end
            end
            
            if ~isempty(sessions_spont{mi})
                for j = 1 : length(sessions_spont{mi})
                    sessionName = sprintf('spont%d',sessions_spont{mi}(j));
                    cd(whisker_d)
                    if exist([mouseName, sessionName],'dir')
                        buildWT_2pad_error(mouseName, sessionName, whisker_d, [], videoFreq, ppm, barRadius)
                    end
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

                        buildWST_2pad_error(mouseName, sessionName, whisker_d, b_session, ppm)
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

                    buildWST_2pad_error(mouseName, sessionName, whisker_d, b_session, ppm)
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
                        buildWST_2pad_error(mouseName, sessionName, whisker_d, b_session, ppm)
                    else
                        buildWST_2pad_error(mouseName, sessionName, whisker_d, [], ppm)
                    end
                end
            end
            
            cd(whisker_d)
            sn_spont = dir([mice{mi},'spont*']);
            for si = 1 : length(sn_spont)
                cd(whisker_d)
                if sn_spont(si).isdir
                    [mouseName, sessionName] = strtok(sn_spont(si).name,'spont');
                    buildWST_2pad_error(mouseName, sessionName, whisker_d, [], ppm)
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
                        buildWST_2pad_error(mouseName, sessionName, whisker_d, b_session, ppm)
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

                        buildWST_2pad_error(mouseName, sessionName, whisker_d, b_session, ppm)
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
                        buildWST_2pad_error(mouseName, sessionName, whisker_d, [], ppm)
                    end
                end
            end

        end
    end        
end

