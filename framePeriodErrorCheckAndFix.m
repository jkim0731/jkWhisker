% mice = {'JK025','JK027','JK030','JK036','JK037','JK038','JK039','JK041','JK052','JK053','JK054','JK056'};
videoloc = 'D:\WhiskerVideo\';
% if strcmp(videoloc(end),filesep)
%     whisker_d = videoloc;
% else
%     whisker_d = ([videoloc filesep]);
% end
% 
% sessions = {[1:100],[1:100],[1:100],[1:100],[1:100],[1:100],[1:100],[1:100],[1:100],[1:100],[1:100],[1:100]};
% sessions_pre = {[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3],[1:3]};
% sessions_piezo = {[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5]};
% sessions_spont = {[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5],[1:5]};
% 
framePeriodErrorSessionsWT = {};
framePeriodErrorSessionsWSTtype1 = {};
framePeriodErrorSessionsWSTtype2 = {};
framePeriodErrorSessionsWSTtype3 = {};
framePeriodErrorSessionsWLtype1 = {};
framePeriodErrorSessionsWLtype2 = {};
framePeriodErrorSessionsWLtype3 = {};
framePeriodErrorSessionsW3type1 = {};
framePeriodErrorSessionsW3type2 = {};
framePeriodErrorSessionsW3type3 = {};
cd(videoloc)
dirlist = dir('JK*');
for i = 1 : length(dirlist)
    if dirlist(i).isdir
        cd([videoloc, dirlist(i).name])
        wtlist = dir('*_WT.mat');
        if ~isempty(wtlist)
            load(wtlist(1).name)
            if w.framePeriodInSec > 1
                framePeriodErrorSessionsWT{end+1} = dirlist(i).name;
                for j = 1 : length(wtlist)
                    load(wtlist(j).name)
                    w.framePeriodInSec = 1 / w.framePeriodInSec;
                    save(wtlist(j).name, 'w')
                end
            end
        end
        
        wstlist = dir('*_WST.mat');
        if ~isempty(wstlist)
            load(wstlist(1).name)
            if ws.framePeriodInSec > 1
                if mean(diff(ws.time{1})) > 10
                    framePeriodErrorSessionsWSTtype1{end+1} = dirlist(i).name;
                    for j = 1 : length(wstlist)
                        load(wstlist(j).name)
                        ws.time{1} = ws.time{1}/(ws.framePeriodInSec^2);
                        ws.time{2} = ws.time{2}/(ws.framePeriodInSec^2);
                        ws.framePeriodInSec = 1 / ws.framePeriodInSec;
                        save(wstlist(j).name, 'ws');
                    end
                else
                    framePeriodErrorSessionsWSTtype2{end+1} = dirlist(i).name;
                    for j = 1 : length(wstlist)
                        load(wstlist(j).name)
                        ws.framePeriodInSec = 1 / ws.framePeriodInSec;
                        save(wstlist(j).name, 'ws');
                    end
                end
            else
                if mean(diff(ws.time{1})) > 10
                    framePeriodErrorSessionsWSTtype3{end+1} = dirlist(i).name;
                    for j = 1 : length(wstlist)
                        load(wstlist(j).name)
                        ws.time{1} = ws.time{1}*(ws.framePeriodInSec^2);
                        ws.time{2} = ws.time{2}*(ws.framePeriodInSec^2);
                        save(wstlist(j).name, 'ws');
                    end
                end
            end
        end
        
        wllist = dir('*_WL_2pad.mat');
        if ~isempty(wllist)
            load(wllist(1).name)
            if wl.framePeriodInSec > 1
                if mean(diff(wl.time{1})) > 10
                    framePeriodErrorSessionsWLtype1{end+1} = dirlist(i).name;
                    for j = 1 : length(wllist)
                        load(wllist(j).name)
                        wl.time{1} = wl.time{1}/(wl.framePeriodInSec^2);
                        wl.time{2} = wl.time{2}/(wl.framePeriodInSec^2);
                        wl.framePeriodInSec = 1 / wl.framePeriodInSec;
                        save(wllist(j).name, 'wl')
                    end
                else
                    framePeriodErrorSessionsWLtype2{end+1} = dirlist(i).name;
                    for j = 1 : length(wllist)
                        load(wllist(j).name)
                        wl.framePeriodInSec = 1 / wl.framePeriodInSec;
                        save(wllist(j).name, 'wl')
                    end
                end
            else
                if mean(diff(wl.time{1})) > 10
                    framePeriodErrorSessionsWLtype3{end+1} = dirlist(i).name;
                    for j = 1 : length(wllist)
                        load(wllist(j).name)
                        wl.time{1} = wl.time{1}*(wl.framePeriodInSec^2);
                        wl.time{2} = wl.time{2}*(wl.framePeriodInSec^2);
                        save(wllist(j).name, 'wl')
                    end
                end
            end
        end
            
        w3list = dir('*_W3_2pad.mat');
        if ~isempty(w3list)
            load(w3list(1).name)
            if w3.framePeriodInSec > 1
                if mean(diff(w3.time)) > 10
                    framePeriodErrorSessionsW3type1{end+1} = dirlist(i).name;
                    for j = 1 : length(w3list)
                        load(w3list(j).name)
                        w3.time = w3.time/(w3.framePeriodInSec^2);                        
                        w3.framePeriodInSec = 1 / w3.framePeriodInSec;
                        save(w3list(j).name, 'w3')
                    end
                else
                    framePeriodErrorSessionsW3type2{end+1} = dirlist(i).name;
                    for j = 1 : length(w3list)
                        load(w3list(j).name)
                        w3.framePeriodInSec = 1 / w3.framePeriodInSec;
                        save(w3list(j).name, 'w3')
                    end
                end
            else
                if mean(diff(w3.time)) > 10
                    framePeriodErrorSessionsW3type3{end+1} = dirlist(i).name;
                    for j = 1 : length(w3list)
                        load(w3list(j).name)
                        w3.time = w3.time*(w3.framePeriodInSec^2);
                        save(w3list(j).name, 'w3')
                    end
                end
            end
        end
    end
end
% for mi = 1 : length(mice)
%     for si = sessions{mi}
%         sessionName = sprintf('S%02d',sessions{mi}(si));
%         try
%             cd([videoloc, mice{mi}, sessionName])
%             wtflist = dir('*_
%             wsa = Whisker.WhiskerSignalTrialArray_2pad(pwd);            
%             for wi = 1 : length(wsa.trials)
%                 ws = wsa.trials{wi};
% %                 noNaNInd = intersect(find(~isnan(sum(ws.whiskerEdgeCoord,2))), ws.poleUpFrames);
%                 if length(intersect(ws.time{1}, ws.time{2})) < ws.nof * 0.9
%                     trackingErrors{mi,si} = [trackingErrors{mi,si}, str2double(ws.trackerFileName)];
%                     trackingErrorsProps{mi,si} = [trackingErrorsProps{mi,si}, length(intersect(ws.time{1}, ws.time{2})) / ws.nof];
%                 end
%             end
%         catch
%         end
%     end
%     for si = sessions_pre{mi}
%         sessionName = sprintf('pre%d',sessions_pre{mi}(si));
%         trackingErrors{mi,si+31} = [];
%         trackingErrorsProps{mi,si} = [];
%         try
%             cd([videoloc, mice{mi}, sessionName])
%             wsa = Whisker.WhiskerSignalTrialArray_2pad(pwd);
%             for wi = 1 : length(wsa.trials)
%                 ws = wsa.trials{wi};
% %                 noNaNInd = intersect(find(~isnan(sum(ws.whiskerEdgeCoord,2))), ws.poleUpFrames);
%                 if length(intersect(ws.time{1}, ws.time{2})) < ws.nof * 0.9
%                     trackingErrors{mi,si+31} = [trackingErrors{mi,si+31}, str2double(ws.trackerFileName)];
%                     trackingErrorsProps{mi,si} = [trackingErrorsProps{mi,si}, length(intersect(ws.time{1}, ws.time{2})) / ws.nof];
%                 end
%             end
%         catch
%         end
%     end
%     for si = sessions_piezo{mi}
%         sessionName = sprintf('piezo%d',sessions_piezo{mi}(si));
%         trackingErrors{mi,si+31} = [];
%         trackingErrorsProps{mi,si} = [];
%         try
%             cd([videoloc, mice{mi}, sessionName])
%             wsa = Whisker.WhiskerSignalTrialArray_2pad(pwd);
%             for wi = 1 : length(wsa.trials)
%                 ws = wsa.trials{wi};
% %                 noNaNInd = intersect(find(~isnan(sum(ws.whiskerEdgeCoord,2))), ws.poleUpFrames);
%                 if length(intersect(ws.time{1}, ws.time{2})) < ws.nof * 0.9
%                     trackingErrors{mi,si+31} = [trackingErrors{mi,si+31}, str2double(ws.trackerFileName)];
%                     trackingErrorsProps{mi,si} = [trackingErrorsProps{mi,si}, length(intersect(ws.time{1}, ws.time{2})) / ws.nof];
%                 end
%             end
%         catch
%         end
%     end
%     for si = sessions_spont{mi}
%         sessionName = sprintf('spont%d',sessions_spont{mi}(si));
%         trackingErrors{mi,si+31} = [];
%         trackingErrorsProps{mi,si} = [];
%         try
%             cd([videoloc, mice{mi}, sessionName])
%             wsa = Whisker.WhiskerSignalTrialArray_2pad(pwd);
%             for wi = 1 : length(wsa.trials)
%                 ws = wsa.trials{wi};
% %                 noNaNInd = intersect(find(~isnan(sum(ws.whiskerEdgeCoord,2))), ws.poleUpFrames);
%                 if length(intersect(ws.time{1}, ws.time{2})) < ws.nof * 0.9
%                     trackingErrors{mi,si+31} = [trackingErrors{mi,si+31}, str2double(ws.trackerFileName)];
%                     trackingErrorsProps{mi,si} = [trackingErrorsProps{mi,si}, length(intersect(ws.time{1}, ws.time{2})) / ws.nof];
%                 end
%             end
%         catch
%         end
%     end    
% end