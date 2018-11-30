mice = {'JK025','JK027','JK030','JK036','JK037','JK038','JK039','JK041','JK052','JK053','JK054','JK056'};
videoloc = 'D:\WhiskerVideo\';
if strcmp(videoloc(end),filesep)
    whisker_d = videoloc;
else
    whisker_d = ([videoloc filesep]);
end
sessions_pre = 1:3;

sessions = 1:31;

trackingErrors = cell(length(mice),34);
trackingErrorsProps = cell(length(mice),34);

for mi = 1 : length(mice)
    for si = sessions
        sessionName = sprintf('S%02d',sessions(si));
        trackingErrors{mi,si} = [];
        trackingErrorsProps{mi,si} = [];
        try
            cd([videoloc, mice{mi}, sessionName])
            wsa = Whisker.WhiskerSignalTrialArray_2pad(pwd);            
            for wi = 1 : length(wsa.trials)
                ws = wsa.trials{wi};
%                 noNaNInd = intersect(find(~isnan(sum(ws.whiskerEdgeCoord,2))), ws.poleUpFrames);
                if length(intersect(ws.time{1}, ws.time{2})) < ws.nof * 0.9
                    trackingErrors{mi,si} = [trackingErrors{mi,si}, str2double(ws.trackerFileName)];
                    trackingErrorsProps{mi,si} = [trackingErrorsProps{mi,si}, length(intersect(ws.time{1}, ws.time{2})) / ws.nof];
                end
            end
        catch
        end
    end
    for si = sessions_pre
        sessionName = sprintf('pre%d',sessions_pre(si));
        trackingErrors{mi,si+31} = [];
        trackingErrorsProps{mi,si} = [];
        try
            cd([videoloc, mice{mi}, sessionName])
            wsa = Whisker.WhiskerSignalTrialArray_2pad(pwd);
            for wi = 1 : length(wsa.trials)
                ws = wsa.trials{wi};
%                 noNaNInd = intersect(find(~isnan(sum(ws.whiskerEdgeCoord,2))), ws.poleUpFrames);
                if length(intersect(ws.time{1}, ws.time{2})) < ws.nof * 0.9
                    trackingErrors{mi,si+31} = [trackingErrors{mi,si+31}, str2double(ws.trackerFileName)];
                    trackingErrorsProps{mi,si} = [trackingErrorsProps{mi,si}, length(intersect(ws.time{1}, ws.time{2})) / ws.nof];
                end
            end
        catch
        end
    end
end
                

