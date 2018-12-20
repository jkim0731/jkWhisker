classdef WhiskerTrialLite_2padArray < handle
    properties
        mouseName = '';
        sessionName = '';
        trials = {};
        trialNums = [];
    end
    
    methods
        function obj = WhiskerTrialLite_2padArray(d)
            cd(d)
            fnlist = dir('*_WL_2pad.mat');
            tns = zeros(length(fnlist),1);
            for i = 1 : length(fnlist)
                tns(i) = str2double(strtok(fnlist(i).name,'_'));
            end
            load(fnlist(1).name)
            obj.mouseName = wl.mouseName;
            obj.sessionName = wl.sessionName;
            
            [obj.trialNums, inds] = sort(tns);
            for i = 1 : length(inds)
                currInd = inds(i);
                load(fnlist(currInd).name)
                disp(['Loading ', fnlist(currInd).name, ', ', num2str(i), ' of ', num2str(length(fnlist)), ' from ', obj.mouseName, ' ', obj.sessionName])
                obj.trials{i} = wl;                
            end
        end
    end
end
