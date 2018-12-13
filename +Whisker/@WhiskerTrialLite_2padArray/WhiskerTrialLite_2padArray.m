classdef WhiskerTrialLite_2padArray < handle
    properties
        mouseName = '';
        sessionName = '';
        trials = {};
        trialNums = [];
    end
    
    methods
        function obj = WhiskerTrialLite_2padArray(d, mouseName, sessionName)
            cd(d)
            obj.mouseName = mouseName;
            obj.sessionName = sessionName;
            fnlist = dir('*_WL_2pad.mat');
            tns = zeros(length(fnlist),1);
            for i = 1 : length(fnlist)
                tns(i) = str2double(strtok(fnlist(i).name,'_'));
            end
            
            [obj.trialNums, inds] = sort(tns);
            for i = 1 : length(inds)
                currInd = inds(i);
                load(fnlist(currInd).name)
                disp(['Loading ', fnlist(currInd).name, ', ', num2str(i), ' of ', num2str(length(fnlist)), ' from ', mouseName, ' ', sessionName])
                obj.trials{i} = wl;                
            end
        end
    end
end
