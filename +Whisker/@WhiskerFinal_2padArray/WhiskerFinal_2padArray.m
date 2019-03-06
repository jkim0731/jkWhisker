classdef WhiskerFinal_2padArray < handle
    properties
        mouseName = '';
        sessionName = '';
        trials = {};
        trialNums = [];
    end
    
    methods
        function obj = WhiskerFinal_2padArray(d)
            cd(d)
            fnlist = dir('*_WF_2pad.mat');
            
            load(fnlist(end).name)
            obj.mouseName = wf.mouseName;
            if strcmp(wf.sessionName, 'S91')
                wf.sessionName = 'S01';
            end
            obj.sessionName = wf.sessionName;
            
            tns = zeros(length(fnlist),1);
            for i = 1 : length(fnlist)
                tns(i) = str2double(strtok(fnlist(i).name,'_'));
            end
            
            [obj.trialNums, inds] = sort(tns);
            for i = 1 : length(inds)
                currInd = inds(i);
                load(fnlist(currInd).name)
                disp(['Loading ', fnlist(currInd).name, ', ', num2str(i), ' of ', num2str(length(fnlist))])
                obj.trials{i} = wf;
            end
        end
    end
end