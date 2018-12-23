classdef Whisker3D_2padArray < handle
    properties
        mouseName = '';
        sessionName = '';
        trials = {};
        trialNums = [];
    end
    
    methods
        function obj = Whisker3D_2padArray(d)
            cd(d)
            fnlist = dir('*_W3_2pad.mat');
            
            load(fnlist(1).name)
            obj.mouseName = w3.mouseName;
            obj.sessionName = w3.sessionName;
            
            tns = zeros(length(fnlist),1);
            for i = 1 : length(fnlist)
                tns(i) = str2double(strtok(fnlist(i).name,'_'));
            end
            
            [obj.trialNums, inds] = sort(tns);
            for i = 1 : length(inds)
                currInd = inds(i);
                load(fnlist(currInd).name)
                disp(['Loading ', fnlist(currInd).name, ', ', num2str(i), ' of ', num2str(length(fnlist))])
                obj.trials{i} = w3;
            end
        end
    end
end