classdef WhiskerTrialLite_2padArray < handle
    properties
        mouseName = '';
        sessionName = '';
        trials = {};
    end
    
    methods
        function obj = WhiskerTrialLite_2padArray(mouseName, sessionName)
            obj.mouseName = mouseName;
            obj.sessionName = sessionName;
            fnlist = dir('*_WL_2pad.mat');
            for i = 1 : length(fnlist)
                load(fnlist(i).name)
                disp(['Loading ', fnlist(i).name, ', ', num2str(i), ' of ', num2str(length(fnlist))])
                obj.trials{i} = wl;                
            end
        end
    end
end
