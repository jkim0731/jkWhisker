classdef Whisker3D_2padArray < handle
    properties
        mouseName = '';
        sessionName = '';
        trials = {};
    end
    
    methods
        function obj = Whisker3D_2padArray(mouseName, sessionName)
            obj.mouseName = mouseName;
            obj.sessionName = sessionName;
            fnlist = dir('*_W3_2pad.mat');
            for i = 1 : length(fnlist)
                load(fnlist(i).name)
                disp(['Loading ', fnlist(i).name, ', ', num2str(i), ' of ', num2str(length(fnlist))])
                obj.trials{i} = w3;                
            end
        end
    end
end