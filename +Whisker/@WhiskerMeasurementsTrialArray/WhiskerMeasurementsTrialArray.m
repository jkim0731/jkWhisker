classdef WhiskerMeasurementsTrialArray < handle
%
%
%   WhiskerMeasurementsTrialArray < handle
%
%
% DHO, 11/09.
%
%
    properties
        mouseName = '';
        sessionName = '';
        trials = {};
    end
 
    properties (Dependent = true)
        trialNums
        trialTypes
    end

    methods (Access = public)    
        
        % Make another constructer that takes simply a cell array of strings, each of which gives file name
        % (minus extensions) for each trial of trial array.
        
        function obj = WhiskerMeasurementsTrialArray(tracker_file_name_base, file_nums, behav_trial_nums,varargin)
           % USAGE:
           %    obj = WhiskerMeasurementsTrialArray(cell_array_of_file_name_bases_with_last_element_traj_ids)
           %    obj = WhiskerMeasurementsTrialArray(tracker_file_name_base, file_nums, behav_trial_nums, trajectory_nums, varargin)
           %             
           %
           % trajectory_nums: Can be given in two forms: (a) a vector of 
           %           trajectory numbers (e.g., [0 1 2]); or (b) a cell
           %           array of two elements with a vector of trajectory numbers
           %           in the first element, and a cell array of corresponding
           %           whisker names in the second element (e.g., {[0 1 2],{'D4,'D3','D2'}}).
           %           In the second case the number of trajectory IDs and whisker
           %           whisker names must match.
           % 
           % varargin{1}: optional mouse_name string. 
           % varargin{2}: optional session_name string. 
           %
           % varargin{3}: optional vector of integers specifying trial type (e.g., 0 for NoGo,
           %               1 for Go). Must be of same length as behav_trial_nums and be matched
           %               element for element.
           %
           % Still need to add some error checking on inputs.
            
            if nargin==0
                return
            end
            
            %------------------------------------------------
            % Single argument: must be cell array of strings giving file names bases (file name minus extension), 
            % except for last element of cell array which must be a vector of trajectory IDs.
            %
            if nargin==1 
                if ~iscell(tracker_file_name_base)
                    error('Argument is not a cell array of strings giving file name bases.')
                elseif isempty(tracker_file_name_base)
                    error('Cell array argument cannot be empty.')
                end
                traj_nums = tracker_file_name_base{end};
                if ~isnumeric(traj_nums)
                    error('Last element of (cell array) argument must be a vector of trajectory IDs.')
                end
                
                ntrials = numel(tracker_file_name_base)-1;
                obj.trials = cell(1, ntrials);
                for k=1:ntrials
                    filename = [tracker_file_name_base{k} '.measurements'];
                    obj.trials{k} = Whisker.WhiskerMeasurementsTrial(filename);
                end
                return
            end
            %------------------------------------------------
            
            
            if nargin > 3
                obj.mouseName = varargin{1};
            end
            if nargin > 4
                obj.sessionName = varargin{2}; 
            end
           
            if nargin > 5
                trial_types = double(varargin{3});
                if ~isnumeric(trial_types)
                    error('Trial type argument (varargin{3}) must be an integer.')
                end
                if numel(trial_types) ~= numel(behav_trial_nums)
                    error('Inputs varargin{3} (trial_types) and behav_trial_nums must be the same length.')
                end
            else
                trial_types = [];
            end
            
            if numel(file_nums) ~= numel(behav_trial_nums)
                error('Inputs file_nums and behav_trial_nums must be the same length.')
            end
            
            obj.trials = cell(1, length(file_nums));
            
            for k=1:length(file_nums)
                fn = file_nums(k);
                if fn < 10
                    filler = '00';
                elseif fn < 100
                    filler = '0';
                elseif fn < 1000
                    filler = '';
                end

                filename = [tracker_file_name_base filler int2str(fn) '.measurements'];
                
                obj.trials{k} = Whisker.WhiskerMeasurementsTrial(filename);
                obj.trials{k}.trialNum = behav_trial_nums(k);
                obj.trials{k}.mouseName = obj.mouseName;
                obj.trials{k}.sessionName = obj.sessionName;
                if ~isempty(trial_types)
                    obj.trials{k}.trialType = trial_types(k);  
                end
            end
        end
                              
        function r = length(obj)
            r = length(obj.trials);
        end
 
    end
    
   
    methods % Dependent property methods; cannot have attributes.
        function value = get.trialNums(obj)
            if ~isempty(obj.trials)
                value = cellfun(@(x) x.trialNum, obj.trials);
            else
                value = [];
            end
        end
        
        function value = get.trialTypes(obj)
            if ~isempty(obj.trials)
                value = cellfun(@(x) x.trialType, obj.trials);
            else
                value = [];
            end
        end
        
    end
    

end

