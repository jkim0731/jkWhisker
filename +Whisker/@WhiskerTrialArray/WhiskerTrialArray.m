classdef WhiskerTrialArray < handle
%
%
%   WhiskerTrialArray < handle
%
%
% DHO, 5/08.
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
        
        function obj = WhiskerTrialArray(tracker_file_name_base, file_nums, behav_trial_nums, trajectory_nums, varargin)
           % USAGE:
           %    obj = WhiskerTrialArray(cell_array_of_file_name_bases_with_last_element_traj_ids)
           %    obj = WhiskerTrialArray(tracker_file_name_base, file_nums, behav_trial_nums, trajectory_nums, varargin)
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
                    filename = tracker_file_name_base{k};
                    obj.trials{k} = Whisker.WhiskerTrial(filename, NaN, traj_nums);
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
           
            if nargin > 6
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

                filename = [tracker_file_name_base filler int2str(fn)];

                trial_num = behav_trial_nums(k);
                
                if isempty(trial_types)
                    obj.trials{k} = Whisker.WhiskerTrial(filename, trial_num, trajectory_nums, obj.mouseName, obj.sessionName);
                else
                    obj.trials{k} = Whisker.WhiskerTrial(filename, trial_num, trajectory_nums, obj.mouseName, obj.sessionName, trial_types(k));
                end
            end
        end
        
        function r = saveobj(obj)
            % 
            %   r = saveobj(obj)
            %
            % Workaround for some documented MATLAB bugs with hierarchical file saving
            %   and for excessively sized files in -V7.3 saving.
            %
            r.mouseName = obj.mouseName;
            r.sessionName = obj.sessionName;
            r.trials = cellfun(@(x) x.saveobj, obj.trials,'UniformOutput',false);
        end
                       
        function r = length(obj)
            r = length(obj.trials);
        end
                
        function r = get_whisker_y_range(obj,tid,varargin)
            %
            %   r = get_whisker_y_range(obj,tid,varargin)
            %
            % Returns a cell array giving in each element, for each trial in this
            % WhiskerTrialArray, the minimum and maximum y-values
            % for whisker tid.
            %
            % r: [ymin ymax], where ymax is toward bottom of
            %      original image, ymin toward top.
            %
            % tid: Either trajectory ID (as an integer) or whisker name 
            %      (as a string).
            %
            % varargin{1} - x-limits in form [xmin xmax].
            %               Can be empty ([]) to allow access to varargin{2}. 
            %
            % varargin{2}: Optional 2 x 1 vector giving starting and ending times (in seconds) to include in
            %               range, inclusive, starting with 0.
            %               Of format: [startTimeInSec endTimeInSec].
            %
            %
            %
            if ~isempty(obj.trials)
                if nargin > 2
                    XValLims = varargin{1};
                else
                    XValLims = [];
                end
                if nargin > 3
                    restrictTime = varargin{2};
                else
                    restrictTime = [];
                end
                r = cellfun(@(x) x.get_whisker_y_range(tid,XValLims,restrictTime), obj.trials,'UniformOutput',false);
            else
                r = [];
            end
        end
        
        function obj = fit_polys(obj)
            % See 'help Whisker.WhiskerTrial.fit_polys'
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                for k=1:ntrials
                    disp(['Fitting for trial number ' int2str(k) '/' int2str(ntrials)])
                    obj.trials{k}.fit_polys;
                end
            end
        end
        
        function obj = fit_polys_roi(obj,rad_lims)
            % See 'help Whisker.WhiskerTrial.fit_polys_roi'
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                for k=1:ntrials
                    disp(['Fitting for trial number ' int2str(k) '/' int2str(ntrials)])
                    obj.trials{k}.fit_polys_roi(rad_lims);
                end
            end
        end
        
        function obj = set_whisker_names_all(obj, names)
            if ~isempty(obj.trials)
                for k=1:length(obj.trials)
                    if numel(names) ~= numel(obj.trials{k}.trajectoryIDs)
                        error(['Number of elements in argument ''names'' and number of trajectoryIDs ' ...
                            'are not equal for trialNum ' int2str(obj.trials{k}.trialNum)])
                    else
                        obj.trials{k}.whiskerNames = names;
                    end
                end
            end
        end
        
        function obj = set_faceSideInImage(obj, face_side)
            %
            %   obj = set_faceSideInImage(obj, face_side)
            %
            % INPUTS:
            %   face_side:  A string; one of: 'top','bottom','left','right'. 
            %               This must be set properly for the follicle-side of the
            %               traced whisker to be determined, and to compute theta properly.
            %
            if isempty(strmatch(face_side,{'top','bottom','left','right'}))
                error('Invalid value of argument ''face_side''.')
            end
            if ~isempty(obj.trials)
                for k=1:length(obj.trials)
                    obj.trials{k}.faceSideInImage = face_side;
                end
            end
        end

        function obj = set_pxPerMm(obj, px_per_mm)
            %
            %   obj = set_pxPerMm(obj, px_per_mm)
            %
            % INPUTS:
            %    A scaler giving the number of pixels per mm for 
            %    the videographic conditions of this trial.
            %
            if ~isempty(obj.trials)
                for k=1:length(obj.trials)
                    obj.trials{k}.pxPerMm = px_per_mm;
                end
            end
        end
        
        function obj = set_framePeriodInSec(obj, frame_period_in_sec)
            %
            %   obj = set_framePeriodInSec(obj, frame_period_in_sec)
            %
            % INPUTS:
            %    A scaler giving the frame period in sec for
            %    the videographic conditions of these trials.
            %
            if ~isempty(obj.trials)
                for k=1:length(obj.trials)
                    obj.trials{k}.framePeriodInSec = frame_period_in_sec;
                end
            end
        end
        
        function obj = set_imagePixelDimsXY(obj, image_dims)
            %
            %   obj = set_imagePixelDimsXY(obj, image_dims)
            %
            % INPUTS:
            %   image_dims:  A 1x2 vector giving [x,y] dimensions of videographic
            %                image used in the whisker tracker.
            %               This must be set properly for the follicle-side of the
            %               traced whisker to be determined.
            %
            if ~isnumeric(image_dims)
                error('Argument image_dims must be a 1x2 vector.')
            end
            if size(image_dims,1) ~= 1 || size(image_dims,2) ~= 2
                error('Argument image_dims must be a 1x2 vector.')
            end
            if ~isempty(obj.trials)
                for k=1:length(obj.trials)
                    obj.trials{k}.imagePixelDimsXY = image_dims;
                end
            end
        end
        
        function obj = set_barRadius(obj, bar_radius_in_pix)
            %
            % obj = set_barRadius(obj, bar_radius_in_pix)
            %
            % Sets barRadius property of each WhiskerTrial in this WhiskerTrialArray.
            %
            if ~isempty(obj.trials)
                for k=1:length(obj.trials)
                    obj.trials{k}.barRadius = bar_radius_in_pix;
                end
            end
        end
        
        function obj = set_bar_offset(obj,dx,dy)
            %
            % obj = set_bar_offset(obj,dx,dy)
            %
            % dx: Number of pixels to offset bar center in x.
            % dy: Number of pixels to offset bar center in y.
            %
            % dx,dy must be the same size.
            %
            % If dx,dy are scalers, the bar offset for each trial will
            % be set to dx,dy.
            % 
            % If dx,dy each has as many elements as this WhiskerTrialArray,
            % each trial will be giving the corresponding dx,dy bar offset.
            %
            %
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                if length(dx)==1 && length(dy)==1
                    dx = repmat(dx,ntrials,1);
                    dy = repmat(dy,ntrials,1);
                elseif length(dx)~=ntrials || length(dy)~=ntrials
                    error('Arguments dx and/or dy are not the correct length.')
                end
                
                for k=1:ntrials
                    obj.trials{k}.set_bar_offset(dx(k),dy(k));
                end
            end
        end
        
        function obj = set_mask_from_points(obj,tid,x,y)
            %
            % For each WhiskerTrial in this WhiskerTrialArray,
            % sets obj.polyFitsMask in order
            % to create a mask defined by the points in x and y.
            %
            % Allows only setting identical masks for each trial for a
            % given trajectory ID.
            %
            % tid: Trajectory ID. Can be a vector with multiple trajectory
            %       IDs. In this case all will be set to have same mask.
            %
            % x: Row vector of x coordinates to define mask.
            % y: Row vector of y coordinates to define mask.
            %
            % If N points are selected, mask will be the (N-1)-th
            % degree polynomial fit to the points for N < 6. For N >= 6
            % the polynomial will be 5-th degree.
            %
            qnum = length(x);
            if length(x) ~= length(y)
                error('Inputs x and y must be of equal length.')
            end
            
            % Make x, y row vectors:
            if size(x,1) > size(x,2)
                x = x';
            end
            if size(y,1) > size(y,2)
                y = y';
            end
            
            if qnum < 2
                error('Must define at least 2 points.')
            elseif qnum < 6
                polyDegree = qnum-1;
            else
                polyDegree = 5;
            end
            
            q = (0:(qnum-1))./(qnum-1); % [0,1]
            
            % polyfit() gives warnings that indicate that we don't need such a high degree
            % polynomials. Turn off.
            warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
            px = polyfit(q,x,polyDegree);
            py = polyfit(q,y,polyDegree);
            warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');
            
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                for h=1:ntrials
                    if isempty(obj.trials{h}.polyFitsMask)
                        obj.trials{h}.polyFitsMask = cell(1,length(obj.trials{h}.trajectoryIDs));
                    end
                    
                    for k=1:length(tid)
                        ind = obj.trials{h}.trajectoryIDs==tid(k);
                        if max(ind) < 1
                            error('Trajectory ID was not found.')
                        end
                        obj.trials{h}.polyFitsMask{ind} = {px,py};
                    end
                end
            end
        end
        
    end
    
    methods (Static = true)
        function obj = loadobj(r)
            %
            %   obj = loadobj(r)
            %
            obj = Whisker.WhiskerTrialArray;
            obj.mouseName = r.mouseName;
            obj.sessionName = r.sessionName;
            obj.trials = r.trials;
            ntrials = length(r.trials);
            for k=1:ntrials
                t = Whisker.WhiskerTrial;
                obj.trials{k} = t.loadobj(r.trials{k});
            end
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

