classdef WhiskerTrial_nx < Whisker.WhiskerTrial
    %
    %   subclass of Whisker.WhiskerTrial
    %
    %   - NX 2/2010
    % 
    %
    %
    %
    properties
        touch_on = {};
        touch_off = {};
        n_touch = [];
        whiskerPadOrigin_nx = [145 367];
        bar_pos_trial = []; % Single coordinate of bar center for the current trial
                            %     In contrast to bar_pos, containing bar
                            %     coordinates for all frames.
    end
    
    properties (Dependent = true)
%         numFramesEachTrajectory
%         allFrameNums
%         numUniqueFrames
%         whiskerPadOrigin
    end
    
    methods (Access = public)
        function obj = WhiskerTrial_nx(tracker_file_name, trial_num, trajectory_nums, varargin)
            %
            %   obj = WhiskerTrial(tracker_file, trial_num, trajectory_nums)
            %   obj = WhiskerTrial(tracker_file, trial_num, trajectory_nums, mouse_name, session_name)
            %
            % trajectory_nums: Can be given in two forms: (a) a vector of
            %           trajectory numbers (e.g., [0 1 2]); or (b) a cell
            %           array of two elements with a vector of trajectory numbers
            %           in the first element, and a cell array of corresponding
            %           whisker names in the second element (e.g., {[0 1 2],{'D4,'D3','D2'}}).
            %           In the second case the number of trajectory IDs and whisker
            %           whisker names must match.
            %
            % varargin{1}: bar coordinate for all frames of the current
            %           trial. 
            %           
            %
            % varargin{2}: Optional mouse_name and session_name strings. If
            %           one is present both must be.
            %
            % varargin{3}: Optional integer-specified trial type (e.g., 0 for NoGo,
            %               1 for Go).
            %
            %
            %
            if nargin==0
                return
            end
            if ~ischar(tracker_file_name)
                error('First argument must be valid file name given as a string.')
            end
            if ~isnumeric(trial_num)
                error('Second argument must be an integer trial number or NaN placeholder.')
            elseif isnan(trial_num)
                % do nothing
            elseif mod(trial_num,1) ~= 0
                error('Second argument must be an integer trial number or NaN placeholder.')
            end
            if nargin > 3
                if ~isnumeric(varargin{1})
                    error('the 4th argument must be numemric.')
                end
                
                obj.bar_pos_trial = varargin{1};
            end
            if nargin > 4
                obj.mouseName = varargin{2}{1};
                obj.sessionName = varargin{2}{2};
            end
            
            if nargin > 5
                obj.trialType = varargin{3};
                if ~isnumeric(obj.trialType)
                    error('Trial type argument (varargin{3}) must be an integer.')
                end
            end
            
            if iscell(trajectory_nums)
                obj.trajectoryIDs = trajectory_nums{1};
                obj.whiskerNames = trajectory_nums{2};
                if numel(obj.trajectoryIDs) ~= numel(obj.whiskerNames)
                    error('Unequal number of whisker names and trajectory IDs in argument trajectory_nums')
                end
            elseif isnumeric(trajectory_nums)
                obj.trajectoryIDs = trajectory_nums;
            else
                error('Argument trajectory_nums is an invalid type.')
            end
            
            obj.trackerFileName = tracker_file_name;
            obj.trialNum = trial_num;
%           %%%% Redefine the properties from the superclass
            obj.faceSideInImage = 'bottom'; % 'top', 'bottom', 'left','right'
            obj.imagePixelDimsXY = [560 360]; % [NumberOfXPixels NumberOfYPixels]
        
            obj.pxPerMm = 25.7; %22.68;
%             obj.whiskerPadOrigin_nx = [193 318];
            try
%                 [r, obj.trackerFileFormat] = Whisker.load_segments([tracker_file_name '.whiskers']);
                [r, obj.trackerFileFormat] = Whisker.load_whiskers_file([tracker_file_name '.whiskers']);
               
                % .measurements file is newer replacement for .trajectories.  If there's a .measurements
                % file, choose that.  If not, choose the .trajectories.  Give a message to alert user
                % if both are found.
                if exist([tracker_file_name '.measurements'],'file')
                    M = Whisker.read_whisker_measurements([tracker_file_name '.measurements']);
                    trajectory_ids = M(:,1); 
                    frame_nums = M(:,2);
                    segment_nums = M(:,3);
                    total_whisker_num = length(obj.trajectoryIDs);
                    for k = 1:total_whisker_num
                        % find the annotated frame for touch onset
                        obj.touch_on{k} = sort(M(M(:,1)== obj.trajectoryIDs(k) + 5, 2));
                        % find the annotated frame for touch offset
                        obj.touch_off{k} = sort(M(M(:,1)== obj.trajectoryIDs(k) - 5, 2));
                    end
                    
                    if exist([tracker_file_name '.trajectories'],'file')
                       disp(['For ' tracker_file_name 'found both .measurements and .trajectories files---using .measurements.'])
                    end
                else 
                    % .measurements file not found; choose .trajectories file.
                    [trajectory_ids, frame_nums, segment_nums] = Whisker.load_trajectories([tracker_file_name '.trajectories']);
                end
                
                if exist([tracker_file_name '.bar'],'file')
                    [bar_f, bar_x, bar_y] = Whisker.load_bar([tracker_file_name '.bar']);
                    obj.barPos = [bar_f, bar_x, bar_y];
                end
            catch ME
                disp(ME)
                error(['Cannot load tracker files for: ' tracker_file_name])
            end
            
            sFrameNums = cellfun(@(x) x{1},r);
            
            D = cell(1,length(trajectory_ids));
            
            for k=1:length(trajectory_ids)
                frame = frame_nums(k);
                trajectory = trajectory_ids(k);
                segment = segment_nums(k);
                
                indFrame = find(sFrameNums==frame);
                
                trialSegs = r{indFrame}{2};
                indSeg = find(trialSegs==segment);
               
                if isempty(indSeg)
                    xdat = single([]); ydat = single([]);
                else
                    if strcmp(obj.trackerFileFormat,'whisker0')
                        xdat = single(r{indFrame}{3}(indSeg,:));
                    else
                        xdat = single(r{indFrame}{3}{indSeg});
                    end
                    ydat = single(r{indFrame}{4}{indSeg});
                end
                
                D{k} = {trajectory, frame, segment, xdat, ydat};
            end
            
            %--- .trajectories files don't necessarily have frameNums in order,
            % so sort them here:
            f =  cellfun(@(x) x{2}, D);
            [tmp,ind] = sort(f,2,'ascend');
            D = D(ind);
            %-----------------------------------------------------------
            
            tr = cellfun(@(x) x{1}, D);
            ind = ismember(tr, obj.trajectoryIDs);
            D = D(ind); % restrict to trajectories specified in trajectory_nums input parameter.
            tr = tr(ind);
            
            % Package by trajectory:
            obj.trackerData = cell(1,length(obj.trajectoryIDs));
            for k=1:length(obj.trajectoryIDs)
                obj.trackerData{k} = D(tr==obj.trajectoryIDs(k));
            end
        end
        
%          function r = saveobj(obj)
%          % Call superclass saveobj method
%             r = saveobj@super(obj); 
%          % Perform subclass save operations
%          end
        
    end
    
    methods (Static = true)
%         function obj = loadobj(r)
%             %
%             %   obj = loadobj(r)
%             %
%             %
%             obj = Whisker.WhiskerTrial;
%             obj.trialNum = r.trialNum;
%             obj.trialType = r.trialType;
%             obj.whiskerNames = r.whiskerNames;
%             obj.trajectoryIDs = r.trajectoryIDs;
%             obj.trackerData = r.trackerData;
%             obj.faceData = r.faceData;
%             obj.framePeriodInSec = r.framePeriodInSec;
%             obj.mouseName = r.mouseName;
%             obj.sessionName = r.sessionName;
%             obj.trackerFileName = r.trackerFileName;
%             obj.useFlag = r.useFlag;
% %             obj.pxPerMm = r.pxPerMm;
% %             obj.polyFitsMask = r.polyFitsMask;
%             % For properties below, need to check if properties exist,
%             % for backwards compatability, since in early (~2008) versions
%             % files may have been saved before these properties existed.
%             % At some point can delete these if-else statements as older
%             % objects become unused.
%             if isfield(r, 'pxPerMm')
%                  obj.pxPerMm = r.pxPerMm;
%             else
%                 obj.pxPerMm = 25.7;
%             end
%             if isfield(r, 'polyFitsMask')
%                 obj.polyFitsMask = r.polyFitsMask;
%             else
%                 obj.polyFitsMask = r.polyFitsMask;
% 
%             if isfield(r,'imagePixelDimsXY')
%                 obj.imagePixelDimsXY = r.imagePixelDimsXY;
%             else
%                 obj.imagePixelDimsXY = [150 200];
%             end
%             if isfield(r,'polyFits')
%                 obj.polyFits = r.polyFits;
%             else
%                 obj.polyFits = {};
%             end
%             if isfield(r,'polyFitsROI')
%                 obj.polyFitsROI = r.polyFitsROI;
%             else
%                 obj.polyFitsROI = {};
%             end
%             if isfield(r,'maskTreatment')
%                 obj.maskTreatment = r.maskTreatment;
%             else
%                 obj.maskTreatment = 'maskNaN';
%             end
%             if isfield(r,'trackerFileFormat')
%                 obj.trackerFileFormat = r.trackerFileFormat;
%             else
%                 obj.trackerFileFormat = 'whisker0';
%             end
%             if isfield(r,'faceSideInImage')
%                 obj.faceSideInImage = r.faceSideInImage;
%             else
%                 obj.faceSideInImage = 'right';
%             end
%             if isfield(r,'barPos')
%                 obj.barPos = r.barPos;
%             else
%                 obj.barPos = [];
%             end
%             if isfield(r,'barRadius')
%                 obj.barRadius = r.barRadius;
%             else
%                 obj.barRadius = 17;
%             end
%         end
%     end
    
end
%   
end


