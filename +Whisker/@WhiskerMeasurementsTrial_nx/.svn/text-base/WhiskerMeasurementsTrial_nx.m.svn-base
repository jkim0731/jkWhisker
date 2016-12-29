classdef WhiskerMeasurementsTrial_nx < Whisker.WhiskerMeasurementsTrial
    
    properties
        touch_on = {};
        touch_off = {};
        n_touch = [];
%         wsk_num = [];
        bar_center = []; % 1x2 coordinate of bar center
    end
    
    methods (Access = public)
        function obj = WhiskerMeasurementsTrial_nx(fn, total_whisker_num)
            if nargin == 0
                fn = '';
                total_whisker_num = 4; % default maximum number of whiskers
            end
            if nargin < 2
                total_whisker_num = 4; % default maximum number of whiskers
            end
            if ~ischar(fn)
               error('Argument must be a file name string.')
            end
            
            if ~exist(fn, 'file')
                error('Input file not found.')
            end
            
%             obj.trialNum = w.trialNum;
%             obj.trialType = w.trialType;
%             obj.whiskerNames = w.whiskerNames;
%            obj.mouseName = w.mouseName;
%            obj.sessionName = w.sessionName;

            obj.trackerFileName = fn;

            % fn = 'JF8632_062308_DO79_084.measurements';
            M = Whisker.read_whisker_measurements(fn);
            
            traj = unique(M(:,1));
%             if min(abs(traj)) ~= 0
%                 traj = traj - min(traj(traj>0)); % if traj id is not starting from 0 then shift
%                 M(:,1) = M(:,1) - min(traj(traj>0));
%             end
%             obj.trajectoryIDs = sort(traj(traj >= 0 & traj < total_whisker_num));
            obj.trajectoryIDs = min(traj(traj>=0)):min(traj(traj>=0))+total_whisker_num-1;
            ntraj = length(obj.trajectoryIDs);
 
            obj.time = cell(1,ntraj);
            obj.theta = cell(1,ntraj);
            obj.kappa = cell(1, ntraj);
            obj.tracingScore = cell(1,ntraj);
            obj.whiskerLength = cell(1,ntraj);
            obj.folliclePosX = cell(1,ntraj);
            obj.folliclePosY = cell(1,ntraj);
            obj.tipPosX = cell(1,ntraj);
            obj.tipPosY = cell(1,ntraj);
            obj.distanceToBarCenter = cell(1,ntraj);
        
            for k=1:ntraj
                tid = obj.trajectoryIDs(k);
                trajInd = M(:,1)==tid;
                
                frameNum = M(trajInd,2);
                [obj.time{k},ind] = sort(frameNum * obj.framePeriodInSec);
                
                if min(abs(traj))~=0
                    temp1 = find(M(:,1)== tid+10);
                    temp2 = find(M(:,1)== tid-10);
                else
                    temp1 = find(M(:,1)== tid+10 | M(:,1)== tid+5);
                    temp2 = find(M(:,1)== tid-10 | M(:,1)== tid-5);
                end
                % find the annotated frame for touch onset 
                obj.touch_on{k} = sort(M(temp1, 2));
                % find the annotated frame for touch offset
                obj.touch_off{k} = sort(M(temp2, 2));
                
                obj.theta{k} = M(trajInd,6); obj.theta{k} = obj.theta{k}(ind);
                obj.kappa{k} = M(trajInd,7); obj.kappa{k} = obj.kappa{k}(ind);
                obj.tracingScore{k} = M(trajInd,5); obj.kappa{k} = obj.kappa{k}(ind);
                obj.whiskerLength{k} = M(trajInd,4); obj.whiskerLength{k} = obj.whiskerLength{k}(ind);
                obj.folliclePosX{k} = M(trajInd,8); obj.folliclePosX{k} = obj.folliclePosX{k}(ind);
                obj.folliclePosY{k} = M(trajInd,9); obj.folliclePosY{k} = obj.folliclePosY{k}(ind);
                obj.tipPosX{k} = M(trajInd,10); obj.tipPosX{k} = obj.tipPosX{k}(ind);
                obj.tipPosY{k} = M(trajInd,11); obj.tipPosY{k} = obj.tipPosY{k}(ind);
                if size(M,2) > 11
                    obj.distanceToBarCenter{k} = M(trajInd,12); obj.distanceToBarCenter{k} = obj.distanceToBarCenter{k}(ind);
                end
            end
        end
            
        
        function [y,t] = get_velocity(obj,tid,varargin)
            %
            %   [y,t] = get_velocity(obj,tid,varargin)
            %
            % Angular velocity in degrees per second.
            %
            % tid: Trajectory ID.
            %
            % varargin{1}: Optional smoothing window, in frames,
            %              **for position (theta) signal.** Velocity is not
            %              separately smoothed. May want to smooth theta to
            %              eliminate noise due to whisker tracking artifacts.
            %               If not specified there is no smoothing. Should be
            %               an odd number (see 'help smooth').
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            %             t = obj.time{ind};
            if nargin > 2
                [theta,t] = obj.get_position(tid,varargin{1});
            else
                [theta,t] = obj.get_position(tid);
            end
%             y = [0 diff(theta)] ./ [0 diff(t)]; % in degrees/sec
              y = [0; diff(theta)] ./ [0; diff(t)]; % in degrees/sec  
        end
        
        function [y,t] = get_velocity_medfilt(obj,tid,varargin)
            %
            %   [y,t] = get_velocity_medfilt(obj,tid,varargin)
            %
            % Angular velocity in degrees per second, after filtering
            % position signal with a median filter.
            %
            % tid: Trajectory ID.
            %
            % varargin{1}: Optional smoothing window, in frames,
            %              **for position (theta) signal.** Velocity is not
            %              separately filtered. May want to filter theta to
            %              eliminate noise due to whisker tracking artifacts.
            %              Default is 3. Should be
            %               an odd number (see help medfilt1).
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            %             t = obj.time{ind};
            if nargin > 2
                span = varargin{1};
                if mod(span,2)==0
                    disp('Varargin{1}, should be odd; will be rounded down.')
                end
            else
                span = 3;
            end
            
            [theta,t] = obj.get_position(tid);
            
            % Check if frames are evenly spaced. If not, must interpolate theta for
            % missing frames prior to applying filter. Will do that, then
            % take only smoothed theta values at the original (non-interpolated)
            % frames.
            frames = t ./ obj.framePeriodInSec;
            if length(unique(diff(frames))) > 1
                newframes = min(frames):max(frames);
                newy = interp1(frames,theta,newframes,'linear');
                yy = medfilt1(newy,span)';
                y = interp1(newframes,yy,frames,'linear'); % could use ismember instead
            else
                y = medfilt1(theta,span);
            end
            % Sanity check:
            if length(y) ~= length(t)
                error('y and t are of unequal lengths.')
            end
            
            y = [0; diff(y)] ./ [0; diff(t)]; % in degrees/sec
        end
        
        function [y,t] = get_acceleration(obj,tid,varargin)
            %
            %   [y,t] = get_acceleration(obj,tid,varargin)
            %
            % Angular acceleration in degrees per second^2.
            %
            % tid: Trajectory ID.
            %
            % varargin{1}: Optional moving average smoothing window, in frames,
            %              **for position (theta) signal.** Acceleration is not
            %              separately smoothed. May want to smooth theta to
            %              eliminate noise due to whisker tracking artifacts.
            %               If not specified there is no smoothing. Should be
            %               an odd number (see 'help smooth').
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            %             t = obj.time{ind};
            
            if nargin > 2
                [velocity,t] = obj.get_velocity(tid,varargin{1});
            else
                [velocity,t] = obj.get_velocity(tid);
            end
            y = [0; diff(velocity)] ./ [0; diff(t)]; % in degrees/sec
        end
        
        
    end
    
end