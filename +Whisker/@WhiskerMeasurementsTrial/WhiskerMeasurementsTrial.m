classdef WhiskerMeasurementsTrial < handle
    %
    %
    %   WhiskerMeasurementsTrial < handle
    %
    %
    % DHO, 11/09.
    %
    %
    properties
        trialNum = [];
        trialType = NaN;
        whiskerNames = {};  % whiskerNames and trajectoryIDs must be of same length
        trajectoryIDs = []; % with matching elements.
        framePeriodInSec = 0.002;
        mouseName = '';
        sessionName = '';
        trackerFileName = '';

        kappa = {}; % mean curvature
        theta = {}; % angle at follicle
        time = {};
        tracingScore = {};
        whiskerLength = {};
        folliclePosX = {};
        folliclePosY = {};
        tipPosX = {};
        tipPosY = {};
        distanceToBarCenter = {};
        
        useFlag = 1;
    end
    
    methods (Access = public)
        function obj = WhiskerMeasurementsTrial(fn)
            %
            %
            %
            %   obj = WhiskerMeasurementsTrial(fn)
            %
            %   fn: file name of a .measurements file.
            %
            
            if nargin==0
                return
            end
            
            if nargin > 2
                error('Too many input arguments.')
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
            obj.trajectoryIDs = sort(traj(traj >= 0));
            
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
        
        
        function r = curvatureDot(obj,varargin)
            %
            %   r = curvatureDot(obj,varargin)
            %
            % varargin{1}: vector of trajectory IDs.
            %
            % If only a single trajectory is specified, r is a vector.
            % If multiple trajectories are specified, r is a cell array
            % of vectors.
            %
            if nargin > 1
                tid = varargin{1};
            else
                tid = obj.trajectoryIDs;
            end
            if isempty(tid)
                tid = obj.trajectoryIDs;
            end
            
            ntraj = length(tid);
            if ntraj > 1
                r = cell(1,ntraj);
                for k=1:ntraj
                    ind = obj.trajectoryIDs==tid(k);
                    if max(ind) < 1
                        error('Trajectory ID was not found.')
                    end
                    t = obj.time{ind};
                    r{k} = [0 diff(obj.kappa{ind})] ./ [0 diff(t)];
                end
            else
                ind = obj.trajectoryIDs==tid;
                if max(ind) < 1
                    error('Trajectory ID was not found.')
                end
                t = obj.time{ind};
                r = [0 diff(obj.kappa{ind})] ./ [0 diff(t)];
            end
        end
        
        function [y,t] = get_position(obj,tid,varargin)
            %
            %   [y,t] = get_position(obj,tid,varargin)
            %
            % tid: Trajectory ID.
            %
            % varargin{1}: Optional smoothing span, in frames.
            %               If not specified there is no smoothing. Should be
            %               an odd number (see 'help smooth').
            %
            %
            % Returns angle with respect to y-axis of coordinate plane
            % in which slope is computed (i.e., of the high-speed video
            % image). Whisker angle is with respect to mouse midline if
            % y-axis of image is parallel with mouse midline.
            %
            % If a *decrease* in slope in image coordinates cooresponds to
            % whisker protraction, then this function returns increasing
            % angles during protraction.
            %
            % Angle is 0 deg when whisker is perpendicular to midline,
            % increases to +90 deg with protraction all the way to parallel
            % to the midline, and decreases to -90 deg with retraction all
            % the way to the midline.
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.theta{ind};
            
            if nargin > 2
                span = varargin{1};
                if span < 3
                    error('Smoothing window in varargin{2} should be an odd number > 3.')
                end
                if mod(span,2)==0
                    disp('Smoothing window in varargin{2} should be an odd number. Smooth() will round it.')
                end
                % Check if frames are evenly spaced. If not, must interpolate theta for
                % missing frames prior to applying a moving average. Will do that, then
                % take only smoothed theta values at the original (non-interpolated)
                % frames.
                frames = t ./ obj.framePeriodInSec;
                if length(unique(diff(frames))) > 1
                    newframes = min(frames):max(frames);
                    newy = interp1(frames,y,newframes,'linear');
                    yy = smooth(newy,span,'moving')';
                    y = interp1(newframes,yy,frames,'linear'); % could use ismember instead
                else
                    y = smooth(y,span,'moving')';
                end
                % Sanity check:
                if length(y) ~= length(t)
                    error('y and t are of unequal lengths.')
                end
                
            end
            
        end
        
        function [y,t] = get_mean_position(obj)
            %
            %   [y,t] = get_mean_position(obj)
            %
            %   y is the ***mean position of all whiskers*** (tids) in trial.
            %   This is for use, e.g., in determining whether there is overall
            %   whisking after fully-automated tracking.            
            %
            % varargin{1}: Optional smoothing span, in frames.
            %               If not specified there is no smoothing. Should be
            %               an odd number (see 'help smooth').
            %
            %
            % Returns angle of average whisker with respect to y-axis of coordinate plane
            % in which slope is computed (i.e., of the high-speed video
            % image). Whisker angle is with respect to mouse midline if
            % y-axis of image is parallel with mouse midline.
            %
            % If a *decrease* in slope in image coordinates cooresponds to
            % whisker protraction, then this function returns increasing
            % angles during protraction.
            %
            % Angle is 0 deg when whisker is perpendicular to midline,
            % increases to +90 deg with protraction all the way to parallel
            % to the midline, and decreases to -90 deg with retraction all
            % the way to the midline.
            %
            %
            T = cell2mat(obj.time'); Y = cell2mat(obj.theta');
            r = Shared.tapply([T Y]);
            t = r(:,1); y = r(:,2);
        end 
        
        function [y,t] = get_curvature(obj,tid)
            %
            %   [y,t] = get_curvature(obj,tid)
            %
            %   t: time in seconds.
            %   y: whisker curvature in units of pixels^(-1).
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.kappa{ind};
        end
        
        function [y,t] = get_curvatureChange(obj,tid,varargin)
            %
            %   [y,t] = get_curvatureChange(obj,tid,varargin)
            %
            %   varargin{1}: Optional period to use for computing
            %       baseline curvature. Starts at 0 but user
            %       specifies endpoint in seconds. Default is 0.05 s.
            %
            %   t: time in seconds.
            %   y: whisker curvature change in units of pixels^(-1).
            %
            
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.kappa{ind};
            
            if nargin > 2
                baselinePeriodEnd = varargin{1};
            else
                baselinePeriodEnd = 0.05;
            end
            
            baseline = mean(y(t<=baselinePeriodEnd));
            y = y-baseline;
        end
        
        function [y,t] = get_curvatureDot(obj,tid)
            %
            %   [y,t] = get_curvatureDot(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.curvatureDot(tid);
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
%             y = [0; diff(y)] ./ [0; diff(t)]; % in degrees/sec 
            y = [0 diff(y)] ./ [0 diff(t)]; % in degrees/sec
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
        
        function plot_whisker_angle(obj,tid,varargin)
            %
            %   plot_whisker_angle(obj,tid,varargin)
            %
            % tid: A single trajectory ID.
            %
            % varargin{1}: Optional plot color/symbol string specifier.
            %              Can be empty ([]) to allow access to varargin{2}.
            %
            % varargin{2}: Optional moving average smoothing window, in frames.
            %               If not specified there is no smoothing. Should be
            %               an odd number (see 'help smooth').
            %
            % If a *decrease* in slope in image coordinates cooresponds to
            % whisker protraction, then this function plots increasing
            % angles during protraction.
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            
            if nargin==2
                plotString = 'k.-';
                [y,t] = obj.get_position(tid);
            elseif nargin==3
                plotString = varargin{1};
                if isempty(plotString)
                    plotString = 'k.-';
                end
                [y,t] = obj.get_position(tid);
            elseif nargin==4
                plotString = varargin{1};
                if isempty(plotString)
                    plotString = 'k.-';
                end
                span = varargin{2};
                [y,t] = obj.get_position(tid,span);
            end
            
            plot(t,y,plotString);
            
            set(gca,'TickDir','out','box','off')
            xlabel('Sec')
        end
        
        function plot_whisker_curvature(obj,tid,varargin)
            %
            %   plot_whisker_curvature(obj,tid,varargin)
            %
            % varargin{1}: Optional plot color/symbol string specifier.
            %
            %   Plots whisker curvature (units of pixels^(-1)) against
            %   time (in seconds).
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.kappa{ind};
            if nargin > 2
                plot(t,y, varargin{1});
            else
                plot(t,y, 'k.-');
            end
            set(gca, 'TickDir','out','box','off')
            xlabel('Sec')
        end
        
        function plot_whisker_curvatureDot(obj,tid,varargin)
            %
            %   plot_whisker_curvatureDot(obj,tid,varargin)
            %
            % varargin{1}: Optional plot color/symbol string specifier.
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.curvatureDot(tid);
            if nargin > 2
                plot(t,y, varargin{1});
            else
                plot(t,y, 'k.-');
            end
        end
        
     
        
        function [t,theta,kappa] = mean_theta_and_kappa(obj,tid,radial_window_theta,radial_window_kappa)
            %
            %  [t,theta,kappa] = mean_theta_and_kappa(obj,tid,radial_window_theta,radial_window_kappa)
            %
            %   The whisker is parameterized as c(q) = (x(q),y(q)), where q has length(x)
            %   and is in [0,1].
            %
            % INPUTS:
            %
            %   tid: Whisker trajectory ID.
            %
            %   radial_window_kappa: 2x1 vector giving arc length region of whisker to
            %   average over for mean kappa measurment, in format [startDistance stopDistance].
            %   Values are inclusive and in units of pixels. If empty ([]), averages over
            %   the whole whisker.
            %
            %   radial_window_theta: 2x1 vector giving arc length region of whisker to
            %   average over for mean kappa measurment, in format [startDistance stopDistance].
            %   Values are inclusive and in units of pixels. If empty ([]), averages over
            %   the whole whisker, which is not likely useful for theta.
            %
            %   For this function to work properly, object property 'radialDirection'
            %   must be correctly set.  See help Whisker.WhiskerTrial.set_radialDirection.
            %
            % RETURNS:
            %
            %   t:  Time in seconds corresponding to each frame.
            %
            %   theta: Angle of the line tangent to the whisker (i.e., to c(q)), averaged over
            %           radial distances (arc length) between and including radial_window_theta(1)
            %           and radial_window_theta(2).
            %
            %   kappa:  Signed curvature averaged over radials distances (arc length) between and including
            %           radial_window_kappa(1) and radial_window_kappa(2), for each frame.
            %           Units of 1/pixels. Abs(kappa(q)) is 1/X where X is
            %           the radius in pixels of the osculating circle at c(q).
            %
            %
            %
            
            [R,THETA,KAPPA] = obj.arc_length_theta_and_kappa(tid);
            
            if isempty(radial_window_kappa)
                kappa = cellfun(@mean, KAPPA);
            else
                %                 kappa = cellfun(@(x,y) mean(x(y >= radial_window_kappa(1) & y <= radial_window_kappa(2))),KAPPA,R);
                nframes = length(R);
                kappa = nan(1,nframes);
                for k=1:length(R)
                    if ~isnan(R{k}(1))
                        kappa(k) = mean(KAPPA{k}( R{k} >= radial_window_kappa(1) & R{k} <= radial_window_kappa(2) ));
                    end
                end
            end
            
            
            if isempty(radial_window_theta)
                theta = cellfun(@mean, THETA);
            else
                %                 theta = cellfun(@(x,y) mean(x(y >= radial_window_theta(1) & y <= radial_window_theta(2))),THETA,R);
                nframes = length(R);
                theta = nan(1,nframes);
                for k=1:length(R)
                    if ~isnan(R{k}(1))
                        theta(k) = mean(THETA{k}( R{k} >= radial_window_theta(1) & R{k} <= radial_window_theta(2) ));
                    end
                end
            end
            
            t = obj.get_time(tid);
            
            % Interpolate to fill missing (NaN) values, arising for instance if there
            % weren't enough pixels to do curve fitting. NaN also arise if obj.maskTreatment
            % is set to 'maskNaN'.
            missing = isnan(kappa);
            if sum(missing)==length(kappa)
                error('Less than 2 non-NaN values for mean kappa timeseries. May be a bad mask, that no whiskers cross.')
            end
            kappa = interp1(t(~missing),kappa(~missing),t,'linear','extrap');
            
            missing = isnan(theta);
            if sum(missing)==length(theta)
                error('Less than 2 non-NaN values for mean theta timeseries. May be a bad mask, that no whiskers cross.')
            end
            theta = interp1(t(~missing),theta(~missing),t,'linear','extrap');
        end
        
        function t = get_time(obj,tid)
            %
            %   t = get_time(obj,tid)
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %
            %   OUTPUTS:
            %       t: time in seconds for each sample in this WhiskerSignalTrial
            %       for given trajectory ID or whisker name.
            %
            if isnumeric(tid) % Trajectory ID specified.
                ind = find(obj.trajectoryIDs == tid);
            elseif ischar(tid) % Whisker name specified.
                ind = strmatch(tid,obj.whiskerNames,'exact');
            else
                error('Invalid type for argument ''tid''.')
            end
            
            if isempty(ind)
                error('Could not find specified trajectory ID.')
            end
            t = obj.time{ind};
        end

        
    end
    
end
