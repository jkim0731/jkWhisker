classdef WhiskerTrialLite < handle
    %
    %
    %   WhiskerTrialLite < handle
    %
    %
    % Contains just final timeseries measurements from whiskers (kappa, theta, etc)
    % and meta-information in order to match with other data (ephys, imaging, etc).
    %
    %
    % DHO, 3/10.
    %
    properties
        trialNum = [];
        trialType = NaN;
        whiskerNames = {};  % whiskerNames and trajectoryIDs must be of same length
        trajectoryIDs = []; % with matching elements.
        framePeriodInSec = 0.002;
        pxPerMm = []; %  Inherited from WhiskerSignalTrial.
        mouseName = '';
        sessionName = '';
        trackerFileName = '';
        useFlag = 1;
        
        rInMm = [];
        
        time = {};
        deltaKappa = {};  
        thetaAtBase = {}; 
        follicleCoordsX = {}; 
        follicleCoordsY = {}; 
        
        % The following require pole tracking:
        thetaAtContact = {};
        distanceToPoleCenter = {};
        Fnorm = {};
        Faxial = {};
        Flateral = {};
        M0 = {}; % Moment at the follicle.
        meanKappa = {};
        
        barPos = []; %  Inherited from WhiskerSignalTrial. [frameNum XPosition YPosition]
        barPosOffset = []; % Inherited from WhiskerSignalTrial. [x y], either 1X2 or nframesX2
        barRadius = []; % Inherited from WhiskerSignalTrial.  In pixels. Must be radius of bar tracked by the bar tracker.

      
    end
    
    properties (Dependent = true)
        theta % For compatability, make this 'alias' to refer to thetaAtBase
        kappa % For compatability, make this 'alias' to refer to deltaKappa
    end
    
        
    methods (Access = public)
        function obj = WhiskerTrialLite(w,varargin)
            %
            %
            %
            %   obj = WhiskerTrialLite(w)
            %   obj = WhiskerTrialLite(w, 'Parameter',value, etc) (see below for parameter/value pairs).
            %
            % INPUTS:
            %   w: a WhiskerSignalTrial object.
            %
            %   Optional parameter/value pair arguments:
            %   
            %   'r_in_mm': The arc-length along whisker at which to measure kappa. Units of mm. Defaults to 1 mm.
            %   
            %   'calc_forces': Either true or false. Requires the pole position
            %           to be tracked (i.e., barPos property of WhiskerSignalTrial must
            %           not be empty). Default is false. If true, will calculate the following timeseries:
            %               -M0:  Moment at the follicle. In Newton-meters.
            %               -Faxial: Axial force into follice. In Newtons.
            %               -Flateral: Lateral force, orthogonal to Faxial. In Newtons.
            %               -deltaKappa: Change from baseline curvature, at point specified by r_point. In 1/mm.
            %               -Fnorm: The force on the whisker normal to the contacted object. In Newtons.
            %               -thetaAtBase: The whisker angle nearest the follicle. In degrees.
            %               -thetaAtContact: The whisker angle nearest the point of contact. I.e., nearest the center of the pole. In degrees.
            %               -distanceToPoleCenter: The closest distance between the whisker and the center of the pole. In mm.
            %               -meanKappa: The mean of kappa over the entire secondary polynomial fitted ROI. In 1/mm.
            %
            %   The following optional parameter/value pair arguments are ignored if 'calc_forces'
            %   is not true:
            %
            % 	'whisker_radius_at_base': Given in microns. Defaults is 33.5 microns.
            %
            % 	'whisker_length': Given in mm. Default is 16 mm.
            %
            % 	'youngs_modulus': In Pa. Default is 5e9 Pa.
            %
            %   'baseline_time_or_kappa_value': Either (1) a 1x2 vector giving starting and stopping times (inclusive) for measuring baseline whisker curvature, in sec; 
            %                                   or (2) a scaler giving a baseline kappa value (measured by the user separately) to directly subtract from kappa
            %                                   timeseries, in 1/mm. Default is [0 0.1].
            %
            %   'proximity_threshold':  Optional argument giving distance from nearest point on whisker
            %                           to bar center, in units of bar radius, beyond which
            %                           the whisker will be extrapolated along the last theta in
            %                           order to determine distance between whisker and bar. Default is not to use.
            %
            % NOTES:
            %   Still need make these arguments settable on a whisker-by-whisker (trajectory ID by trajectory ID) basis.   
            %
            %
            %
            if nargin==0
                return
            end
            
            p = inputParser;
            p.addRequired('w', @(x) isa(x,'Whisker.WhiskerSignalTrial'));                      
            p.addParamValue('r_in_mm', 1, @(x) isnumeric(x) && numel(x)==1);
            p.addParamValue('calc_forces', false, @islogical);
       
            p.addParamValue('whisker_radius_at_base', 33.5, @isnumeric);
            p.addParamValue('whisker_length', 16, @isnumeric);
            p.addParamValue('youngs_modulus', 5e9, @isnumeric);
            p.addParamValue('baseline_time_or_kappa_value', [0 0.1], @isnumeric);
            p.addParamValue('proximity_threshold', -1, @isnumeric);
            
            p.parse(w,varargin{:});
           
              
            obj.trialNum = w.trialNum;
            obj.trialType = w.trialType;
            obj.whiskerNames = w.whiskerNames;
            obj.trajectoryIDs = w.trajectoryIDs;
            obj.framePeriodInSec = w.framePeriodInSec;
            obj.pxPerMm = w.pxPerMm;
            obj.mouseName = w.mouseName;
            obj.sessionName = w.sessionName;
            obj.trackerFileName = w.trackerFileName;
            obj.follicleCoordsX = w.follicleCoordsX;
            obj.follicleCoordsY = w.follicleCoordsY;  
            obj.useFlag = w.useFlag;
            obj.rInMm = p.Results.r_in_mm;
            
            ntraj = length(obj.trajectoryIDs);
                        
            obj.time = cell(1,ntraj);
            obj.deltaKappa = cell(1,ntraj);  
            obj.thetaAtBase = cell(1,ntraj);     
            obj.thetaAtContact = cell(1,ntraj);
            obj.distanceToPoleCenter = cell(1,ntraj);
            obj.Fnorm = cell(1,ntraj);
            obj.Faxial = cell(1,ntraj);
            obj.Flateral = cell(1,ntraj);
            obj.M0 = cell(1,ntraj); % Moment at the follicle.
            obj.meanKappa = cell(1,ntraj); % Mean kappa over the ROI.
          
            obj.barPos = w.barPos; %  Inherited from WhiskerSignalTrial. [frameNum XPosition YPosition]
            obj.barPosOffset = w.barPosOffset; % Inherited from WhiskerSignalTrial. [x y], either 1X2 or nframesX2
            obj.barRadius = w.barRadius; % Inherited from WhiskerSignalTrial.  In pixels. Must be radius of bar tracked by the bar tracker.
        
            for k=1:ntraj
                tid = obj.trajectoryIDs(k);
                
                obj.time{k} = w.get_time(tid);
                
                if all(isnan(obj.time{k})) % For instance with TIDs used for contact annotation, may not have any observations.
                    obj.deltaKappa{k} = NaN;
                    obj.thetaAtBase{k} = NaN;
                    if p.Results.calc_forces == true
                        obj.thetaAtContact{k} = NaN;
                        obj.distanceToPoleCenter{k} = NaN;
                        obj.M0{k} = NaN;
                        obj.Faxial{k} = NaN;
                        obj.Flateral{k} = NaN;
                        obj.Fnorm{k} = NaN;
                        obj.meanKappa{k} = NaN;
                    end
                    continue
                end
               
                if p.Results.calc_forces == true
                   [obj.M0{k},obj.Faxial{k},t,obj.deltaKappa{k},obj.Fnorm{k},...
                       obj.thetaAtBase{k},obj.thetaAtContact{k},obj.distanceToPoleCenter{k}, obj.meanKappa{k}, obj.Flateral{k}] = ...
                       w.calc_M0_Faxial(tid,p.Results.r_in_mm,p.Results.whisker_radius_at_base, p.Results.whisker_length,...
                       p.Results.youngs_modulus,p.Results.baseline_time_or_kappa_value,p.Results.proximity_threshold); 
                   
                else
                    % Should consolidate into single function to optimize the following: 
                    [theta,obj.deltaKappa{k},y,x,t] = w.get_theta_kappa_at_roi_point(tid,p.Results.r_in_mm);
                    [obj.thetaAtBase{k},kappa0,t] = w.get_theta_kappa_at_base(tid);
                end
                
            end
        end
        
        function tid = name2tid(obj, whisker_name)
            if ~ischar(whisker_name)
                error('Argument whisker_name must be a string.')
            end
            if numel(obj.whiskerNames) ~= numel(obj.trajectoryIDs)
                error('This WhiskerTrialLite does not have matching whiskerNames and trajectoryIDs.')
            end
            tid = obj.trajectoryIDs( strmatch(whisker_name, obj.whiskerNames) );
        end
        
        function whisker_name = tid2name(obj, trajectory_id)
            if ~isnumeric(trajectory_id)
                error('Argument trajectory_id must be an integer.')
            end
            if numel(trajectory_id) > 1
                error('Only one trajectory_id is allowed.')
            end
            whisker_name = obj.whiskerNames(obj.trajectoryIDs==trajectory_id);
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
            T = cell2mat(obj.time); Y = cell2mat(obj.theta);
            r = Shared.tapply([T' Y']);
            t = r(:,1); y = r(:,2);
        end
              
        function [y,t] = get_follicle_translation(obj,tid)
            %
            %    [y,t] = get_follicle_translation(obj,tid)
            %
            %  INPUTS:
            %
            %   tid: Trajectory ID.
            %
            %
            %  RETURNS:
            %
            % y: The distance the follicle has translated from the previous
            %    frame. Units of pixels.
            %
            % t: The time of each observation in y.
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Could not find specified trajectory ID.')
            end
            
            if isempty(obj.follicleCoordsX) || isempty(obj.follicleCoordsX)
                error(['obj.follicleCoordsX or obj.follicleCoordsY is empty. ' ...
                    'Must run obj.recompute_cached_follicle_coords before this method.'])
            end
            
            t = obj.get_time(tid);
            
            dx = [0 diff(obj.follicleCoordsX{ind})];
            dy = [0 diff(obj.follicleCoordsY{ind})];
            
            y = sqrt(dx.^2 + dy.^2);
        end
        
        function [y,x,t] = get_cached_follicle_coords(obj,tid)
            %
            %    [y,x,t] = get_cached_follicle_coords(obj,tid)
            %
            %  INPUTS:
            %
            %   tid: Trajectory ID.
            %
            %
            %  RETURNS:
            %
            % y: The y coordinate in image pixels of the follicle for each
            %    time point (i.e. frame).
            %
            % x: The x coordinate in image pixels of the follicle for each
            %    time point.
            %
            % t: The time of each observation in x,y.
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Could not find specified trajectory ID.')
            end
            
            if isempty(obj.follicleCoordsX) || isempty(obj.follicleCoordsX)
                error(['obj.follicleCoordsX or obj.follicleCoordsY is empty. ' ...
                    'Must run obj.recompute_cached_follicle_coords in WhiskerSignalTrial object before '...
                    'creating this WhiskerTrialLite.'])
            end
            
            t = obj.get_time(tid);
            
            x = obj.follicleCoordsX{ind};
            y = obj.follicleCoordsY{ind};
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

        function [y,t] = get_follicleCoordsX(obj,tid)
            %
            %   [y,t] = get_follicleCoordsX(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            if isempty(obj.follicleCoordsX)
                y = nan(size(t));
            else
                y = obj.follicleCoordsX{ind};
            end
        end
        
        function [y,t] = get_follicleCoordsY(obj,tid)
            %
            %   [y,t] = get_follicleCoordsY(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            if isempty(obj.follicleCoordsY)
                y = nan(size(t));
            else
                y = obj.follicleCoordsY{ind};
            end
        end
        
        function [y,t] = get_thetaAtContact(obj,tid)
            %
            %   [y,t] = get_thetaAtContact(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.thetaAtContact{ind};
        end

        function [y,t] = get_distanceToPoleCenter(obj,tid)
            %
            %   [y,t] = get_distanceToPoleCenter(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.distanceToPoleCenter{ind};
        end
        
        function [y,t] = get_thetaAtBase(obj,tid)
            %
            %   [y,t] = get_thetaAtBase(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.thetaAtBase{ind};
        end
        
        function [y,t] = get_deltaKappa(obj,tid)
            %
            %   [y,t] = get_deltaKappa(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.deltaKappa{ind};
        end
        
        function [y,t] = get_Fnorm(obj,tid)
            %
            %   [y,t] = get_Fnorm(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.Fnorm{ind};
        end
        
        function [y,t] = get_Faxial(obj,tid)
            %
            %   [y,t] = get_Faxial(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.Faxial{ind};
        end
        
        function [y,t] = get_Flateral(obj,tid)
            %
            %   [y,t] = get_Flateral(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.Flateral{ind};
        end
        
        function [y,t] = get_M0(obj,tid)
            %
            %   [y,t] = get_M0(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.M0{ind};
        end
        
        function [y,t] = get_meanKappa(obj,tid)
            %
            %   [y,t] = get_meanKappa(obj,tid)
            %
            %
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            t = obj.time{ind};
            y = obj.meanKappa{ind};
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
            y = [0 diff(theta)] ./ [0 diff(t)]; % in degrees/sec
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
            y = [0 diff(velocity)] ./ [0 diff(t)]; % in degrees/sec
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
            y = obj.kappa{ind};% * obj.pxPerMm;
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
        
        function t = get_time(obj,tid)
            %
            %   t = get_time(obj,tid)
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %
            %   OUTPUTS:
            %       t: time in seconds for each sample in this WhiskerTrialLite
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
    
    methods % Dependent property methods; cannot have attributes.
        
        function value = get.theta(obj)
            value = obj.thetaAtBase;
        end
        
        function value = get.kappa(obj)
            value = obj.deltaKappa;
        end
        
    end
       
end
