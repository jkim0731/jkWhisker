classdef WhiskerSignalTrial < handle
    %
    %
    %   WhiskerSignalTrial < handle
    %
    % A class that deals with polynomial fits to tracked whisker (x,y) data.
    % The original (x,y) data are not preserved in WhiskerSignalTrial.
    %
    % DHO, 8/08.
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
        
        % polyFits: Inherited from WhiskerTrial.
        % Cell array of length length(trajectoryIDs), of format:
        % {{XPolyCoeffs_tid0, YPolyCoeffs_tid0},...,{XPolyCoeffs_tidN, YPolyCoeffs_tidN}};
        polyFits = {};
        
        polyFitsROI = {}; % polyFitsROI: Inherited from WhiskerTrial.
                          % Same format as polyFits but polynomials are fitted only to a 
                          % constant region of arc length, and in addition to x and y coefficients
                          % there is stored the "q" values, i.e. the points along the normalized full
                          % whisker ([0,1]) that the ROI fitting begins,as well as the two corresponding
                          % values in units of pixels.  We store both for speed later. If a mask will be specified to define
                          % the arc-length origin it must be applied prior to populating polyFitsROI.
                          % Populated by method fit_polys_roi().

        % polyFitsMask:
        % Cell array of length length(trajectoryIDs), of format:
        %
        % There is generally a noisy edge to the tracked
        % whiskers on the side of whisker pad, which can interfere
        % with proper measurement of radial distances. For each whisker
        % for each trial, can specify here a polynomial in image coordinate
        % space to "mask out" the noisy edge. That is, radial distance for
        % purposes of mean theta and mean kappa measurements will be measured starting
        % at the intersection of the tracked whisker with this masking polynomial
        % if the whisker in fact crosses the masking polynomial. I.e., the radial
        % distance is r_new = r - r_intersection where r_new is the
        % new radial distance used in mean theta and mean kappa measurements, r is the
        % original radial distance, and r_intersection is the point of intersection
        % between the fitted whisker and the masking polynomial. If there is no intersection,
        % then r_new = r.  Also, if polyFitsMask is empty (or is empty for a given whisker)
        % then r_new = r.
        %
        % Ultimately, may want to do this separately for every frame, perhaps after
        % face tracking.
        %
        % Polynomials can be of any order, and are reconstructed based on the number
        % of coefficients.
        % If polyFitsMask{k}{1} and polyFitsMask{k}{2} are NxM matrices where
        % M is the polynomial degree + 1 and N is the number of frames, then each frame
        % has its own mask. For instance, this could be used after face tracking.
        % If instead polyFitsMask{k}{1} and polyFitsMask{k}{2} are 1xM vectors
        % where M is the polynomial degree + 1, then the same mask is used for all
        % frames.
        %         polyFitsMask = {{[25 120],[157.5 31.5]},{[25 120],[157.5 31.5]},{[25 120],[157.5 31.5]}};
        polyFitsMask = {};
        
        %   maskTreatment: String describing treatment of mask. Or, can be cell array
        %                   of strings, of same length as obj.trajectoryIDs and with
        %                   matching entries, in order to set maskTreatment differently
        %                   for different trajectory IDs.
        %               Values: 'none', 'mask', 'maskNaN'.
        %                   none: Ignore the mask.
        %                   mask: Subtract from each radial distance in R
        %                      the radial distance at the intersection of
        %                      each fitted whisker with the mask.  If there is
        %                      no intersection for a given whisker, make no change
        %                      in the radial distance measurement: i.e. 0 is still
        %                      at the end. If obj.polyFitsMask is empty ({}), make
        %                      no change.
        %                   maskNaN: Same as mask except that if the whisker does
        %                      not intersect the mask in a given frame, set all its
        %                      values in R to NaN.
        maskTreatment = 'maskNaN'; % Alternatively, {'maskNaN','none','mask'} to set three (or more)
        % different trajectories to have three (or more) different treatments.
        
        kappa = {}; % Can always recompute from polyFits and q, but cache for speed.
        theta = {}; % Can always recompute from polyFits and q, but cache for speed.
        follicleExtrapDistInPix = 0; % Sets extrapolated distance past end of whisker, or past intersection of whisker 
                                     % and mask if a mask is defined, at which follicle is estimated to be. Is not used until
                                     % recompute_cached_follicle_coords() is called.  
        follicleCoordsX = {}; % Can always recompute, but cache for speed.
        follicleCoordsY = {};
        barPos = []; %  Inherited from WhiskerTrial. [frameNum XPosition YPosition]
        barPosOffset = []; % [x y], either 1X2 or nframesX2
        barRadius = []; % Inherited from WhiskerTrial.  In pixels. Must be radius of bar tracked by the bar tracker.
        time = {};
        pxPerMm = 22.68; %  Inherited from WhiskerTrial, but give default value.
        faceSideInImage = 'top'; % Inherited from WhiskerTrial, but give default value.
        % Can be: 'top', 'bottom', 'left','right'.
        % May need to make this a cell array of strings, one per trajectory ID
        % to handle case of tracked whiskers on both sides of head.
        protractionDirection = 'rightward';  % Inherited from WhiskerTrial, but give default value.
                                            % Can be: 'downward','upward','rightward','leftward'.

        useFlag = 1;
    end
    
    properties (Dependent = true)
        barPosClean % Bar position after processing to clean up bar tracker errors/limitations.
    end
    
    methods (Access = public)
        function obj = WhiskerSignalTrial(w, varargin)
            %
            % USAGE:
            %
            %   obj = WhiskerSignalTrial(w)
            %   obj = WhiskerSignalTrial(w, 'polyRoiInPix',[roiMin, roiMax])
            %   obj = WhiskerSignalTrial(w, 'polyRoiInPix',{[trajectoryIDs],[roiMin1, roiMax1], ...
            %                                     [roiMin2, roiMax2], ...
            %                                     [roiMin3, roiMax3], ...
            %                                     [roiMinN, roiMaxN]})
            %
            % INPUTS:
            %   w: a WhiskerTrial object.
            %
            %   Optional argument polyRoiInPix:
            %       Sets arc-length limits (in pixels) on which to perform secondary curve fitting.
            %       This argument can be given in two forms:
            %           (1) an 1x2 vector that gives the ROI for *all* whiskers; or
            %           (2) a cell array where first element is a vector of trajectory IDs
            %               (of length N) and subsequent elements comprise N 1x2 vectors
            %               giving ROIs for the trajectory IDs specified in the first
            %               element (respectively).
            %       Limits are inclusive.
            %   
            %   If WhiskerSignalTrial(w) is called without 'polyRoiInPix'
            %   argument, and w.polyFitsROI is empty (i.e.,
            %   fit_polys_roi() method of WhiskerTrial w was not called earlier), then
            %   an error is given.  A WhiskerSignalTrial requires
            %   polynomials to be fitted to the whiskers.
            %
            %
            p = inputParser;
            p.addOptional('w', @(x) isa(x,'Whisker.WhiskerTrial'));                      
            p.addParamValue('polyRoiInPix', NaN); 
            p.parse(varargin{:});
            
            if nargin==0
                return
            end
            
            if isempty(w.polyFitsROI) && any(isnan(p.Results.polyRoiInPix))
                error(['If ''fit_polys_roi()'' has not already been called on ' ...
                    'WhiskerSignalTrial argument, then argument ''polyRoiInPix'' must be given.'])
            end

            if ~iscell(p.Results.polyRoiInPix) % Single ROI given; not specified individually for different trajectories.
                tidList = w.trajectoryIDs;
                roiAll = cell(1,length(tidList));
                % Copy ROI for each trajectory:
                for k=1:length(roiAll)
                    roiAll{k} = p.Results.polyRoiInPix;
                end
            else
                tidList = p.Results.polyRoiInPix{1};
                roiAll = p.Results.polyRoiInPix(2:end);
            end
            
            if length(tidList) ~= length(roiAll)
                error(['When given as a cell array, the first element of argument ''polyRoiInPix must be 1xN vector of trajectory IDs' ...
                    'and subsequent elements must comprise N 1x2 vectors giving the ROI for each trajectory.']);
            end
            
            if length(tidList) ~= length(w.trajectoryIDs)
                error('Number of trajectory IDs specified in varagin{1} does not match number of trajectory IDs in WhiskerTrial argument w.')
            end
            
            obj.trialNum = w.trialNum;
            obj.trialType = w.trialType;
            obj.whiskerNames = w.whiskerNames;
            obj.trajectoryIDs = w.trajectoryIDs;
            obj.framePeriodInSec = w.framePeriodInSec;
            obj.mouseName = w.mouseName;
            obj.sessionName = w.sessionName;
            obj.trackerFileName = w.trackerFileName;
            obj.faceSideInImage = w.faceSideInImage;
            obj.protractionDirection = w.protractionDirection;
            obj.pxPerMm = w.pxPerMm;
            obj.barPos = w.barPos;
            obj.barRadius = w.barRadius;
            obj.barPosOffset = w.barPosOffset;  
            obj.polyFitsMask = w.polyFitsMask;
            
            ntraj = length(obj.trajectoryIDs);
            
            obj.theta = cell(1,ntraj);
            obj.kappa = cell(1,ntraj);
            obj.time = cell(1,ntraj);
            
            for k=1:ntraj
                tid = obj.trajectoryIDs(k);
                ind = find(tidList==tid);
                if numel(ind) ~= 1
                    error('Trajectory either not found or found multiple times; argument varargin{1} must be incorrect.')
                end
                
                disp(['Fitting polys for TID = ' int2str(tid)])
                
                if ~isnan(roiAll{ind})
                    w.fit_polys_roi(roiAll{ind});
                end
                
                obj.polyFits = w.polyFits; % Put this here so that if WhiskerTrial had empty polyFits property
                                            % before calling w.mean_theta_and_kappa, it will now be transferred
                                            % to WhiskerSignalTrial.
                obj.polyFitsROI = w.polyFitsROI;
                obj.time{k} = w.get_time(tid);

            end
        end
        
        function obj = recompute_cached_mean_theta_kappa(obj, varargin)
            %
            %  obj = recompute_cached_mean_theta_kappa(obj, varargin)
            %
            % Recompute obj.theta, obj.kappa, and obj.time.
            %
            %   obj = recompute_cached_mean_theta_kappa(obj)
            %   obj = recompute_cached_mean_theta_kappa(obj, [xmin_position, xmax_position, xmin_curv, xmax_curv])
            %   obj = recompute_cached_mean_theta_kappa(obj, {[trajectoryIDs],[xmin_position1, xmax_position1, xmin_curv1, xmax_curv1], ...
            %                                     [xmin_position2, xmax_position2, xmin_curv2, xmax_curv2], etc})
            %
            %
            %   varargin{1}: pixel x-limits on which compute position and curvature.
            %       Sets an x-dimension ROI. This argument can be given in two forms:
            %           (1) an 1x4 vector that gives x-coordinates for *all* whiskers; or
            %           (2) a cell array where first element is a vector of trajectory IDs
            %               (of length N) and subsequent elements comprise N 1x4 vectors
            %               giving x-coordinates for the trajectory IDs specified in the first
            %               element (respectively).
            %
            %  UPDATE DESCRIPTION OF LIMITS WITH FOLLOWING:
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
            if nargin > 2
                error('Too many input arguments.')
            end
            if nargin==1
                tidList = obj.trajectoryIDs;
                XValLimsAll = cell(1,length(tidList)); % If no x-coordinates given for region of interest, all elts will be left empty.
            else
                xroi = varargin{1};
                if ~iscell(xroi) % Single set of x-coordinates given; not specified individually for different trajectories.
                    tidList = obj.trajectoryIDs;
                    XValLimsAll = cell(1,length(tidList));
                    % Copy single set of x-coordinates for each trajectory:
                    for k=1:length(XValLimsAll)
                        XValLimsAll{k} = xroi;
                    end
                else
                    tidList = xroi{1};
                    XValLimsAll = xroi(2:end);
                end
            end
            
            if length(tidList) ~= length(XValLimsAll)
                error(['First element of varargin{1} must be 1xN vector of trajectory IDs' ...
                    'and subsequent elements must comprise N 1x4 vectors giving x-coordinates']);
            end
            
            if length(tidList) ~= length(obj.trajectoryIDs)
                error('Number of trajectory IDs specified in varagin{1} does not match number of trajectory IDs in WhiskerTrial argument w.')
            end
            
            ntraj = length(obj.trajectoryIDs);
            
            obj.theta = cell(1,ntraj);
            obj.kappa = cell(1, ntraj);
            
            for k=1:ntraj
                tid = obj.trajectoryIDs(k);
                ind = find(tidList==tid);
                if numel(ind) ~= 1
                    error('Trajectory either not found or found multiple times; argument varargin{1} must be incorrect.')
                end
                XValLims = XValLimsAll{ind};
                if isempty(XValLims)
                    XValLimsPosition = [];
                    XValLimsCurv = [];
                else
                    XValLimsPosition = XValLims(1:2);
                    XValLimsCurv = XValLims(3:4);
                end
                
                disp(['Traj=' int2str(tid)])
                disp(['XValLimsPosition = ' num2str(XValLimsPosition)])
                disp(['XValLimsCurv = ' num2str(XValLimsCurv)])
                
                [t, theta, kappa] = obj.mean_theta_and_kappa(tid, XValLimsPosition, XValLimsCurv);
                
                obj.theta{k} = theta;
                obj.kappa{k} = kappa;
            end
            
        end
        
        function obj = recompute_cached_follicle_coords(obj, extrap_distance, varargin)
            %
            %  obj = recompute_cached_follicle_coords(obj, extrap_distance, varargin)
            %  obj = recompute_cached_follicle_coords(obj, extrap_distance)
            %  obj = recompute_cached_follicle_coords(obj, extrap_distance, tidList)
            %  obj = recompute_cached_follicle_coords(obj, extrap_distance, tidList, roi)
            %
            % Recomputes obj.follicleCoords.
            %
            %  INPUTS:
            %
            %   extrap_distance: The distance in units of pixels
            %        to extrapolate the whisker in order to get the
            %        (x,y) coordinates of the follicle. Extrapolation is
            %        based on the angle theta at the base of the whisker
            %        (radial distance = 0) if no mask is defined or if the
            %        whisker does not cross the mask.  Otherwise
            %        extrapolation is based on theta at the intersection of
            %        the whisker and the mask.
            %
            %   varargin{1}: Optional vector of trajectoryIDs. If varargin{1}
            %       is not given or is empty ([]), all trajectoryIDs in obj.trajectoryIDs
            %       are used. 
            %
            %   varargin{2}: Optional arc-length region of interest in which to
            %        fit a straight line to the whisker. This straight line is
            %        then extrapolated by extrap_distance to estimate the follicle
            %        position. Units of pixels.
            %        ** IMPORTANT**: Giving this varargin{2} argument
            %        fundamentally changes how the follicle position is estimated.
            %        If varargin{2} is *not* given, the most proximal theta value is used
            %        to set the direction of extrapolation.  See obj.get_follicle_coords().
            %
            %
            if nargin > 4
                error('Too many input arguments.')
            end
            if nargin > 3
                roi = varargin{2};
                fitLine = 1;
            else
                fitLine = 0;
            end
            if nargin==2
                tidList = obj.trajectoryIDs;
            else
                tidList = varargin{1};
                if isempty(tidList)
                    tidList = obj.trajectoryIDs;
                end
            end
            
            if extrap_distance <= 0 || length(extrap_distance) > 1
                error('Argument extrap_distance must be a positive scaler.')
            end
            
            obj.follicleExtrapDistInPix = extrap_distance; % Cache this argument in an object property.
                                                           % We need it later when computing moment of inertia.
            
            ntraj = length(obj.trajectoryIDs);
            obj.follicleCoordsX = cell(1,ntraj);
            obj.follicleCoordsY = cell(1,ntraj);
            if isempty(obj.time) % May have already been populated by recompute_cached_mean_theta_kappa()
                obj.time = cell(1,ntraj);
            end
            
            for k=1:length(tidList)
                tid = tidList(k);
                if ~ismember(tid,obj.trajectoryIDs)
                    error(['Trajectory ID ' int2str(tid) ' not found.'])
                end
                
                disp(['Computing follicle coordinates for TID = ' int2str(tid)])
                
                if fitLine==1
                    [y,x,t] = obj.get_follicle_coords(tid,extrap_distance,roi);
                else
                    [y,x,t] = obj.get_follicle_coords(tid,extrap_distance);
                end
                
                obj.follicleCoordsX{k} = x;
                obj.follicleCoordsY{k} = y;
                if isempty(obj.time{k}) % May have already been populated by recompute_cached_mean_theta_kappa()
                    obj.time{k} = t;
                end
                
            end
            
        end
        
        function tid = name2tid(obj, whisker_name)
            if ~ischar(whisker_name)
                error('Argument whisker_name must be a string.')
            end
            if numel(obj.whiskerNames) ~= numel(obj.trajectoryIDs)
                error('This WhiskerSignalTrial does not have matching whiskerNames and trajectoryIDs.')
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
        
        function  [pkVal,pkTime,pkBinarySignal] = get_whisk_peaks(obj,tid)
            %
            %   [pkVal,pkTime,pkBinarySignal] = get_whisk_peaks(obj,tid)
            %
            %
            %
            [y,t] = obj.get_position(tid);
            thresh = .01; % Should make argument. ***
            
            bandPassCutOffsInHz = [10 25];
            sampleRate = 500; % Should read from data.***
            W1 = bandPassCutOffsInHz(1) / (sampleRate/2);
            W2 = bandPassCutOffsInHz(2) / (sampleRate/2);
            [b,a]=butter(2,[W1 W2]);
            
            y = filtfilt(b,a,y);
            
            [pkVal,pkInd] = Whisker.lmax(y,0);
            
            ind = abs(pkVal) > thresh;
            pkVal = pkVal(ind);
            pkInd = pkInd(ind);
            
            pkTime = t(pkInd);
            
            pkBinarySignal = zeros(size(t));
            pkBinarySignal(pkInd) = 1;
        end
        
        function [f,Y] = get_FFT_position(obj,tid)
            %
            %   [f,Y] = get_FFT_position(obj,tid)
            %
            %
            Fs = 1 / obj.framePeriodInSec;
            [y,t] = obj.get_position(tid);
            
            % Check if frames are evenly spaced. If not, must interpolate theta for
            % missing frames prior to spectral analysis:
            frames = t ./ obj.framePeriodInSec;
            if length(unique(diff(frames))) > 1
                newframes = min(frames):max(frames);
                y = interp1(frames,y,newframes,'linear');
            end
            
            L = length(y);
            % Mean-subtract signal to eliminate DC term of transform:
            y = y - mean(y);
            NFFT = 2^nextpow2(L);
            %             h=hamming(L,'periodic')';
            %             y = y.*h;
            Y = fft(y,NFFT)/L;
            Y = 2*abs(Y(1:NFFT/2));
            f = Fs/2*linspace(0,1,NFFT/2);
        end
        
        function r = get_spectral_power_band_position(obj,tid,band)
            %
            %   r = get_spectral_power_band_position(obj,tid,band)
            %
            % band: [minFreqInclusive maxFreqInclusive]
            %
            [f,Y] = get_FFT_position(obj,tid);
            ind = f>=band(1) & f<=band(2);
            if isempty(ind)
                error('Invalid frequency band')
            end
            r = sum(Y(ind));
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
        
        function [th,t] = get_fitted_line_angle_in_roi(obj,tid,roi)
            %
            %    [th,t] = get_fitted_line_angle_in_roi(obj,tid,roi)
            %
            %  INPUTS:
            %
            %   tid: Trajectory ID.
            %
            %   roi: Arc-length region of interest in which to
            %        fit a straight line to the whisker. In units of pixels.
            %        The location of arc-length 0 depends on the mask settings.
            %        See documentation for object property maskTreatment.
            %
            %   RETURNS:
            %
            %   th: Whisker angle, theta, in degrees. Angle of 0 means that
            %       the whisker is parallel to the medial-lateral axis.
            %       Protraction corresponds to increasing angle.
            %       The value of obj.faceSideInImage must be set correctly to
            %       get correct angles from this method.
            %
            if numel(roi) ~= 2
                error('Argument roi must be a 2-element vector giving an arc-length ROI.')
            elseif roi(2) < roi(1)
                error('Argument roi  must give an arc-length ROI in format: [startInPix stopInPix].')
            end
            
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Could not find specified trajectory ID.')
            end
            
            if isempty(obj.polyFits)
                error('obj.polyFits is empty.')
            end
            
            if ~ischar(obj.faceSideInImage)
                error('(x,y) coordinate specification of face location not yet implemented.')
            end
            
            t = obj.get_time(tid);
            
            if all(isnan(t))
                disp(['All time values for trajectoryID=' int2str(tid) ' are NaN; perhaps '...
                    'nothing is tracked for this trial. Returning NaN.'])
                th = NaN;
                t = NaN;
                return
            end
            
            [R,THETA,KAPPA] = obj.arc_length_theta_and_kappa(tid);
            rind = cellfun(@(x) find(x >= roi(1) & x <= roi(2)), R,'UniformOutput',false);
            

            nframes = length(t);
            
            th = zeros(1,nframes);
            
            q = linspace(0,1);
            
            fittedX = obj.polyFits{ind}{1};
            fittedY = obj.polyFits{ind}{2};
            
            for k=1:nframes
                %                 disp(['Frame=' int2str(k)])
                px = fittedX(k,:);
                py = fittedY(k,:);
                
                xall = polyval(px,q);
                yall = polyval(py,q);
                
                % If maskTreatment property is set to 'maskNaN', and the whisker
                % does not cross the mask, then all the radial distance values in R
                % for that whisker for that frame will be set to NaN.  In that case,
                % there will be an empty element in rind.  Propagate the NaN here:
                if isempty(rind{k})
                    th(k) = NaN;
                    continue
                end
                
                x = xall(rind{k});
                y = yall(rind{k});
                
                if numel(x) < 2 || numel(y) < 2
                    th(k) = NaN;
                    continue
                end
                
                if strcmp(obj.faceSideInImage,'right') || strcmp(obj.faceSideInImage,'left')
                    if numel(unique(x)) < 2
                        disp('Found whisker at either +180 deg or -180 deg; ambiguous--setting to NaN')
                    end
                end
                
                if strcmp(obj.faceSideInImage,'right') && strcmp(obj.protractionDirection,'upward') %
                    p = polyfit(x,y,1);
                    th(k) = atand(p(1));
                elseif strcmp(obj.faceSideInImage,'right') && strcmp(obj.protractionDirection,'downward') %
                    p = polyfit(x,y,1);
                    th(k) = -atand(p(1));
                elseif strcmp(obj.faceSideInImage,'left') && strcmp(obj.protractionDirection,'upward') %
                    p = polyfit(x,y,1);
                    th(k) = -atand(p(1));
                elseif strcmp(obj.faceSideInImage,'left') && strcmp(obj.protractionDirection,'downward') %
                    p = polyfit(x,y,1);
                    th(k) = atand(p(1));
                elseif strcmp(obj.faceSideInImage,'top') && strcmp(obj.protractionDirection,'rightward') % 
                    p = polyfit(y,x,1);
                    th(k) = atand(p(1));
                elseif strcmp(obj.faceSideInImage,'top') && strcmp(obj.protractionDirection,'leftward') %
                    p = polyfit(y,x,1);
                    th(k) = -atand(p(1));
                elseif strcmp(obj.faceSideInImage,'bottom') && strcmp(obj.protractionDirection,'rightward') %
                    p = polyfit(y,x,1);
                    th(k) = -atand(p(1));
                elseif strcmp(obj.faceSideInImage,'bottom') && strcmp(obj.protractionDirection,'leftward') %
                    p = polyfit(y,x,1);
                    th(k) = atand(p(1));
                else
                    error('Invalid value of obj.faceSideInImage.')
                end
            end
            
            ind = isinf(th);
            if sum(ind) > 0
                disp('Found frames with whisker at either +180 deg or -180 deg; ambiguous--setting to NaN')
                th(ind) = NaN;
            end
            
        end
        
        function [y,x,t] = get_follicle_coords(obj,tid,extrap_distance,varargin)
            %
            %    [y,x,t] = get_follicle_coords(obj,tid,extrap_distance,varargin)
            %
            %  INPUTS:
            %
            %   tid: Trajectory ID.
            %
            %   extrap_distance: The distance in units of pixels
            %        to extrapolate the whisker in order to get the
            %        (x,y) coordinates of the follicle. Extrapolation is
            %        based on the angle theta at the base of the whisker
            %        (radial distance = 0) if no mask is defined or if the
            %        whisker does not cross the mask.  Otherwise
            %        extrapolation is based on theta at the intersection of
            %        the whisker and the mask.
            %
            %   varargin{1}: Optional arc-length region of interest in which to
            %        fit a straight line to the whisker. This straight line is
            %        then extrapolated by extrap_distance to estimate the follicle
            %        position. Units of pixels.
            %        ** IMPORTANT**: Giving this varargin{1} argument
            %        fundamentally changes how the follicle position is estimated.
            %        If varargin{1} is *not* given, the most proximal theta value is used
            %        (as described above) to set the direction of extrapolation.
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
            npoints = 100;  % Number of points to use in reconstructing whisker from fitted polynomials.
                            % Sets the resolution at which the follicle coordinates can be estimated when.
                            % A higher value is better but runs more slowly. Doesn't have any effect unless
                            % a mask is being used. 
            if nargin < 4
                fitLine = 0;
            else
                fitLine = 1;
                roiLine = varargin{1};
                if numel(roiLine) ~= 2
                    error('varargin{1} must be a 2-element vector giving an arc-length ROI.')
                elseif roiLine(2) < roiLine(1)
                    error('varargin{1} must give an arc-length ROI in format: [startInPix stopInPix].')
                end
            end
            
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Could not find specified trajectory ID.')
            end
            
            if isempty(obj.polyFits)
                error('obj.polyFits is empty.')
            end
            
%             if all(isnan(obj.theta{ind}))
%                 disp(['All theta values for trajectoryID=' int2str(tid) ' are NaN; perhaps '...
%                     'nothing is tracked for this trial. Setting face follicle coordinates to NaN.'])
%                 x = NaN;
%                 y = NaN;
%                 t = NaN;
%                 return
%             end
            
            [R,THETA,KAPPA] = obj.arc_length_theta_and_kappa(tid,npoints);
            rind = cellfun(@(x) find(x==0), R,'UniformOutput',false);
            emp = find(cellfun(@isempty, rind));
            
            if fitLine==1 % Get angle based on line fitted to ROI.
                [th,t] = obj.get_fitted_line_angle_in_roi(tid,roiLine);
                if ~isempty(emp)
                    for q=1:length(emp)
                        th(emp(q)) = NaN;
                    end
                end
            else % Use theta at point where arc-length = 0.
                th = cellfun(@(x,y) x(y), THETA, rind,'UniformOutput',false);
                if ~isempty(emp)
                    % If maskTreatment property is set to 'maskNaN', and the whisker
                    % does not cross the mask, then all the radial distance values in R
                    % for that whisker for that frame will be set to NaN.  In that case,
                    % there will be an empty element in rind. Propagate the NaN here:
                    for q=1:length(emp)
                        th{emp(q)} = NaN;
                    end
                end
                th = cell2mat(th);
                t = obj.get_time(tid);
            end
            
            nframes = length(t);
            
            if nframes==0
                disp(['Number of frames for trajectoryID=' int2str(tid) ' is 0; perhaps '...
                    'nothing is tracked for this trial. Setting face follicle coordinates to NaN.'])
                x = NaN;
                y = NaN;
                t = NaN;
                return
            end

            x = zeros(1,nframes);
            y = zeros(1,nframes);
            
            q = linspace(0,1,npoints);
            
            fittedX = obj.polyFits{ind}{1};
            fittedY = obj.polyFits{ind}{2};
            
            for k=1:nframes
                %                 disp(['Frame=' int2str(k)])
                px = fittedX(k,:);
                py = fittedY(k,:);
                
                xall = polyval(px,q);
                yall = polyval(py,q);
                
                
                % If maskTreatment property is set to 'maskNaN', and the whisker
                % does not cross the mask, then all the radial distance values in R
                % for that whisker for that frame will be set to NaN.  In that case,
                % there will be an empty element in rind.  Propagate the NaN here:
                if isempty(rind{k})
                    x(k) = NaN;
                    y(k) = NaN;
                    % Alternatively, could just use (x,y,theta) values for point closest to
                    % face by uncommenting next two lines:
                    % xall(end);
                    % xall(end);
                else
                    x(k) = xall(rind{k});
                    y(k) = yall(rind{k});
                end
            end
            
            
            % Protraction means theta is increasing. 
            % Theta is 0 when perpendicular to the midline of the mouse.
            
            if strcmp(obj.faceSideInImage,'top') && strcmp(obj.protractionDirection,'rightward')
                deltax = -extrap_distance*sind(th);
                deltay = -extrap_distance*cosd(th);
                
            elseif strcmp(obj.faceSideInImage,'top') && strcmp(obj.protractionDirection,'leftward')
                deltax = extrap_distance*sind(th);
                deltay = -extrap_distance*cosd(th);
                
            elseif strcmp(obj.faceSideInImage,'left') && strcmp(obj.protractionDirection,'downward')
                deltax = -extrap_distance*cosd(th);
                deltay = -extrap_distance*sind(th);
            
            elseif strcmp(obj.faceSideInImage,'left') && strcmp(obj.protractionDirection,'upward')
                deltax = -extrap_distance*cosd(th);
                deltay = extrap_distance*sind(th);
            
            elseif strcmp(obj.faceSideInImage,'right') && strcmp(obj.protractionDirection,'upward')
                deltax = extrap_distance*cosd(th);
                deltay = extrap_distance*sind(th);
                
            elseif strcmp(obj.faceSideInImage,'right') && strcmp(obj.protractionDirection,'downward')
                deltax = extrap_distance*cosd(th);
                deltay = -extrap_distance*sind(th);            
                
            elseif strcmp(obj.faceSideInImage,'bottom') && strcmp(obj.protractionDirection,'rightward')
                deltax = -extrap_distance*sind(th);
                deltay = extrap_distance*cosd(th);  
            
            elseif strcmp(obj.faceSideInImage,'bottom') && strcmp(obj.protractionDirection,'leftward')
                deltax = extrap_distance*sind(th);
                deltay = extrap_distance*cosd(th);    
            else
                error('Invalid value of property ''faceSideInImage'' or ''protractionDirection''')
            end
            
            x = x + deltax;
            y = y + deltay;
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
                    'Must run obj.recompute_cached_follicle_coords before this method.'])
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
        
        function obj = set_mask_from_points(obj,tid,x,y)
            %
            % Sets obj.polyFitsMask in order
            % to create a mask defined by the points in x and y.
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
            
            if isempty(obj.polyFitsMask)
                obj.polyFitsMask = cell(1,length(obj.trajectoryIDs));
            end
            
            for k=1:length(tid)
                ind = obj.trajectoryIDs==tid(k);
                if max(ind) < 1
                    error('Trajectory ID was not found.')
                end
                obj.polyFitsMask{ind} = {px,py};
            end
        end
        
        function obj = set_bar_offset(obj,dx,dy)
            %
            % obj = set_bar_offset(obj,dx,dy)
            %
            % dx: Number of pixels to offset bar center in x.
            % dy: Number of pixels to offset bar center in y.
            %
            %
            %
            if length(dx) ~= 1 || length(dy) ~= 1
                error('Arguments dx and dy must both be scaler.')
            end
            
            obj.barPosOffset = [dx dy];
            
            % Should also add capability to set different offsets for
            % different frames.
        end
        
        function plot_mask(obj,tid,varargin)
            %
            % Plots the polynomial mask defined by obj.polyFitsMask.
            %
            % tid: A single trajectory ID.
            %
            % varargin{1}: Plot symbol string, e.g., 'k-'.  A string
            %               that can be given as an argument to plot().
            %
            % varargin{2}: Tracked frame number, from 1 to the number of
            %        frames tracked. Not necessarily frame number from
            %        the original movie (unless all frames were tracked).
            %        If the mask is the same for all frames, this argument
            %        is ignored.
            %
            %
            %
            if isempty(obj.polyFitsMask)
                disp('obj.polyFitsMask is empty; nothing to plot.')
                return
            end
            
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Trajectory ID was not found.')
            end
            
            if isempty(obj.polyFitsMask{ind})
                disp(['obj.polyFitsMask for tid ' int2str(tid) 'is empty; nothing to plot.'])
                return
            end
            
            if length(tid) > 1
                error('Only a single trajectory ID allowed.')
            end
            
            if nargin > 4
                error('Too many input arguments.')
            end
            
            if nargin > 2
                plotString = varargin{1};
            else
                plotString = 'k-';
            end
            
            if nargin > 3
                frame = varargin{2};
                if frame > length(obj.time{ind})
                    error('varargin{2}, the frame number, exceeds the number of tracked frames.')
                end
            else
                frame = 1;
            end
            
            px = obj.polyFitsMask{ind}{1};
            py = obj.polyFitsMask{ind}{2};
            
            if size(px,1) > 1
                px = px(frame,:);
                py = py(frame,:);
            end
            
            q = linspace(0,1);
            
            x = polyval(px,q);
            y = polyval(py,q);
            
            plot(x,y,plotString,'LineWidth',2)
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
        
%         function [y,t] = get_follicle_arc_length_position(obj,tid)
%             %
%             % [y,t] = get_follicle_arc_length_position(obj,tid)
%             %
%             % Fits polynomials to the the follicle x and y coordinates, parameterizing
%             % the path they take along the face as c(q) = (x(q),y(q)), where q is in [0,1]. 
%             % Arc length at each point along this path is then computed. For each frame
%             % the closest point on the curve to the follicle is determined. The arc-length
%             % location of this point is returned in r.
%             %
%             % obj.follicleCoordsX and obj.follicleCoordsY must be non-empty (i.e.
%             % obj.recompute_cached_follicle_coords() must have been called.
%             %
%             %
%             %
%             
%             
%             
%         end

        function [R,THETA,KAPPA,varargout] = arc_length_theta_and_kappa(obj,tid,varargin)
            %
            % [R,THETA,KAPPA] = arc_length_theta_and_kappa(obj,tid)
            % [R,THETA,KAPPA,Y,X] = arc_length_theta_and_kappa(obj,tid)
            % [R,THETA,KAPPA] = arc_length_theta_and_kappa(obj,tid,npoints)
            %
            %   The whisker is parameterized as c(q) = (x(q),y(q)), where q 
            %   is in [0,1].
            %
            %   For this function to work properly, object properties 'faceSideInImage'
            %   and 'imagePixelDimsXY' must be correctly set.
            %
            % INPUTS:
            %
            %   tid: Whisker trajectory ID.
            %
            %   varargin{1}: Optional, integer giving number of points (values of q) to 
            %        use in reconstructing each whisker. Default is 100 points.  
            %    
            %
            % RETURNS:
            %
            %   R:  A cell array where each element is the arc length, computed moving outward from
            %       whisker follicle along the whisker for a single frame. Units of pixels.
            %
            %   THETA: A cell array where each element theta is the angle of the line tangent to the
            %           whisker (i.e., to c(q)) at each value of q. In
            %           degrees.
            %
            %   KAPPA: A cell array where each element kappa is the signed curvature at each point
            %           on the whisker (i.e., for each value of q). Units of 1/pixels. Abs(kappa(q)) is 1/X where X is
            %           the radius in pixels of the osculating circle at c(q).
            %
            %   Optionally also:
            %
            %   X,Y: Cell arrays containing the x and y image (pixel) coordinates corresponding to the
            %        values in R, THETA, and KAPPA.
            %
            %
            %
            %   kappa(q) = (x'y'' - y'x'') / (x'^2 + y'^2)^(3/2)
            %   theta(q) = atand(y'/x')
            %   arc_length(q) = cumsum(sqrt(x'^2 + y'^2))
            %
            %
            if nargin < 3
                npoints = 100;
            else
                npoints = varargin{1};
                if isempty(npoints)
                    npoints = 100;
                end
            end
            
            ind = find(obj.trajectoryIDs == tid);
            if isempty(ind)
                error('Could not find specified trajectory ID.')
            end
            
            nframes = size(obj.polyFits{ind}{1},1);
            
            if isempty(obj.polyFits)
                error('obj.polyFits is empty.')
            end
            
            R = cell(1,nframes);
            THETA = cell(1,nframes);
            KAPPA = cell(1,nframes);
            
            if nargout>3
                X = cell(1,nframes);
                Y = cell(1,nframes);
            end
            
            fittedX = obj.polyFits{ind}{1};
            fittedY = obj.polyFits{ind}{2};
            
            q = linspace(0,1,npoints);
            
            for k=1:nframes
                
                px = fittedX(k,:);
                py = fittedY(k,:);
                
                pxDot = polyder(px);
                pxDoubleDot = polyder(pxDot);
                
                pyDot = polyder(py);
                pyDoubleDot = polyder(pyDot);
                
                xDot = polyval(pxDot,q);
                xDoubleDot = polyval(pxDoubleDot,q);
                
                yDot = polyval(pyDot,q);
                yDoubleDot = polyval(pyDoubleDot,q);
                
                dq = [0 diff(q)];
                
                % Arc length as a function of q, after integration below:
                R{k} = cumsum(sqrt(xDot.^2 + yDot.^2) .* dq); % arc length segments, in pixels, times dq.
                
                
                % Angle (in degrees) as a function of q:
                % Protraction means theta is increasing.
                % Theta is 0 when perpendicular to the midline of the mouse.
                if strcmp(obj.faceSideInImage,'top') && strcmp(obj.protractionDirection,'rightward')
                    THETA{k} = atand(xDot ./ yDot);
                elseif strcmp(obj.faceSideInImage,'top') && strcmp(obj.protractionDirection,'leftward')
                    THETA{k} = -atand(xDot ./ yDot);
                elseif strcmp(obj.faceSideInImage,'left') && strcmp(obj.protractionDirection,'downward')
                    THETA{k} = atand(yDot ./ xDot); 
                elseif strcmp(obj.faceSideInImage,'left') && strcmp(obj.protractionDirection,'upward')
                    THETA{k} = -atand(yDot ./ xDot); 
                elseif strcmp(obj.faceSideInImage,'right') && strcmp(obj.protractionDirection,'upward')
                    THETA{k} = atand(yDot ./ xDot); 
                elseif strcmp(obj.faceSideInImage,'right') && strcmp(obj.protractionDirection,'downward')
                    THETA{k} = -atand(yDot ./ xDot);
                elseif strcmp(obj.faceSideInImage,'bottom') && strcmp(obj.protractionDirection,'rightward')
                    THETA{k} = -atand(xDot ./ yDot);
                elseif strcmp(obj.faceSideInImage,'bottom') && strcmp(obj.protractionDirection,'leftward')
                    THETA{k} = atand(xDot ./ yDot);
                else
                    error('Invalid value of property ''faceSideInImage'' or ''protractionDirection''')
                end
                
                
                % Signed curvature as a function of q:
                KAPPA{k} = (xDot.*yDoubleDot - yDot.*xDoubleDot) ./ ((xDot.^2 + yDot.^2).^(3/2)); % SIGNED CURVATURE, in 1/pixels.
                %                 KAPPA{k} = abs(xDot.*yDoubleDot - yDot.*xDoubleDot) ./ ((xDot.^2 + yDot.^2).^(3/2)); % CURVATURE, in 1/pixels.
                
                if nargout>3
                    X{k} = polyval(px,q);
                    Y{k} = polyval(py,q);
                    varargout{1} = Y;
                    varargout{2} = X;
                end
            end
            
            % Apply any mask:
            if iscell(obj.maskTreatment)
                if length(obj.maskTreatment) ~= length(obj.trajectoryIDs)
                    error('obj.maskTreatment and obj.trajectoryIDs must have the same length and matching entries.')
                end
                mask_treatment = obj.maskTreatment{ind};
            else
                mask_treatment = obj.maskTreatment;
            end
            
            if strcmp(mask_treatment,'none')
                return
            end
            
            % If there is a polynomial mask specified (see documentation for
            % object property 'polyFitsMask') and varargin{1} is given, subtract from each element of
            % R the radial distance at the intersection, if any, of the fitted
            % whisker and the polynomial mask.
            if isempty(obj.polyFitsMask)
                return
            else
                pm = obj.polyFitsMask{ind};
                if isempty(pm)
                    return
                end
            end
            
            fittedXMask = obj.polyFitsMask{ind}{1};
            fittedYMask = obj.polyFitsMask{ind}{2};
            
            q = linspace(0,1,npoints);
            
            for k=1:nframes
                
                px = fittedX(k,:);
                py = fittedY(k,:);
                
                if size(fittedXMask,1) > 1
                    pxm = fittedXMask(k,:);
                    pym = fittedYMask(k,:);
                else
                    pxm = fittedXMask;
                    pym = fittedYMask;
                end
                
                C1 = [polyval(px,q); polyval(py,q)];
                C2 = [polyval(pxm,q); polyval(pym,q)];
                
                P = Whisker.InterX(C1,C2); % Find points where whisker and mask curves intersect. Slower but more
                                           % accurate version that isn't limited in resolution by the number of
                                           % points whisker and mask are evaluated at.
                if size(P,2) > 1   % Don't need for faster version, which handles this.
                    disp('Found more than 1 intersection of whisker and mask curves; using only first.')
                    P = P(:,1);
                end
                                           
%                 P = Whisker.InterXFast(C1,C2); % Find points where whisker and mask curves intersect. Much faster version
%                                                 % that is limited in resolution by the number of
%                                                 % points whisker and mask are evaluated at (i.e., by number of points in q).
                                                
                if isempty(P)
                    if strcmp(mask_treatment,'maskNaN')
                        R{k} = nan(size(R{k}));
                    end
                else
                    % Find at what q the whisker is at (P(1),P(2)), i.e., q s.t. x(q)=P(1),y(q)=P(2).
                    % Doesn't match exactly (maybe due to roundoff error), so find closest.
                    C = C1 - repmat(P,[1 size(C1,2)]);
                    err = sqrt(C(1,:).^2 + C(2,:).^2);
                    ind2 = err==min(err);
                    R{k} = R{k} - R{k}(ind2);
                end
            end
        end
        
        function [R,THETA,varargout] = arc_length_and_theta(obj,tid,varargin)
            % 
            % POINT OF THIS IS THAT IT'S FASTER THAN arc_length_theta_and_kappa when
            % KAPPA IS NOT NEEDED.****
            %
            % [R,THETA] = arc_length_theta_and_kappa(obj,tid)
            % [R,THETA,Y,X] = arc_length_theta_and_kappa(obj,tid)
            % [R,THETA] = arc_length_theta_and_kappa(obj,tid,npoints)
            %
            %   The whisker is parameterized as c(q) = (x(q),y(q)), where q 
            %   is in [0,1].
            %
            %   For this function to work properly, object properties 'faceSideInImage'
            %   and 'imagePixelDimsXY' must be correctly set.
            %
            % INPUTS:
            %
            %   tid: Whisker trajectory ID.
            %
            %   varargin{1}: Optional, integer giving number of points (values of q) to 
            %        use in reconstructing each whisker. Default is 100 points.  
            %    
            %
            % RETURNS:
            %
            %   R:  A cell array where each element is the arc length, computed moving outward from
            %       whisker follicle along the whisker for a single frame. Units of pixels.
            %
            %   THETA: A cell array where each element theta is the angle of the line tangent to the
            %           whisker (i.e., to c(q)) at each value of q. In
            %           degrees.
            %
            %
            %   Optionally also:
            %
            %   X,Y: Cell arrays containing the x and y image (pixel) coordinates corresponding to the
            %        values in R, THETA, and KAPPA.
            %
            %
            %
            %   kappa(q) = (x'y'' - y'x'') / (x'^2 + y'^2)^(3/2)
            %   theta(q) = atand(y'/x')
            %   arc_length(q) = cumsum(sqrt(x'^2 + y'^2))
            %
            %
            if nargin < 3
                npoints = 100;
            else
                npoints = varargin{1};
                if isempty(npoints)
                    npoints = 100;
                end
            end
            
            ind = find(obj.trajectoryIDs == tid);
            if isempty(ind)
                error('Could not find specified trajectory ID.')
            end
            
            nframes = size(obj.polyFits{ind}{1},1);
            
            if isempty(obj.polyFits)
                error('obj.polyFits is empty.')
            end
            
            R = cell(1,nframes);
            THETA = cell(1,nframes);
            
            if nargout>3
                X = cell(1,nframes);
                Y = cell(1,nframes);
            end
            
            fittedX = obj.polyFits{ind}{1};
            fittedY = obj.polyFits{ind}{2};
            
            q = linspace(0,1,npoints);
            
            for k=1:nframes
                
                px = fittedX(k,:);
                py = fittedY(k,:);
                
                pxDot = polyder(px);
                
                pyDot = polyder(py);
                
                xDot = polyval(pxDot,q);
                
                yDot = polyval(pyDot,q);
                
                dq = [0 diff(q)];
                
                % Arc length as a function of q, after integration below:
                R{k} = cumsum(sqrt(xDot.^2 + yDot.^2) .* dq); % arc length segments, in pixels, times dq.
                
                
                % Angle (in degrees) as a function of q:
                % Protraction means theta is increasing.
                % Theta is 0 when perpendicular to the midline of the mouse.
                if strcmp(obj.faceSideInImage,'top') && strcmp(obj.protractionDirection,'rightward')
                    THETA{k} = atand(xDot ./ yDot);
                elseif strcmp(obj.faceSideInImage,'top') && strcmp(obj.protractionDirection,'leftward')
                    THETA{k} = -atand(xDot ./ yDot);
                elseif strcmp(obj.faceSideInImage,'left') && strcmp(obj.protractionDirection,'downward')
                    THETA{k} = atand(yDot ./ xDot); 
                elseif strcmp(obj.faceSideInImage,'left') && strcmp(obj.protractionDirection,'upward')
                    THETA{k} = -atand(yDot ./ xDot); 
                elseif strcmp(obj.faceSideInImage,'right') && strcmp(obj.protractionDirection,'upward')
                    THETA{k} = atand(yDot ./ xDot); 
                elseif strcmp(obj.faceSideInImage,'right') && strcmp(obj.protractionDirection,'downward')
                    THETA{k} = -atand(yDot ./ xDot);
                elseif strcmp(obj.faceSideInImage,'bottom') && strcmp(obj.protractionDirection,'rightward')
                    THETA{k} = -atand(xDot ./ yDot);
                elseif strcmp(obj.faceSideInImage,'bottom') && strcmp(obj.protractionDirection,'leftward')
                    THETA{k} = atand(xDot ./ yDot);
                else
                    error('Invalid value of property ''faceSideInImage'' or ''protractionDirection''')
                end
                
                if nargout>3
                    X{k} = polyval(px,q);
                    Y{k} = polyval(py,q);
                    varargout{1} = Y;
                    varargout{2} = X;
                end
            end
            
            % Apply any mask:
            if iscell(obj.maskTreatment)
                if length(obj.maskTreatment) ~= length(obj.trajectoryIDs)
                    error('obj.maskTreatment and obj.trajectoryIDs must have the same length and matching entries.')
                end
                mask_treatment = obj.maskTreatment{ind};
            else
                mask_treatment = obj.maskTreatment;
            end
            
            if strcmp(mask_treatment,'none')
                return
            end
            
            % If there is a polynomial mask specified (see documentation for
            % object property 'polyFitsMask') and varargin{1} is given, subtract from each element of
            % R the radial distance at the intersection, if any, of the fitted
            % whisker and the polynomial mask.
            if isempty(obj.polyFitsMask)
                return
            else
                pm = obj.polyFitsMask{ind};
                if isempty(pm)
                    return
                end
            end
            
            fittedXMask = obj.polyFitsMask{ind}{1};
            fittedYMask = obj.polyFitsMask{ind}{2};
            
            q = linspace(0,1,npoints);
            
            for k=1:nframes
                
                px = fittedX(k,:);
                py = fittedY(k,:);
                
                if size(fittedXMask,1) > 1
                    pxm = fittedXMask(k,:);
                    pym = fittedYMask(k,:);
                else
                    pxm = fittedXMask;
                    pym = fittedYMask;
                end
                
                C1 = [polyval(px,q); polyval(py,q)];
                C2 = [polyval(pxm,q); polyval(pym,q)];
                
                P = Whisker.InterX(C1,C2); % Find points where whisker and mask curves intersect. Slower but more
                                           % accurate version that isn't limited in resolution by the number of
                                           % points whisker and mask are evaluated at.
                if size(P,2) > 1   % Don't need for faster version, which handles this.
                    disp('Found more than 1 intersection of whisker and mask curves; using only first.')
                    P = P(:,1);
                end
                                           
%                 P = Whisker.InterXFast(C1,C2); % Find points where whisker and mask curves intersect. Much faster version
%                                                 % that is limited in resolution by the number of
%                                                 % points whisker and mask are evaluated at (i.e., by number of points in q).
                                                
                if isempty(P)
                    if strcmp(mask_treatment,'maskNaN')
                        R{k} = nan(size(R{k}));
                    end
                else
                    % Find at what q the whisker is at (P(1),P(2)), i.e., q s.t. x(q)=P(1),y(q)=P(2).
                    % Doesn't match exactly (maybe due to roundoff error), so find closest.
                    C = C1 - repmat(P,[1 size(C1,2)]);
                    err = sqrt(C(1,:).^2 + C(2,:).^2);
                    ind2 = err==min(err);
                    R{k} = R{k} - R{k}(ind2);
                end
            end
        end
        
        function [R,THETA,KAPPA,varargout] = arc_length_theta_and_kappa_in_roi(obj,tid,varargin)
            %
            % [R,THETA,KAPPA] = arc_length_theta_and_kappa_in_roi(obj,tid)
            % [R,THETA,KAPPA,Y,X] = arc_length_theta_and_kappa_in_roi(obj,tid)
            % [R,THETA,KAPPA] = arc_length_theta_and_kappa_in_roi(obj,tid,npoints)
            %
            %   The whisker is parameterized as c(q) = (x(q),y(q)), where q 
            %   is in [0,1]. However, here we deal with a second representation of the
            %   whisker, where polynomials were fitted over a constant region of arc-length,
            %   [q0,q1], q0,q1 both within [0,1]. 
            %
            %   For this function to work properly, object properties 'faceSideInImage'
            %   and 'imagePixelDimsXY' must be correctly set.
            %
            % INPUTS:
            %
            %   tid: Whisker trajectory ID.
            %
            %   varargin{1}: Optional, integer giving number of points (values of q) to 
            %        use in reconstructing each whisker. Default is 100 points.  
            %    
            %
            % RETURNS:
            %
            %   R:  A cell array where each element is the arc length, computed moving outward from
            %       whisker follicle along the whisker for a single frame, but only those values
            %       falling within the ROI. Units of pixels.
            %
            %   THETA: A cell array where each element theta is the angle of the line tangent to the
            %           whisker (i.e., to c(q)) at each value of q in the ROI. In
            %           degrees.
            %
            %   KAPPA: A cell array where each element kappa is the signed curvature at each point
            %           on the whisker within the ROI. Units of 1/pixels. Abs(kappa(q)) is 1/X where X is
            %           the radius in pixels of the osculating circle at c(q). Uses the secondary polynomials
            %           fitted within the ROI.
            %
            %   Optionally also:
            %
            %   X,Y: Cell arrays containing the x and y image (pixel) coordinates corresponding to the
            %        values in R, THETA, and KAPPA.
            %
            %
            %
            %
            if nargin < 3
                npoints = 100;
            else
                npoints = varargin{1};
                if isempty(npoints)
                    npoints = 100;
                end
            end
            
            ind = find(obj.trajectoryIDs == tid);
            if isempty(ind)
                error('Could not find specified trajectory ID.')
            end
            
            nframes = size(obj.polyFitsROI{ind}{1},1);
            
            if isempty(obj.polyFitsROI)
                error('obj.polyFitsROI is empty.')
            end
            
            R = cell(1,nframes);
            THETA = cell(1,nframes);
            KAPPA = cell(1,nframes);
            
            if nargout>3
                X = cell(1,nframes);
                Y = cell(1,nframes);
            end
            
            fittedX = obj.polyFitsROI{ind}{1};
            fittedY = obj.polyFitsROI{ind}{2};
            fittedQ = obj.polyFitsROI{ind}{3};
            
            for k=1:nframes
                
                px = fittedX(k,:);
                py = fittedY(k,:);
                pq = fittedQ(k,:);
                
                q = linspace(pq(1),pq(2),npoints);
                
                pxDot = polyder(px);
                pxDoubleDot = polyder(pxDot);
                
                pyDot = polyder(py);
                pyDoubleDot = polyder(pyDot);
                
                xDot = polyval(pxDot,q);
                xDoubleDot = polyval(pxDoubleDot,q);
                
                yDot = polyval(pyDot,q);
                yDoubleDot = polyval(pyDoubleDot,q);
                
                dq = [0 diff(q)];
                
                % Arc length as a function of q, after integration below:
                R{k} = cumsum(sqrt(xDot.^2 + yDot.^2) .* dq); % arc length segments, in pixels, times dq.
                
                
                R{k} = R{k} + pq(3); % Add the whole-whisker arc-length of the first point in the ROI to make this quantity the arc-length
                                     % measured over the whole fitted whisker, not just arc-length over the ROI.
                                     % Note that pq(3) already accounts for mask, so do not need to do so in this method.
                

                % Angle (in degrees) as a function of q:
                % Protraction means theta is increasing.
                % Theta is 0 when perpendicular to the midline of the mouse.
                if strcmp(obj.faceSideInImage,'top') && strcmp(obj.protractionDirection,'rightward')
                    THETA{k} = atand(xDot ./ yDot);
                elseif strcmp(obj.faceSideInImage,'top') && strcmp(obj.protractionDirection,'leftward')
                    THETA{k} = -atand(xDot ./ yDot);
                elseif strcmp(obj.faceSideInImage,'left') && strcmp(obj.protractionDirection,'downward')
                    THETA{k} = atand(yDot ./ xDot); 
                elseif strcmp(obj.faceSideInImage,'left') && strcmp(obj.protractionDirection,'upward')
                    THETA{k} = -atand(yDot ./ xDot); 
                elseif strcmp(obj.faceSideInImage,'right') && strcmp(obj.protractionDirection,'upward')
                    THETA{k} = atand(yDot ./ xDot); 
                elseif strcmp(obj.faceSideInImage,'right') && strcmp(obj.protractionDirection,'downward')
                    THETA{k} = -atand(yDot ./ xDot);
                elseif strcmp(obj.faceSideInImage,'bottom') && strcmp(obj.protractionDirection,'rightward')
                    THETA{k} = -atand(xDot ./ yDot);
                elseif strcmp(obj.faceSideInImage,'bottom') && strcmp(obj.protractionDirection,'leftward')
                    THETA{k} = atand(xDot ./ yDot);
                else
                    error('Invalid value of property ''faceSideInImage'' or ''protractionDirection''')
                end
                
                % Signed curvature as a function of q:
                KAPPA{k} = (xDot.*yDoubleDot - yDot.*xDoubleDot) ./ ((xDot.^2 + yDot.^2).^(3/2)); % SIGNED CURVATURE, in 1/pixels.
                %                 KAPPA{k} = abs(xDot.*yDoubleDot - yDot.*xDoubleDot) ./ ((xDot.^2 + yDot.^2).^(3/2)); % CURVATURE, in 1/pixels.
                
                if nargout>3
                    X{k} = polyval(px,q);
                    Y{k} = polyval(py,q);
                    varargout{1} = Y;
                    varargout{2} = X;
                end
            end
            
            % Do not need to apply any mask, because pq(3) above already accounts for mask.

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
        
        function [rNearest,thetaNearest,kappaNearest,yNearest,xNearest,dist,t] = get_r_theta_kappa_nearest_bar(obj,tid,proximity_threshold)
            %
            % USAGE:
            %   [rNearest,thetaNearest,kappaNearest,yNearest,xNearest,dist,t] = get_r_theta_kappa_nearest_bar(obj,tid)
            %   [rNearest,thetaNearest,kappaNearest,yNearest,xNearest,dist,t] = get_r_theta_kappa_nearest_bar(obj,tid,proximity_threshold)
            %
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %       proximity_threshold:  Optional argument giving distance from nearest point on whisker
            %               to bar center, in units of bar radius, beyond which
            %               the whisker will be extrapolated along the last theta in
            %               order to determine distance between whisker and bar.
            %
            %   OUTPUTS:
            %       rNearest: Arc-length (radial) distance along whisker to point nearest
            %                 center of the bar. Units of pixels.
            %       thetaNearest: Theta at rNearest.
            %       kappaNearest: Kappa at rNearest.
            %       YNearest: Image-coordinate Y value at rNearest. This estimates point of contact. In pixels.
            %       XNearest: Image-coordinate X value at rNearest. This estimates point of contact. In pixels.
            %       dist: Distance from bar center to nearest point on whisker. Units of pixels.
            %       t: The corresponding times of each observation.
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
            
            if nargin < 3
                proximity_threshold = -1;
            end
            
            t = obj.time{ind};
            f = t / obj.framePeriodInSec;
            nframes = length(f);
            
            % Add any offset to the tracked bar position:
            if isempty(obj.barPosOffset);
                bp = obj.barPosClean;
            else
                bp = obj.barPosClean;
                if size(obj.barPosOffset,1)==1
                    bp = bp + repmat([0 obj.barPosOffset],size(bp,1),1);
                elseif size(obj.barPosOffset,1)==size(bp,1) && size(obj.barPosOffset,2)==2
                    bp = bp + obj.barPosOffset;
                else
                    error('Size of obj.barPosOffset is not valid; must be either 1 x 2 or nframes x 2.')
                end
            end
            
            bar_f = bp(:,1);
            bar_x = bp(:,2);
            bar_y = bp(:,3);
            
            [R,THETA,KAPPA] = obj.arc_length_theta_and_kappa(tid);
            
            rNearest = zeros(1,nframes);
            thetaNearest = zeros(1,nframes);
            kappaNearest = zeros(1,nframes);
            dist = zeros(1,nframes);
            xNearest = zeros(1,nframes);
            yNearest = zeros(1,nframes);
            
            q = linspace(0,1); % This must be same as that used in obj.arc_length_theta_and_kappa().
            
            fittedX = obj.polyFits{ind}{1};
            fittedY = obj.polyFits{ind}{2};
            
            for k=1:nframes
                %                 disp(['Computing nearest-point values for frame=' int2str(k)])
                px = fittedX(k,:);
                py = fittedY(k,:);
                
                x = polyval(px,q);
                y = polyval(py,q);
                
%                 ind = bar_f == f(k);
                ind = abs(bar_f-f(k)) < 1e-12;
                
                if ~any(ind) || any(isnan(x)) || any(isnan(y))
                    rNearest(k) = NaN;
                    thetaNearest(k) = NaN;
                    kappaNearest(k) = NaN;
                    xNearest(k) = NaN;
                    yNearest(k) = NaN;
                    dist(k) = NaN;
                    continue
                end
                bx = bar_x(ind); % center of bar -- offset if appropriate
                by = bar_y(ind);
                
                % Find closest point iteratively:
                npoints = length(x);
                d = zeros(npoints,1);
                for j=1:npoints
                    d(j) = sqrt((x(j)-bx)^2 + (y(j)-by)^2);
                end
                ind = find(d==min(d));
                ind = ind(1); % in case there are points of equal distance, take first; Later this should be first along arc-length of whisker*****
                
                if proximity_threshold < 0
                    dist(k) = d(ind) - obj.barRadius;
                elseif d(ind) < obj.barRadius*proximity_threshold
                    thetaW = THETA{k}(ind); % Angle of whisker at last point on whisker.
                    
                    if strcmp(obj.faceSideInImage,'top')
                        a = [tand(thetaW) 1]; % Vector extending from last point on whisker along angle of whisker.
                    elseif strcmp(obj.faceSideInImage,'bottom')  % DHO, 12oct11: Don't have test data for these conditions, remains untested
                        a = [tand(thetaW) -1]; % Vector extending from last point on whisker along angle of whisker.
                    elseif strcmp(obj.faceSideInImage,'left')
                        a = [1 tand(thetaW)]; % Vector extending from last point on whisker along angle of whisker.
                    else strcmp(obj.faceSideInImage,'right') % DHO, 12oct11: Don't have test data for these conditions, remains untested
                        a = [-1 tand(thetaW)]; % Vector extending from last point on whisker along angle of whisker.
                    end
                    
                    b = [bx-x(ind) by-y(ind)]; % Vector from last point on whisker to pole center.
                    a_norm = norm(a);
                    b_norm = norm(b);
                    thetaWP = abs(acosd(dot(a,b)/(a_norm*b_norm))); % Angle between whisker and pole. Don't care about sign, so take absolute value.
                    
                    dist(k) = b_norm*sind(thetaWP);
                    dist(k) = dist(k) - obj.barRadius;
                else
                    dist(k) = d(ind) - obj.barRadius;
                end
                
                
                rNearest(k) = R{k}(ind);
                thetaNearest(k) = THETA{k}(ind);
                kappaNearest(k) = KAPPA{k}(ind);
                xNearest(k) = x(ind);
                yNearest(k) = y(ind);
            end
        end
        
        function [theta0,kappa0,t] = get_theta_kappa_at_base(obj,tid)
            %
            %   [theta0,kappa0,t] = get_theta_kappa_at_base(obj,tid)
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %
            %   OUTPUTS:
            %
            %       theta0: Theta at radial distance 0. Radial distance 0 is determined
            %                     in part by the mask, if present. In degrees.
            %
            %       kappa0: kappa at radial distance 0. Radial distance 0 is determined
            %                     in part by the mask, if present. Units of 1/pixels.
            %
            %       t: The corresponding times of each observation.
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
            
            if all(isnan(t))
                disp(['Nothing tracked for tid ' int2str(tid) '; setting theta0,kappa0,t to NaN.'])
                theta0 = NaN;
                kappa0 = NaN;
                t = NaN;
                return
            end
            
            
            f = t / obj.framePeriodInSec;
            nframes = length(f);
            
            [R,THETA,KAPPA] = obj.arc_length_theta_and_kappa(tid);
            
            theta0 = zeros(1,nframes);
            kappa0 = zeros(1,nframes);
            
            for k=1:nframes
                r = R{k};
                rval = min(r(r >= 0)); % Take the minimum value >= 0.
                if isempty(rval)
                    theta0(k) = NaN;
                    kappa0(k) = NaN;
                else
                    ind = find(r==rval,1,'first');
                    theta0(k) = THETA{k}(ind);
                    kappa0(k) = KAPPA{k}(ind);
                end
            end
        end
        
        function [thetap,kappap,y,x,t] = get_theta_kappa_at_point(obj,tid,r_in_mm)
            %
            %   [thetap,kappap,y,x,t] = get_theta_kappa_at_point(obj,tid,r_in_mm)
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %
            %       r_in_mm: Distance along whisker at which to measure theta and kappa.
            %                Note that this calculation does not extrapolate back
            %                to the follicle. Alternatively, this argument can be given
            %                as the string 'max', in which case theta, kappa, y, x are
            %                given for the furthest arc-length distance along the whisker.
            %
            %
            %
            %   OUTPUTS:
            %
            %       thetap: Theta at radial distance specified by r_in_mm.
            %                     Radial distance is measured outward from the
            %                     intersection of the whisker and the mask, if present.
            %                     In degrees.
            %
            %       kappap: Kappa at radial distance specified by r_in_mm.
            %               Radial distance is determined outward from the intersection
            %               of the whisker and the mask, if present. Units of 1/pixels.
            %
            %       x,y: The image (pixel) coordinates of the point at r_in_mm.
            %
            %       t: The corresponding times of each observation.
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
            f = t / obj.framePeriodInSec;
            nframes = length(f);
            
            [R,THETA,KAPPA,Y,X] = obj.arc_length_theta_and_kappa(tid);
            
            thetap = zeros(1,nframes);
            kappap = zeros(1,nframes);
            x = zeros(1,nframes);
            y = zeros(1,nframes);
            
            if ischar(r_in_mm)
                if ~strcmp(r_in_mm,'max')
                    error('Invalid value for argument ''r_in_mm''')
                end
                for k=1:nframes
                    if isempty(R{k})
                        thetap(k) = NaN;
                        kappap(k) = NaN;
                        y(k) = NaN;
                        x(k) = NaN;
                    else
                        thetap(k) = THETA{k}(end);
                        kappap(k) = KAPPA{k}(end);
                        y(k) = Y{k}(end);
                        x(k) = X{k}(end);
                    end
                end
            else
                for k=1:nframes
                    r = R{k} / obj.pxPerMm;
                    rval = min(r(r >= r_in_mm)); % Take the minimum value >= r_in_mm.
                    if isempty(rval)
                        thetap(k) = NaN;
                        kappap(k) = NaN;
                        y(k) = NaN;
                        x(k) = NaN;
                    else
                        ind = find(r==rval,1,'first');
                        thetap(k) = THETA{k}(ind);
                        kappap(k) = KAPPA{k}(ind);
                        y(k) = Y{k}(ind);
                        x(k) = X{k}(ind);
                    end
                end
            end
            
            
        end
        
        function [thetap,kappap,y,x,t] = get_theta_kappa_at_roi_point(obj,tid,r_in_mm)
            %
            %   [thetap,kappap,y,x,t] = get_theta_kappa_at_roi_point(obj,tid,r_in_mm)
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %
            %       r_in_mm: Distance along whisker at which to measure theta and kappa.
            %                Note that this calculation does not extrapolate back
            %                to the follicle. Also, r_in_mm **MUST FALL WITHIN THE ROI**
            %                or else a NaN is returned.
            %
            %
            %   OUTPUTS:
            %
            %       thetap: Theta at radial distance specified by r_in_mm.
            %                     Radial distance is measured outward from the
            %                     intersection of the whisker and the mask, if present.
            %                     In degrees.
            %
            %       kappap: Kappa at radial distance specified by r_in_mm.
            %               Radial distance is determined outward from the intersection
            %               of the whisker and the mask, if present. Units of 1/pixels.
            %
            %       x,y: The image (pixel) coordinates of the point at r_in_mm.
            %
            %       t: The corresponding times of each observation.
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
            f = t / obj.framePeriodInSec;
            nframes = length(f);
            
            [R,THETA,KAPPA,Y,X] = obj.arc_length_theta_and_kappa_in_roi(tid); 
            
            thetap = zeros(1,nframes);
            kappap = zeros(1,nframes);
            x = zeros(1,nframes);
            y = zeros(1,nframes);
            
            for k=1:nframes
                r = R{k} / obj.pxPerMm;
                rval = min(r(r >= r_in_mm)); % Take the minimum value >= r_in_mm.
                if isempty(rval)
                    disp('r_in_mm not found within fitted whisker ROI; setting to NaN for this frame.')
                    thetap(k) = NaN;
                    kappap(k) = NaN;
                    y(k) = NaN;
                    x(k) = NaN;
                else
                    ind = find(r==rval,1,'first');
                    thetap(k) = THETA{k}(ind);
                    kappap(k) = KAPPA{k}(ind);
                    y(k) = Y{k}(ind);
                    x(k) = X{k}(ind);
                end
            end
        end
        
        function [M0,Faxial,t,varargout] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
                whisker_length,youngs_modulus,baseline_time_or_kappa_value, varargin)
            %
            %   DHO 6/18/10 NOTE: This should be rewritten using inputParser. 
            %
            % USAGE:
            %
            %    function [M0,Faxial,t] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_or_kappa_value)
            %
            %    function [M0,Faxial,t,varargout] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_or_kappa_value)
            %
            %    function [M0,Faxial,t,deltaKappa] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_or_kappa_value)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_or_kappa_value)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm,thetaAtBase] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_or_kappa_value)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm,thetaAtBase,thetaAtContact] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_or_kappa_value)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm,thetaAtBase,thetaAtContact,distanceToPoleCenter] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_or_kappa_value)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm,thetaAtBase,thetaAtContact,distanceToPoleCenter, meanKappa] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_or_kappa_value)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm,thetaAtBase,thetaAtContact,distanceToPoleCenter, meanKappa, Flateral] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_or_kappa_value)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm,thetaAtBase,thetaAtContact,distanceToPoleCenter, meanKappa, Flateral] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_or_kappa_value,proximity_threshold)            
            %
            %
            % INPUTS:
            %
            % 	tid:  trajectory ID or string specifying the whisker the use.
            % 	r_point: radial distance along whisker at which to measure kappa. In mm.
            % 	whisker_radius_at_base: Given in microns.
            % 	whisker_length: Given in mm.
            % 	youngs_modulus: In Pa.
            % 	baseline_time_or_kappa_value: Either (1) a 1x2 vector giving starting and stopping times (inclusive) for measuring baseline whisker curvature, in sec; 
            %                                 or (2) a scaler giving a baseline kappa value (measured by the user separately) to directly subtract from kappa
            %                                 timeseries, in 1/mm.
            %
            %   Optionally,
            %       
            %   proximity_threshold:  Optional argument giving distance from nearest point on whisker
            %                         to bar center, in units of bar radius, beyond which
            %                         the whisker will be extrapolated along the last theta in
            %                         order to determine distance between whisker and bar.             
            %
            % OUTPUTS:
            %
            % 	M0:  Moment at the follicle. In Newton-meters.
            % 	Faxial: Axial force into follice. In Newtons.
            % 	t: The time of each M0, Faxial observation. In sec.
            %
            % 	Optionally,
            %
            % 	deltaKappa = Change from baseline curvature, at point specified by r_point. In 1/mm.
            % 	Fnorm - The force on the whisker normal to the contacted object. In Newtons.
            %   thetaAtBase - The whisker angle nearest the follicle. In degrees.
            %   thetaAtContact - The whisker angle nearest the point of contact. I.e., nearest the center of the pole. In degrees.
            %   distanceToPoleCenter - The closest distance between the whisker and the center of the pole. In mm.
            %   meanKappa - The mean of kappa over the entire secondary polynomial fitted ROI. In 1/mm.
            %   Flateral - Lateral force, orthogonal to Faxial; pushes the whisker against the posterior side of the follicle. In Newtons.
            %
            % REQUIRES:
            % 	-Follicle coordinates are already computed.
            % 	-Property pxPerMm is set correctly for the current videographic conditions.
            % 	-Property faceSideInImage is set correctly.
            % 	-Property barPos is set correctly.
            %
            % Assumptions:
            % 	-Whisker is conical.
            %   -Young's modulus is same everywhere on whisker.
            %   -Whisker cross-section is circular.
            %
            %
            if nargin > 7
                proximity_threshold = varargin{1};
            else
                proximity_threshold = -1; % -1 means unused
            end
            
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
            
            if isempty(obj.follicleCoordsX) || isempty(obj.follicleCoordsX)
                error(['obj.follicleCoordsX or obj.follicleCoordsY is empty. ' ...
                    'Must run obj.recompute_cached_follicle_coords before this method.'])
            end
            
            if isempty(obj.barPos)
                error('obj.barPos is empty. Bar position is required to define point of contact.')
            end
            
            r_point = r_point * 1e-3; % Argument given in mm; convert to meters.
            whisker_radius_at_base = whisker_radius_at_base * 1e-6; % Argument given in micrometers; convert to meters.
            whisker_length = whisker_length * 1e-3;  % Argument given in mm; convert to meters.
            
            pixelsPerMeter = 1e3 * obj.pxPerMm;
            
            % DHO 11/5/10: Need to access follicleExtrapDistanceInPix and then add that (mm/pixel converted) value
            % to r_point in equation for II. Otherwise the base of the cone (whisker) starts not at the follicle but
            % at the first tracked point or at the intersection of the whisker and the mask, if a mask is defined.
            II = pi/4*(whisker_radius_at_base*(1-(r_point+obj.follicleExtrapDistInPix./pixelsPerMeter)/whisker_length))^4; % In Meters^4. Eqn A2,3 from Birdwell et al.
            
%             II = pi/4*(whisker_radius_at_base*(1-r_point/whisker_length))^4; % In Meters^4. Eqn A2,3 from Birdwell et al.
            EI = youngs_modulus*II; % Meter^4-pascals. Bending stiffness at point r_point.
            
            x0 = obj.follicleCoordsX{ind};
            y0 = obj.follicleCoordsY{ind};
            
            if isempty(x0) || isempty(y0)
              disp(['No follicle coordinates computed for tid ' int2str(tid) ', cannot compute forces, setting to NaN.'])
              t = obj.time{ind};
              x0 = nan(size(t));
              y0 = nan(size(t));
            end

            
            % Could speed up by combining next few lines, since several call arc_length_theta_and_kappa().
            % From one function call, need to get:
            %   t, theta0, kappaPoint, yPoint, xPoint, meanKappa (from ROIKAPPA), thetaContact, yContact, xContact, distanceToPoleCenter 
            
            FC = obj.get_force_calc_vals(tid, r_point*1e3, proximity_threshold); % Give r_point argument in mm.
            
            % get_force_calc_vals() combines:
            %            [theta0,tmp,t] = obj.get_theta_kappa_at_base(tid);
            %            [tmp,kappaPoint,yPoint,xPoint,t] = obj.get_theta_kappa_at_roi_point(tid,r_point*1e3)
            %            [tmp1,thetaContact,tmp2,yContact,xContact,distanceToPoleCenter,t] = obj.get_r_theta_kappa_nearest_bar(tid);
            %
            % FC is structure with fields: t (sec), theta0 (degrees), kappaPoint (1/pixels), yPoint (pixels), xPoint (pixels),
            % meanKappa (in ROI; 1/pixels), thetaContact (degrees), yContact (pixels), xContact (pixels), distanceToPoleCenter (pixels).
                        % Construct output structure:
                        
            t = FC.t; % In sec
            theta0 = FC.theta0; % In degrees
            kappaPoint = FC.kappaPoint; % In 1/pixels
            yPoint = FC.yPoint; % In pixels
            xPoint = FC.xPoint; % In pixels
            meanKappa = FC.meanKappa; % In 1/pixels.
            thetaContact = FC.thetaContact; % In 1/pixels
            yContact = FC.yContact; % In pixels
            xContact = FC.xContact; % In pixels
            distanceToPoleCenter = FC.distanceToPoleCenter; % In pixels
            
%            [theta0, t] = obj.get_fitted_line_angle_in_roi(obj,tid,[0 32]) 
            [theta0,tmp,t] = obj.get_theta_kappa_at_point(tid,0.5); %  
            
%             [tmp,kappaPoint,yPoint,xPoint,t] = obj.get_theta_kappa_at_roi_point(tid,r_point*1e3); % Give second argument in mm. kappaPoint in 1/pixels.
            
            if numel(baseline_time_or_kappa_value)==2
                disp(['Using kappa averaged over times ' num2str(baseline_time_or_kappa_value(1)) ' to ' num2str(baseline_time_or_kappa_value(2)) ...
                    ' (inclusive) for baseline kappa for trajectory ID ' int2str(tid) '.'])
                baselineKappa = pixelsPerMeter * nanmean(kappaPoint(t >= baseline_time_or_kappa_value(1) & t <= baseline_time_or_kappa_value(2))); % Now in 1/m.
            elseif numel(baseline_time_or_kappa_value)==1
                disp(['Using user-specified kappa value of ' num2str(baseline_time_or_kappa_value) ' 1/mm for baseline kappa for trajectory ID ' int2str(tid) '.'])
                baselineKappa = baseline_time_or_kappa_value * 1e3; % baseline_time_or_kappa_value given in 1/mm; convert to 1/m.
            else 
                error('Argument ''baseline_time_or_kappa_value'' has wrong number of elements.')
            end
            
            kappaPoint = pixelsPerMeter * kappaPoint - baselineKappa; % In 1/m.
            
%             [tmp,tmp1,ROIKAPPA] = obj.arc_length_theta_and_kappa_in_roi(tid);
%             meanKappa = cellfun(@mean, ROIKAPPA); % In 1/pixels.
            
            meanKappa = pixelsPerMeter * meanKappa - baselineKappa;
            
%             [tmp1,thetaContact,tmp2,yContact,xContact,distanceToPoleCenter,t] = obj.get_r_theta_kappa_nearest_bar(tid);
            
            theta0 = theta0*(2*pi/360); % Convert all angles to radians.
            thetaContact = thetaContact*(2*pi/360);
            
            rPointNorm = sqrt((xContact-xPoint).^2 + (yContact-yPoint).^2) / pixelsPerMeter; % in meters.
            r0Norm = sqrt((xContact-x0).^2 + (yContact-y0).^2) / pixelsPerMeter; % in meters.
            
            %**** CHECK: NEED TO ACCOUNT FROM PROTRACTION DIRECTION IN FOLLOWING? *******
            
            % Get angle of vector from r_point to contact point:
            dx = (xContact-xPoint); dy = (yContact-yPoint);
            if strcmp(obj.faceSideInImage,'top') || strcmp(obj.faceSideInImage,'bottom')
                thetaPoint2Cont = atan(dx./dy);
            else
                thetaPoint2Cont = atan(dy./dx);
            end
            
            % Get angle of vector from follicle to contact point:
            dx = (xContact-x0); dy = (yContact-y0);
            if strcmp(obj.faceSideInImage,'top') || strcmp(obj.faceSideInImage,'bottom')
                thetaFoll2Cont = atan(dx./dy);
            else
                thetaFoll2Cont = atan(dy./dx);
            end
            
            Fnorm = (kappaPoint*EI) ./ (rPointNorm .* sin(pi/2 + thetaPoint2Cont - thetaContact)); % in Newtons
            Faxial = Fnorm.*sin(theta0 - thetaContact); % in Newtons.
            Flateral = Fnorm.*cos(theta0 - thetaContact); % in Newtons.
            M0 = r0Norm.*Fnorm.*sin(pi/2 + thetaFoll2Cont - thetaContact); % in Newton-meters
            
            if nargout > 3
                varargout{1} = kappaPoint / 1e3; % kappaPoint in 1/m; return in 1/mm.
            end
            if nargout > 4
                varargout{2} = Fnorm;
            end
            if nargout > 5
                varargout{3} = theta0 / (2*pi/360); % Whisker angle nearest the follicle. Currently in radians; return in degrees.
            end
            if nargout > 6
                varargout{4} = thetaContact / (2*pi/360); % Whisker angle at point nearest center of pole. Currently in radians; return in degrees.
            end
            if nargout > 7
                varargout{5} = distanceToPoleCenter / obj.pxPerMm; % distanceToPoleCenter in pixels; return in mm.
            end
            if nargout > 8
               varargout{6} = meanKappa / 1e3; % meanKappa in 1/m; return in 1/mm.
            end
            if nargout > 9
               varargout{7} = Flateral; % In Newtons.
            end
        end
        
        function FC = get_force_calc_vals(obj, tid, r_in_mm, varargin)
            %
            % INPUTS:
            %
            % 	tid:  trajectory ID or string specifying the whisker the use.
            % 	r_in_mm: radial distance along whisker at which to measure kappa. In mm.
            %   
            %   Optionally, varargin{1} can be 'proximity_threshold' argument, giving distance from nearest point on whisker
            %    to bar center, in units of bar radius, beyond which
            %    the whisker will be extrapolated along the last theta in
            %    order to determine distance between whisker and bar. 
            % 
            %
            % OUTPUTS:
            %
            % FC: A structure with fields: t (sec), theta0 (degrees), kappaPoint (1/pixels), yPoint (pixels), xPoint (pixels),
            % meanKappa (in ROI; 1/pixels), thetaContact (degrees), yContact (pixels), xContact (pixels), distanceToPoleCenter (pixels).
            %
            % Explanation of fields:
            %       theta0: Theta at radial distance 0. Radial distance 0 is determined
            %                     in part by the mask, if present. In degrees.
            %
            % This method combines features of the following methods, but is more efficient
            % than calling all of the following methods separately:
            %            [theta0,tmp,t] = obj.get_theta_kappa_at_base(tid);
            %            [tmp,kappaPoint,yPoint,xPoint,t] = obj.get_theta_kappa_at_roi_point(tid,r_point*1e3)
            %            [tmp1,thetaContact,tmp2,yContact,xContact,distanceToPoleCenter,t] = obj.get_r_theta_kappa_nearest_bar(tid);
            %
            %
            
            
            %----------------------------------
            
            %   [theta0,kappa0,t] = get_theta_kappa_at_base(obj,tid)
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %
            %   OUTPUTS:
            %
            %       theta0: Theta at radial distance 0. Radial distance 0 is determined
            %                     in part by the mask, if present. In degrees.
            %
            %       t: The corresponding times of each observation.
            %

            if isnumeric(tid) % Trajectory ID specified.
                ind_tid = find(obj.trajectoryIDs == tid);
            elseif ischar(tid) % Whisker name specified.
                ind_tid = strmatch(tid,obj.whiskerNames,'exact');
            else
                error('Invalid type for argument ''tid''.')
            end
            
            if isempty(ind_tid)
                error('Could not find specified trajectory ID.')
            end
            t = obj.time{ind_tid};
            
            if all(isnan(t)) %*** MODIFY?
                disp(['Nothing tracked for tid ' int2str(tid) '; setting theta0,kappa0,t to NaN.'])
                theta0 = NaN;
                t = NaN;
                return
            end
            
            f = t / obj.framePeriodInSec;
            nframes = length(f);
            
            [R,THETA] = obj.arc_length_and_theta(tid); 
            
            theta0 = zeros(1,nframes);
            
            for k=1:nframes
                r = R{k};
                rval = min(r(r >= 0)); % Take the minimum value >= 0.
                if isempty(rval)
                    theta0(k) = NaN;
                else
                    ind = find(r==rval,1,'first');
                    theta0(k) = THETA{k}(ind);
                end
            end
            %----------------------------------
            %   [thetap,kappap,y,x,t] = get_theta_kappa_at_roi_point(obj,tid,r_in_mm)
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %
            %       r_in_mm: Distance along whisker at which to measure theta and kappa.
            %                Note that this calculation does not extrapolate back
            %                to the follicle. Also, r_in_mm **MUST FALL WITHIN THE ROI**
            %                or else a NaN is returned.
            %
            %
            %   OUTPUTS:
            %
            %
            %       kappaPoint: Kappa at radial distance specified by r_in_mm.
            %               Radial distance is determined outward from the intersection
            %               of the whisker and the mask, if present. Units of 1/pixels.
            %
            %       xPoint,yPoint: The image (pixel) coordinates of the point at r_in_mm.
            %
            %       t: The corresponding times of each observation.
            %
            
            [R_ROI,THETA_ROI,KAPPA_ROI,Y,X] = obj.arc_length_theta_and_kappa_in_roi(tid);
            
            kappaPoint = zeros(1,nframes);
            xPoint = zeros(1,nframes);
            yPoint = zeros(1,nframes);
            
            for k=1:nframes
                r = R_ROI{k} / obj.pxPerMm;
                rval = min(r(r >= r_in_mm)); % Take the minimum value >= r_in_mm.
                if isempty(rval)
                    disp('r_in_mm not found within fitted whisker ROI; setting to NaN for this frame.')
                    kappaPoint(k) = NaN;
                    yPoint(k) = NaN;
                    xPoint(k) = NaN;
                else
                    ind = find(r==rval,1,'first');
                    kappaPoint(k) = KAPPA_ROI{k}(ind);
                    yPoint(k) = Y{k}(ind);
                    xPoint(k) = X{k}(ind);
                end
            end
            %----------------------------------
            %   [rNearest,thetaNearest,kappaNearest,dist,t] = get_r_theta_kappa_nearest_bar(obj,tid)
            %   [rNearest,thetaNearest,kappaNearest,dist,t] = get_r_theta_kappa_nearest_bar(obj,tid,proximity_threshold)
            %
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %       proximity_threshold:  Optional argument giving distance from nearest point on whisker
            %               to bar center, in units of bar radius, beyond which
            %               the whisker will be extrapolated along the last theta in
            %               order to determine distance between whisker and bar.
            %
            %   OUTPUTS:
            %       rNearest: Arc-length (radial) distance along whisker to point nearest
            %                 center of the bar. Units of pixels.
            %       thetaNearest: Theta at rNearest.
            %       kappaNearest: Kappa at rNearest.
            %       YNearest: Image-coordinate Y value at rNearest. This estimates point of contact. In pixels.
            %       XNearest: Image-coordinate X value at rNearest. This estimates point of contact. In pixels.
            %       dist: Distance from bar center to nearest point on whisker. Units of pixels.
            %       t: The corresponding times of each observation.
            %
            
            if nargin < 4
                proximity_threshold = -1; 
            else
                proximity_threshold = varargin{1};
            end
            
            % Add any offset to the tracked bar position:
            if isempty(obj.barPosOffset);
                bp = obj.barPosClean;
            else
                bp = obj.barPosClean;
                if size(obj.barPosOffset,1)==1
                    bp = bp + repmat([0 obj.barPosOffset],size(bp,1),1);
                elseif size(obj.barPosOffset,1)==size(bp,1) && size(obj.barPosOffset,2)==2
                    bp = bp + obj.barPosOffset;
                else
                    error('Size of obj.barPosOffset is not valid; must be either 1 x 2 or nframes x 2.')
                end
            end
            
            bar_f = bp(:,1);
            bar_x = bp(:,2);
            bar_y = bp(:,3);
            
            %             [R,THETA,KAPPA] = obj.arc_length_theta_and_kappa(tid);
            
            thetaNearest = zeros(1,nframes);
            dist = zeros(1,nframes);
            xNearest = zeros(1,nframes);
            yNearest = zeros(1,nframes);
            
            q = linspace(0,1); % This must be same as that used in obj.arc_length_theta_and_kappa().
            
            fittedX = obj.polyFits{ind_tid}{1};
            fittedY = obj.polyFits{ind_tid}{2};
            
            for k=1:nframes
                %                 disp(['Computing nearest-point values for frame=' int2str(k)])
                px = fittedX(k,:);
                py = fittedY(k,:);
                
                x = polyval(px,q);
                y = polyval(py,q);
                
                %                 ind = bar_f == f(k);
                ind = abs(bar_f-f(k)) < 1e-12;
                
                if ~any(ind) || any(isnan(x)) || any(isnan(y))
                    thetaNearest(k) = NaN;
                    xNearest(k) = NaN;
                    yNearest(k) = NaN;
                    dist(k) = NaN;
                    continue
                end
                bx = bar_x(ind); % center of bar -- offset if appropriate
                by = bar_y(ind);
                
                % Find closest point iteratively:
                npoints = length(x);
                d = zeros(npoints,1);
                for j=1:npoints
                    d(j) = sqrt((x(j)-bx)^2 + (y(j)-by)^2);
                end
                ind = find(d==min(d));
                ind = ind(1); % in case there are points of equal distance, take first; Later this should be first along arc-length of whisker*****
                
                if proximity_threshold < 0
                    dist(k) = d(ind) - obj.barRadius;
                elseif d(ind) < obj.barRadius*proximity_threshold
                    thetaW = THETA{k}(ind); % Angle of whisker at last point on whisker.
                    
                    if strcmp(obj.faceSideInImage,'top')
                        a = [tand(thetaW) 1]; % Vector extending from last point on whisker along angle of whisker.
                    elseif strcmp(obj.faceSideInImage,'bottom')  % DHO, 12oct11: Don't have test data for these conditions, remains untested
                        a = [tand(thetaW) -1]; % Vector extending from last point on whisker along angle of whisker.
                    elseif strcmp(obj.faceSideInImage,'left')
                        a = [1 tand(thetaW)]; % Vector extending from last point on whisker along angle of whisker.
                    else strcmp(obj.faceSideInImage,'right') % DHO, 12oct11: Don't have test data for these conditions, remains untested
                        a = [-1 tand(thetaW)]; % Vector extending from last point on whisker along angle of whisker.
                    end
                    
                    b = [bx-x(ind) by-y(ind)]; % Vector from last point on whisker to pole center.
                    a_norm = norm(a);
                    b_norm = norm(b);
                    thetaWP = abs(acosd(dot(a,b)/(a_norm*b_norm))); % Angle between whisker and pole. Don't care about sign, so take absolute value.
                    
                    dist(k) = b_norm*sind(thetaWP);
                    dist(k) = dist(k) - obj.barRadius;
                else
                    dist(k) = d(ind) - obj.barRadius;
                end
                
                thetaNearest(k) = THETA{k}(ind);
                xNearest(k) = x(ind);
                yNearest(k) = y(ind);
            end
            
            % Construct output structure:
            FC.t = t; % In sec
            FC.theta0 = theta0; % In degrees
            FC.kappaPoint = kappaPoint; % In 1/pixels
            FC.yPoint = yPoint; % In pixels
            FC.xPoint = xPoint; % In pixels
            FC.meanKappa = cellfun(@mean, KAPPA_ROI); % In 1/pixels.
            FC.thetaContact = thetaNearest; % In 1/pixels
            FC.yContact = yNearest; % In pixels
            FC.xContact = xNearest; % In pixels
            FC.distanceToPoleCenter = dist; % In pixels
            
        end
        
        
        function varargout = plot_fitted_whisker_time_projection(obj, tid, varargin)
            %
            %   varargout = plot_fitted_whisker_time_projection(obj, tid, varargin)
            %
            % USAGE:
            %       plot_fitted_whisker_time_projection(obj, tid, varargin)
            %       [X,Y] = plot_fitted_whisker_time_projection(obj, tid, varargin)
            %
            % INPUTS:
            %
            % tid: trajectory ID.  Only a single ID is allowed.
            %
            % varargin{1}: Plot color/symbol string.
            %
            % varargin{2}: Optional 2 x 1 vector giving starting and ending times (in seconds) to include in
            %               returned image, inclusive, starting with 0.
            %               Of format: [startTimeInSec endTimeInSec]. Can be empty ([]) to allow
            %               access to varargin{3}.
            %
            % varargin{3}: Cell array with two elements. First element is a 1x2 vector giving radial distance
            %              beginning and end points (inclusive) to draw in the color given in
            %              varargin{1}, in format [startRadialDistance stopRadialDistance]. Radial distance
            %              is whisker arc length moving outward from follicle and is in units of pixels.
            %              The second element of the cell array is the plot color/symbol string to use when
            %              plotting the radial segment of the whisker between startRadialDistance and
            %              stopRadialDistance.
            %
            % OUTPUTS:
            %
            % If no output argument is requested, only plotting happens.  If two output argument
            % is requested, two cell arrays of length nframes are returned, X, and Y, each element of which
            % contain the x and y data, respectively, that were plotted for a given frame.
            %
            %
            if numel(tid) > 1
                error('Only a single trajectory ID is alowed.')
            end
            
            if nargin > 2
                plotString = varargin{1};
            else
                plotString = 'k-';
            end
            
            t = obj.get_time(tid);
            
            if nargin > 3
                restrictTime = varargin{2};
                if isempty(restrictTime)
                    restrictTime = [min(t) max(t)];
                elseif restrictTime(2) <= restrictTime(1)
                    error('Invalid format for varargin{2}.')
                elseif max(restrictTime) > max(t)
                    disp('Warning: varargin{2} exceeds max time; setting to max.')
                    restrictTime(restrictTime==max(restrictTime)) = max(t);
                    if restrictTime(1)==restrictTime(2)
                        error('varargin{2}: Both times exceed max time.')
                    end
                elseif min(restrictTime < 0)
                    disp('varargin{2}: times start at 0.')
                end
            else
                restrictTime = [min(t) max(t)];
            end
            
            if nargin > 4
                if ~iscell(varargin{3})
                    error('Wrong format for varargin{3}.')
                end
                plotString2 = varargin{3}{2};
                rstart = varargin{3}{1}(1);
                rstop = varargin{3}{1}(2);
            end
            
            frameInds = find(t >= restrictTime(1) & t <= restrictTime(2));
            
            hold on; axis ij
            
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Could not find specified trajectory ID.')
            end
            
            if isempty(obj.polyFits)
                error('obj.polyFits is empty.')
            end
            
            fittedX = obj.polyFits{ind}{1};
            fittedY = obj.polyFits{ind}{2};
            
            nframes = length(frameInds);
            x = cell(1,nframes);
            y = cell(1,nframes);
            
            q = linspace(0,1);
            
            for k=1:nframes
                f = frameInds(k);
                
                px = fittedX(f,:);
                py = fittedY(f,:);
                x{k} = polyval(px,q);
                y{k} = polyval(py,q);
            end
            
            X = cell(1,nframes);
            Y = cell(1,nframes);
            
            if nargin > 4
                [R,THETA,KAPPA] = obj.arc_length_theta_and_kappa(tid);
                for k=1:nframes
                    xx = x{k};
                    yy = y{k};
                    f = frameInds(k);
                    radDist = R{f};
                    plot(xx,yy,plotString)
                    rind = radDist >= rstart & radDist <= rstop;
                    plot(xx(rind),yy(rind),plotString2)
                    X{k} = xx(rind);
                    Y{k} = yy(rind);
                end
            else
                for k=1:nframes
                    plot(x{k},y{k},plotString)
                    X{k} = x{k};
                    Y{k} = y{k};
                end
            end
            
            if nargout == 1 || nargout > 2
                error('Invalid number of output arguments.')
            elseif nargout == 2
                varargout{1} = X;
                varargout{2} = Y;
            end
            
        end
        
        function plot_fitted_whisker_ROI_time_projection(obj, tid, varargin)
            %
            %   plot_fitted_whisker_time_ROI_projection(obj, tid, varargin)
            %
            % Plots whisker, but only the part within the secondary polynomial
            % fitting ROI.
            %
            % tid: trajectory ID.  Only a single ID is allowed.
            %
            % varargin{1}: Plot color/symbol string.
            %
            % varargin{2}: Optional 2 x 1 vector giving starting and ending times (in seconds) to include in
            %               returned image, inclusive, starting with 0.
            %               Of format: [startTimeInSec endTimeInSec]. Can be empty ([]) to allow
            %               access to varargin{3}.
            %
            %
            %
            if numel(tid) > 1
                error('Only a single trajectory ID is alowed.')
            end
            
            if nargin > 2
                plotString = varargin{1};
            else
                plotString = 'k-';
            end
            
            t = obj.get_time(tid);
            
            if nargin > 3
                restrictTime = varargin{2};
                if isempty(restrictTime)
                    restrictTime = [min(t) max(t)];
                elseif restrictTime(2) <= restrictTime(1)
                    error('Invalid format for varargin{2}.')
                elseif max(restrictTime) > max(t)
                    disp('Warning: varargin{2} exceeds max time; setting to max.')
                    restrictTime(restrictTime==max(restrictTime)) = max(t);
                    if restrictTime(1)==restrictTime(2)
                        error('varargin{2}: Both times exceed max time.')
                    end
                elseif min(restrictTime < 0)
                    disp('varargin{2}: times start at 0.')
                end
            else
                restrictTime = [min(t) max(t)];
            end
            
            frameInds = find(t >= restrictTime(1) & t <= restrictTime(2));
            
            hold on; axis ij
            
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Could not find specified trajectory ID.')
            end
            
            if isempty(obj.polyFitsROI)
                error('obj.polyFitsROI is empty.')
            end
            
            fittedX = obj.polyFitsROI{ind}{1};
            fittedY = obj.polyFitsROI{ind}{2};
            fittedQ = obj.polyFitsROI{ind}{3};
            
            nframes = length(frameInds);
            x = cell(1,nframes);
            y = cell(1,nframes);
                       
            for k=1:nframes
                f = frameInds(k);
                
                px = fittedX(f,:);
                py = fittedY(f,:);
                pq = fittedQ(f,:);
                
                q = linspace(pq(1),pq(2));
                
                x{k} = polyval(px,q);
                y{k} = polyval(py,q);
                
                plot(x{k},y{k},plotString)
            end
        end
        
        function plot_fitted_whisker_movie(obj, tid, varargin)
            %
            %   plot_fitted_whisker_movie(obj, tid, varargin)
            %
            % tid: trajectory ID.  Only a single ID is allowed.
            %
            % varargin{1}: Plot color/symbol string.
            %
            % varargin{2}: Optional 2 x 1 vector giving starting and ending times (in seconds) to include in
            %               returned image, inclusive, starting with 0.
            %               Of format: [startTimeInSec endTimeInSec]. Can be empty ([]) to allow
            %               access to varargin{3}.
            %
            % varargin{3}: Cell array with two elements. First element is a 1x2 vector giving radial distance
            %              beginning and end points (inclusive) to draw in the color given in
            %              varargin{1}, in format [startRadialDistance stopRadialDistance]. Radial distance
            %              is whisker arc length moving outward from follicle and is in units of pixels.
            %              The second element of the cell array is the plot color/symbol string to use when
            %              plotting the radial segment of the whisker between startRadialDistance and
            %              stopRadialDistance.
            %
            %
            %
            if numel(tid) > 1
                error('Only a single trajectory ID is alowed.')
            end
            
            if nargin > 2
                plotString = varargin{1};
            else
                plotString = 'k-';
            end
            
            t = obj.get_time(tid);
            
            if nargin > 3
                restrictTime = varargin{2};
                if isempty(restrictTime)
                    restrictTime = [min(t) max(t)];
                elseif restrictTime(2) <= restrictTime(1)
                    error('Invalid format for varargin{2}.')
                elseif max(restrictTime) > max(t)
                    disp('Warning: varargin{2} exceeds max time; setting to max.')
                    restrictTime(restrictTime==max(restrictTime)) = max(t);
                    if restrictTime(1)==restrictTime(2)
                        error('varargin{2}: Both times exceed max time.')
                    end
                elseif min(restrictTime < 0)
                    disp('varargin{2}: times start at 0.')
                end
            else
                restrictTime = [min(t) max(t)];
            end
            
            if nargin > 4
                if ~iscell(varargin{3})
                    error('Wrong format for varargin{3}.')
                end
                plotString2 = varargin{3}{2};
                rstart = varargin{3}{1}(1);
                rstop = varargin{3}{1}(2);
            end
            
            frameInds = find(t >= restrictTime(1) & t <= restrictTime(2));
            
            hold on; axis ij
            
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Could not find specified trajectory ID.')
            end
            
            if isempty(obj.polyFits)
                error('obj.polyFits is empty.')
            end
            
            fittedX = obj.polyFits{ind}{1};
            fittedY = obj.polyFits{ind}{2};
            
            nframes = length(frameInds);
            x = cell(1,nframes);
            y = cell(1,nframes);
            
            q = linspace(0,1);
            
            for k=1:nframes
                f = frameInds(k);
                px = fittedX(f,:);
                py = fittedY(f,:);
                x{k} = polyval(px,q);
                y{k} = polyval(py,q);
            end
            
            if nargin > 4
                [R,THETA,KAPPA] = obj.arc_length_theta_and_kappa(tid);
                for k=1:nframes
                    xx = x{k};
                    yy = y{k};
                    f = frameInds(k);
                    radDist = R{f};
                    plot(xx,yy,plotString)
                    rind = radDist >= rstart & radDist <= rstop;
                    plot(xx(rind),yy(rind),plotString2)
                    title(['t=' num2str(t(k)) ' s'])
                    pause(.1)
                    cla
                end
            else
                for k=1:nframes
                    plot(x{k},y{k},plotString)
                    title(['t=' num2str(t(k)) ' s'])
                    pause(.1)
                    cla
                end
            end
        end
        
        function plot_follicle_position_time_projection(obj, tid, varargin)
            %
            %   plot_follicle_position_time_projection(obj, tid, varargin)
            %
            % tid: trajectory ID.  Only a single ID is allowed.
            %
            % varargin{1}: Plot color/symbol string.
            %
            % varargin{2}: Optional 2 x 1 vector giving starting and ending times (in seconds) to include in
            %               returned image, inclusive, starting with 0.
            %               Of format: [startTimeInSec endTimeInSec]. Can be empty ([]) to allow
            %               access to varargin{3}.
            %
            %
            %
            if numel(tid) > 1
                error('Only a single trajectory ID is alowed.')
            end
            
            if nargin > 2
                plotString = varargin{1};
            else
                plotString = 'k-';
            end
            
            [y,x,t] = obj.get_cached_follicle_coords(tid);
            
            if nargin > 3
                restrictTime = varargin{2};
                if isempty(restrictTime)
                    restrictTime = [min(t) max(t)];
                elseif restrictTime(2) <= restrictTime(1)
                    error('Invalid format for varargin{2}.')
                elseif max(restrictTime) > max(t)
                    disp('Warning: varargin{2} exceeds max time; setting to max.')
                    restrictTime(restrictTime==max(restrictTime)) = max(t);
                    if restrictTime(1)==restrictTime(2)
                        error('varargin{2}: Both times exceed max time.')
                    end
                elseif min(restrictTime < 0)
                    disp('varargin{2}: times start at 0.')
                end
            else
                restrictTime = [min(t) max(t)];
            end
            
            frameInds = find(t >= restrictTime(1) & t <= restrictTime(2));
            hold on; axis ij
            plot(x(frameInds),y(frameInds),plotString)
            
        end
        
        
        
        
    end
    methods % Dependent property methods; cannot have attributes.
        
        function value = get.barPosClean(obj)
            % Bar position cleaned up by processing to handle bar tracker errors/limitations.
            
            if isempty(obj.barPos)
                value = [];
            else
                threshInMm = 2; % Values exceeding this distance (in mm) from the
                % modal position will be set to the mode.
                
                f = obj.barPos(:,1);
                x = round(obj.barPos(:,2));
                y = round(obj.barPos(:,3));
                ymode = mode(y);
                xmode = mode(x);
                
                ind = abs((y-ymode)/obj.pxPerMm) > threshInMm;
                
                % Replace outliers with mode. Note that this mode
                % has been rounded to the nearest pixel, whereas original
                % data in obj.barPos was sub-pixel resolution.
                y(ind) = ymode;
                x(ind) = xmode;
                
                value = [f x y];
            end
        end
        
    end
    
end
