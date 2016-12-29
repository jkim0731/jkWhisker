classdef WhiskerSignalTrialArray < handle
%
%   WhiskerSignalTrialArray < handle
%
%
% DHO, 8/08.
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

        function obj = WhiskerSignalTrialArray(w, varargin) 
            %
            % USAGE:
            %
            %   obj = WhiskerSignalTrialArray(w)
            %   obj = WhiskerSignalTrialArray(w, 'polyRoiInPix',[roiMin, roiMax])
            %   obj = WhiskerSignalTrialArray(w, 'polyRoiInPix',{[trajectoryIDs],[roiMin1, roiMax1], ...
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
            %   If WhiskerSignalTrialArray(w) is called without 'polyRoiInPix'
            %   argument, and w.fit_polys_roi has not been already called on 
            %   the WhiskerTrialArray argument, then
            %   an error is given.  A WhiskerSignalTrialArray requires
            %   polynomials to be fitted to the whiskers.
            %
            %
            p = inputParser;
            p.addOptional('w', @(x) isa(x,'Whisker.WhiskerTrialArray'));                      
            p.addParamValue('polyRoiInPix', NaN); 
%             p.parse(varargin{:});
            
            %-------------
            if nargin==0
                return
            end
                        
            obj.mouseName = w.mouseName;
            obj.sessionName = w.sessionName;

            ntrials = length(w);

            obj.trials = cell(1, ntrials);
            for k=1:ntrials
                disp(['Trial = ' int2str(k)])
                obj.trials{k} = Whisker.WhiskerSignalTrial(w.trials{k}, varargin{:}); 
            end

        end

        function r = length(obj)
            %
            %   r = length(obj)
            %
            r = length(obj.trials);
        end
        
        function set_all_whiskerNames(obj, whisker_names)
            %
            %   set_all_whiskerNames(obj, whisker_names)
            %
            % whisker_names: Cell array of whisker names, e.g., {'D4','D3','D2'}.
            %
            % Sets whiskerNames property for all WhiskerSignalTrials 
            % in this WhiskerSignalTrialArray.
            %
            if ~iscell(whisker_names)
                error('Argument whisker_names must be a cell array of names.')
            end
            
            ntrials = length(obj);
            for k=1:ntrials
                obj.trials{k}.whiskerNames = whisker_names;
            end
        end
        
        function obj = set_pxPerMm(obj, px_per_mm)
            %
            %   obj = set_pxPerMm(obj, px_per_mm)
            %
            % INPUTS:
            %    A scaler giving the number of pixels per mm for
            %    the videographic conditions of these trials.
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
        
        function set_all_maskTreatment(obj,mask_treatment)
            %
            %   set_all_maskTreatments(obj,mask_treatment)
            %
            %   maskTreatment: String describing treatment of mask. 
            %
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
            %
            % Sets maskTreatment property for all WhiskerSignalTrials
            % in this WhiskerSignalTrialArray.
            %
            if ~ischar(mask_treatment)
                error('Argument mask_treatment must be a string.')
            end
            
            ntrials = length(obj);
            for k=1:ntrials
                obj.trials{k}.maskTreatment = mask_treatment;
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
            % If dx,dy each has as many elements as this WhiskerSignalTrialArray,
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
        
        function obj = set_mask_from_points(obj,tid,x,y)
            %
            % For each WhiskerSignalTrial in this WhiskerSignalTrialArray,
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
        
        function [M0,Faxial,t,varargout] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
                whisker_length,youngs_modulus,baseline_time_end)
            %
            % USAGE:
            %
            %    function [M0,Faxial,t] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_end)
            %
            %    function [M0,Faxial,t,varargout] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_end)
            %
            %    function [M0,Faxial,t,deltaKappa] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_end)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_end)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm,thetaAtBase] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_end)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm,thetaAtBase,thetaAtContact] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_end)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm,thetaAtBase,thetaAtContact,distanceToPoleCenter] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_end)
            %
            %    function [M0,Faxial,t,deltaKappa,Fnorm,thetaAtBase,thetaAtContact,distanceToPoleCenter,meanKappa] = calc_M0_Faxial(obj,tid,r_point,whisker_radius_at_base,...
            %      whisker_length,youngs_modulus,baseline_time_end)
            %
            %
            % INPUTS:
            %
            % 	tid:  trajectory ID or string specifying the whisker the use.
            % 	r_point: radial distance along whisker at which to measure kappa. In mm.
            % 	whisker_radius_at_base: Given in microns.
            % 	whisker_length: Given in mm.
            % 	youngs_modulus: In Pa.
            % 	baseline_time_end: In sec. For measuring baseline whisker curvature.
            %
            % OUTPUTS (Each a cell array, with elements comprising the following for one trial):
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
            %   meanKappa - The mean of kappa over the entire secondary polynomial fitted ROI.
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
            if ~isempty(obj.trials)
                ntrials = length(obj);
                M0 = cell(1,ntrials);
                Faxial = cell(1,ntrials);
                t = cell(1,ntrials);
                deltaKappa = cell(1,ntrials);
                Fnorm = cell(1,ntrials);
                thetaAtBase = cell(1,ntrials);
                thetaAtContact = cell(1,ntrials);
                distanceToPoleCenter = cell(1,ntrials);
                meanKappa = cell(1,ntrials);
                for k=1:ntrials
                    disp(['Trial = ' int2str(k) '/' int2str(ntrials)])
                    [M0{k},Faxial{k},t{k},deltaKappa{k},Fnorm{k},thetaAtBase{k},thetaAtContact{k},distanceToPoleCenter{k},meanKappa{k}] = ...
                        obj.trials{k}.calc_M0_Faxial(tid,r_point,whisker_radius_at_base,...
                        whisker_length,youngs_modulus,baseline_time_end);
                end
            else
                M0 = {};
                Faxial = {};
                t = {};
                deltaKappa = {};
                Fnorm = {};
                thetaAtBase = {};
                thetaAtContact = {};
                distanceToPoleCenter = {};
                meanKappa = {};
            end
            
            if nargout > 3
                varargout{1} = deltaKappa;
            end
            if nargout > 4
                varargout{2} = Fnorm;
            end
            if nargout > 5
                varargout{3} = thetaAtBase;
            end
            if nargout > 6
                varargout{4} = thetaAtContact;
            end
            if nargout > 7
                varargout{5} = distanceToPoleCenter;
            end
            if nargout > 8
                varargout{6} = meanKappa;
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
            if ~isempty(obj.trials)
                ntrials = length(obj);
                if nargin==1
                    for k=1:ntrials
                        obj.trials{k}.recompute_cached_mean_theta_kappa;
                    end
                else
                    for k=1:ntrials
                        disp(['Trial = ' int2str(k) '/' int2str(ntrials)])
                        obj.trials{k}.recompute_cached_mean_theta_kappa(varargin{:});
                    end
                end
            end
        end
        
        function [Y,T] = get_velocity(obj,tid,varargin)
            %
            %   [Y,T] = get_velocity(obj,tid,varargin)
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
            %
            % Returns two cell arrays: Y, each
            % element of which contains the velocity
            % time series for a trial (in degrees/sec);
            % and T, which contains the times of each
            % sample in Y.
            %
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                Y = cell(1,ntrials); T = cell(1,ntrials);
                if nargin > 2
                    for k=1:length(obj.trials)
                        [y,t] = obj.trials{k}.get_velocity(tid,varargin{1});
                        Y{k} = y; T{k} = t;
                    end
                else
                    for k=1:length(obj.trials)
                        [y,t] = obj.trials{k}.get_velocity(tid);
                        Y{k} = y; T{k} = t;
                    end
                end
            else
                Y = {}; T = {};
            end
        end
        
        function [Y,T] = get_velocity_medfilt(obj,tid,varargin)
            %
            %   [Y,T] = get_velocity_medfilt(obj,tid,varargin)
            %
            % Angular velocity in degrees per second, after filtering
            % position signal with a median filter.
            %
            % tid: Trajectory ID.
            %
            %
            % varargin{1}: Optional smoothing window, in frames,
            %              **for position (theta) signal.** Velocity is not
            %              separately filtered. May want to filter theta to
            %              eliminate noise due to whisker tracking artifacts.
            %              Default is 3. Should be
            %               an odd number (see help medfilt1).
            %
            %
            % Returns two cell arrays: Y, each
            % element of which contains the velocity
            % time series for a trial (in degrees/sec);
            % and T, which contains the times of each
            % sample in Y.
            %
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                Y = cell(1,ntrials); T = cell(1,ntrials);
                if nargin > 2
                    for k=1:length(obj.trials)
                        [y,t] = obj.trials{k}.get_velocity_medfilt(tid,varargin{1});
                        Y{k} = y; T{k} = t;
                    end
                else
                    for k=1:length(obj.trials)
                        [y,t] = obj.trials{k}.get_velocity_medfilt(tid);
                        Y{k} = y; T{k} = t;
                    end
                end
            else
                Y = {}; T = {};
            end
        end
        
        function [Y,T] = get_position(obj,tid,varargin)
            %
            %   [Y,T] = get_position(obj,tid,varargin)
            %
            %
            % tid: Trajectory ID.
            %
            % varargin{1}: Optional smoothing span, in frames.
            %               If not specified there is no smoothing. Should be
            %               an odd number (see 'help smooth').
            %
            % Returns two cell arrays: Y, each
            % element of which contains the angular position
            % time series for a trial (in degrees);
            % and T, which contains the times of each
            % sample in Y.
            %
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                Y = cell(1,ntrials); T = cell(1,ntrials);
                if nargin > 2
                    for k=1:length(obj.trials)
                        [y,t] = obj.trials{k}.get_position(tid,varargin{1});
                        Y{k} = y; T{k} = t;
                    end
                else
                    for k=1:length(obj.trials)
                        [y,t] = obj.trials{k}.get_position(tid);
                        Y{k} = y; T{k} = t;
                    end
                end
            else
                Y = {}; T = {};
            end
        end
        
       function [Y,X,T] = get_follicle_coords(obj,tid,extrap_distance,varargin)
            %
            %    [Y,X,T] = get_y0(obj,tid,extrap_distance,varargin)
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
            %        whisker does not cross the maks.  Otherwise
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
            % Y: A cell array where each element contains the y-coordinates 
            %    in image pixel units of the follicle for one trial. Y has
            %    one element for each trial.
            %
            % X: A cell array where each element contains the x-coordinates 
            %    in image pixel units of the follicle for one trial. X has
            %    one element for each trial.
            %
            % T: The corresponding times.
            %
            %
            %
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                Y = cell(1,ntrials); X = cell(1,ntrials); T = cell(1,ntrials);
                for k=1:length(obj.trials)
                    disp(k)
                    if nargin > 3
                        [y,x,t] = obj.trials{k}.get_follicle_coords(tid,extrap_distance,varargin{:});
                    else
                        [y,x,t] = obj.trials{k}.get_follicle_coords(tid,extrap_distance);
                    end
                    Y{k} = y; X{k} = x; T{k} = t;
                end
            else
                Y = {}; X = {}; T = {};
            end
        end

        function obj = recompute_cached_follicle_coords(obj, extrap_distance, varargin)
            %
            %  obj = recompute_cached_follicle_coords(obj, extrap_distance, varargin)
            %  obj = recompute_cached_follicle_coords(obj, extrap_distance)
            %  obj = recompute_cached_follicle_coords(obj, extrap_distance, tidList)
            %  obj = recompute_cached_follicle_coords(obj, extrap_distance, tidList, roi)
            %
            % Recomputes obj.follicleCoords for each WhiskerSignalTrial in this WhiskerSignalTrialArray.
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
            
            if extrap_distance <= 0 || length(extrap_distance) > 1
                error('Argument extrap_distance must be a positive scaler.')
            end
            
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                for k=1:length(obj.trials)
                    disp(['Trial ' int2str(k) '/' int2str(ntrials)])
                    if nargin > 2
                        obj.trials{k} = obj.trials{k}.recompute_cached_follicle_coords(extrap_distance, varargin{:});
                    else
                        obj.trials{k} = obj.trials{k}.recompute_cached_follicle_coords(extrap_distance);
                    end
                end
            end
        end
        

        function [Y,X,T] = get_cached_follicle_coords(obj,tid)
            %
            %    [Y,X,T] = get_cached_follicle_coords(obj,tid)
            %
            %  INPUTS:
            %
            %   tid: Trajectory ID.
            %
            %
            %  RETURNS:
            %
            % Y: A cell array where each element gives y for a trial, 
            %    where y is the y coordinate in image pixels of the follicle for each 
            %    time point (i.e. frame). 
            %
            % X: A cell array where each element gives x for a trial, 
            %    where x is the x coordinate in image pixels of the follicle for each 
            %    time point (i.e. frame). 
            %
            % T: The corresponding times of each observation in X,Y.
            %
            %
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                Y = cell(1,ntrials); X = cell(1,ntrials);  T = cell(1,ntrials);
                for k=1:length(obj.trials)
                    [y,x,t] = obj.trials{k}.get_cached_follicle_coords(tid);
                    Y{k} = y; X{k} = x; T{k} = t;
                end
            else
                Y = {}; X = {}; T = {};
            end
        end
        
        function [Y,T] = get_follicle_translation(obj,tid)
            %
            %    [Y,T] = get_follicle_translation(obj,tid)
            %
            %  INPUTS:
            %
            %   tid: Trajectory ID.
            %
            %
            %  RETURNS:
            %
            % Y: A cell array where each element gives y for a trial, 
            %    where y is the distance 
            %    the follicle has translated from the previous
            %    frame. Units of pixels.
            %
            % T: The corresponding times of each observation.
            %
            %
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                Y = cell(1,ntrials); T = cell(1,ntrials);
                for k=1:length(obj.trials)
                    [y,t] = obj.trials{k}.get_follicle_translation(tid);
                    Y{k} = y; T{k} = t;
                end
            else
                Y = {}; T = {};
            end
        end
        
        function [Y,T] = get_mean_position(obj)
            %
            %   [Y,T] = get_mean_position(obj)
            %
            %  INPUTS:
            %
            %   None.
            %
            %
            %  RETURNS:
            %
            % Y: A cell array where each element gives y for a trial,where
            %   y is the ***mean position of all whiskers*** (tids) in trial.
            %   This is for use, e.g., in determining whether there is overall
            %   whisking after fully-automated tracking.
            %
            %   T: The corresponding times of each observation.
            %
            %
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                Y = cell(1,ntrials); T = cell(1,ntrials);
                for k=1:length(obj.trials)
                    [y,t] = obj.trials{k}.get_mean_position;
                    Y{k} = y; T{k} = t;
                end
            else
                Y = {}; T = {};
            end
        end
        
        function [THETA0,KAPPA0,T] = get_theta_kappa_at_base(obj,tid)
            %
            %   [THETA0,KAPPA0,T] = get_theta_kappa_at_base(obj,tid)
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %
            %   OUTPUTS:
            %
            %       THETA0: Cell array where each elements gives for one trial
            %               theta at radial distance 0. Radial distance 0 is determined
            %               in part by the mask, if present. In degrees.
            %
            %       KAPPA0: Cell array where each elements gives for one trial
            %               kappa at radial distance 0. Radial distance 0 is determined
            %               in part by the mask, if present. Units of 1/pixels.
            %
            %       T: The corresponding times of each observation.
            %
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                THETA0 = cell(1,ntrials);
                KAPPA0 = cell(1,ntrials);
                T = cell(1,ntrials);
                for k=1:length(obj.trials)
                    [THETA0{k},KAPPA0{k},T{k}] = obj.trials{k}.get_theta_kappa_at_base(tid);
                end
            else
                THETA0 = {}; 
                KAPPA0 = {}; 
                T = {}; 
            end
        end
        
        function [THETAP,KAPPAP,T] = get_theta_kappa_at_point(obj,tid,r_in_mm)
            %
            %   [THETAP,KAPPAP,T] = get_theta_kappa_at_point(obj,tid,r_in_mm)
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %
            %       r_in_mm: Distance along whisker at which to measure theta and kappa.
            %                Note that this calculation does not extrapolate back  
            %                to the follicle.
            %
            %
            %   OUTPUTS:
            %
            %       THETAP: Cell array where each element is for one trial
            %                     theta at radial distance specified by r_in_mm.
            %                     Radial distance is measured outward from the
            %                     intersection of the whisker and the mask, if present.
            %                     In degrees.
            %
            %       KAPPAP: Cell array where each element is for one trial
            %               kappa at radial distance specified by r_in_mm.
            %               Radial distance is determined outward from the intersection
            %               of the whisker and the mask, if present. Units of 1/pixels.
            %
            %       T: The corresponding times of each observation.
            %
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                THETAP = cell(1,ntrials);
                KAPPAP = cell(1,ntrials);
                T = cell(1,ntrials);
                for k=1:length(obj.trials)
                    [THETAP{k},KAPPAP{k},T{k}] = obj.trials{k}.get_theta_kappa_at_point(tid,r_in_mm);
                end
            else
                THETAP = {};
                KAPPAP = {};
                T = {};
            end
        end
        
       function [R_NEAREST,THETA_NEAREST,KAPPA_NEAREST,Y_NEAREST,X_NEAREST,DIST,T] = get_r_theta_kappa_nearest_bar(obj,tid,proximity_threshold)
            %
            %   [rNearest,thetaNearest,kappaNearest,dist,t] =
            %           get_r_theta_kappa_nearest_bar(obj,tid,proximity_threshold)
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %       proximity_threshold: Optional argument giving distance from nearest point on whisker
            %               to bar center, in units of bar radius, beyond which
            %               the whisker will be extrapolated along the last theta in
            %               order to determine distance between whisker and
            %               bar. 
            %
            %   INPUTS:
            %       tid:  Trajectory ID (as an integer) or whisker name (as a string).
            %
            %   OUTPUTS:
            %       R_NEAREST: Cell array where each element gives for one trial the arc-length 
            %                 (radial) distance along whisker to point nearest
            %                 center of the bar. Units of pixels.
            %       THETA_NEAREST: Cell array where each element gives for one trial theta at rNearest.
            %       KAPPA_NEAREST: Cell array where each element gives for one trial kappa at rNearest.
            %       Y_NEAREST: Cell array where each element gives for one trial image y-coordinate at rNearest. In pixels.
            %       X_NEAREST: Cell array where each element gives for one trial image x-coordinate at rNearest. In pixels.
            %       DIST: Cell array where each element gives for one trial distance from bar center 
            %             to nearest point on whisker. Units of pixels.
            %       T: Cell array containing the corresponding times.
            %
      
            if ~isempty(obj.trials)
                ntrials = length(obj.trials);
                R_NEAREST = cell(1,ntrials);
                THETA_NEAREST = cell(1,ntrials);
                KAPPA_NEAREST = cell(1,ntrials);
                Y_NEAREST = cell(1,ntrials);
                X_NEAREST = cell(1,ntrials);
                DIST = cell(1,ntrials);
                T = cell(1,ntrials);
                if nargin < 3
                    for k=1:length(obj.trials)
                        [R_NEAREST{k},THETA_NEAREST{k},KAPPA_NEAREST{k},Y_NEAREST, X_NEAREST, DIST{k},T{k}] = obj.trials{k}.get_r_theta_kappa_nearest_bar(tid);
                    end
                else
                    for k=1:length(obj.trials)
                        [R_NEAREST{k},THETA_NEAREST{k},KAPPA_NEAREST{k},Y_NEAREST, X_NEAREST, DIST{k},T{k}] = obj.trials{k}.get_r_theta_kappa_nearest_bar(tid,proximity_threshold);
                    end
                end
            else
                R_NEAREST = {};
                THETA_NEAREST = {};
                KAPPA_NEAREST = {};
                Y_NEAREST = {};
                X_NEAREST = {};
                DIST = {};
                T = {};
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

