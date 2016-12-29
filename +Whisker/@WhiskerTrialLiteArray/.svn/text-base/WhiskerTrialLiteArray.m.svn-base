classdef WhiskerTrialLiteArray < handle
    %
    %   WhiskerTrialLiteArray < handle
    %
    % Contains an array of WhiskerTrialLite objects, each of which contains
    % just final timeseries measurements from whiskers (kappa, theta, etc)
    % and meta-information in order to match with other data (ephys, imaging, etc).
    %
    %
    % type 'help Whisker.WhiskerTrialLiteArray.WhiskerTrialLiteArray' for constructer usage.
    %
    % DHO, 3/10.
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
        
        function obj = WhiskerTrialLiteArray(w_or_d,varargin)
            %
            %-------------
            % USAGE:
            %-------------
            % There are two basic uses of this constructor. In the first case the
            % first argument is a WhiskerSignalTrialArray object. In the second case, the
            % first argument is a string giving the path of a directory containing
            % WhiskerTrialLite objects in individual .mat files. These .mat
            % files must have the suffix/extension '_WL.mat'.
            %
            %   When giving a WhiskerSignalTrialArray as input:
            %
            %   obj = WhiskerTrialLiteArray(w)
            %   obj = WhiskerTrialLiteArray(w, 'polyRoiInPix',[roiMin, roiMax])
            %   obj = WhiskerTrialLiteArray(w, 'polyRoiInPix',{[trajectoryIDs],[roiMin1, roiMax1], ...
            %                                     [roiMin2, roiMax2], ...
            %                                     [roiMin3, roiMax3], ...
            %                                     [roiMinN, roiMaxN]})
            %
            %   obj = WhiskerTrialLiteArray(w, 'polyRoiInPix',[roiMin, roiMax],'r_in_mm',some_number)
            %   obj = WhiskerTrialLiteArray(w, 'r_in_mm',some_number)
            %
            %   When giving a file list of .mat files:
            %
            %   obj = WhiskerTrialLiteArray(d)
            %   obj = WhiskerTrialLiteArray(d,'include_files',{'filename1','filename2','filenameN'})
            %   obj = WhiskerTrialLiteArray(d,'ignore_files',{'filename1','filename2','filenameN'})
            %
            %----------------------------------------------------------------
            % INPUTS - For use when a WhiskerSignalTrialArray is given:
            %----------------------------------------------------------------
            %   w: a WhiskerSignalTrialArray object.
            %
            %   Optional parameter/value pair arguments:
            %
            %           'r_in_mm': The arc-length along whisker at which to measure kappa. Units of mm. Defaults to 1 mm.
            %
            %           'calc_forces': Either true or false. Requires the pole position
            %                   to be tracked (i.e., barPos property of WhiskerSignalTrial must
            %                   not be empty). Default is false. If true, will calculate the following timeseries:
            %                       -M0:  Moment at the follicle. In Newton-meters.
            %                       -Faxial: Axial force into follice. In Newtons.
            %                       -deltaKappa: Change from baseline curvature, at point specified by r_point. In 1/mm.
            %                       -Fnorm: The force on the whisker normal to the contacted object. In Newtons.
            %                       -thetaAtBase: The whisker angle nearest the follicle. In degrees.
            %                       -thetaAtContact: The whisker angle nearest the point of contact. I.e., nearest the center of the pole. In degrees.
            %                       -distanceToPoleCenter: The closest distance between the whisker and the center of the pole. In mm.
            %                       -meanKappa: The mean of kappa over the entire secondary polynomial fitted ROI. In 1/mm.
            %
            %   The following optional parameter/value pair arguments are ignored if 'calc_forces'
            %   is not true:
            %
            %           'whisker_radius_at_base': Given in microns. Defaults is 33.5 microns.
            %
            %           'whisker_length': Given in mm. Default is 16 mm.
            %
            %           'youngs_modulus': In Pa. Default is 5e9 Pa.
            %
            %           'baseline_time_end': In sec. For measuring baseline whisker curvature. Default is 0.1 s.
            %
            %   Note: Still need make these arguments settable on a whisker-by-whisker (trajectory ID by trajectory ID) basis.
            %
            %----------------------------------------------------------------
            % INPUTS - For use when a directory path is given:
            %----------------------------------------------------------------
            %
            %   d: Directory path name as string.
            %
            %
            %   Optional parameter/value pair arguments:
            %
            %           'include_files': Optional cell array of strings giving file name prefixes
            %                           of files in directory 'd' to process. Files will be processed
            %                           in the order they are given within this cell array. *NOTE*: If
            %                           this argument is not given, *all* '_WL.mat' files in directory
            %                           'd' will be processed.
            %
            %           'ignore_files': Optional cell array of strings giving file name prefixes
            %                           (i.e., file names without the '_WL.mat' suffix/extension) to ignore.
            %                           Trumps 'include_files' argument.
            %
            %   Other parameter/value pairs from above are not valid because in this usage WhiskerTrialLite
            %   objects have previously been created and are here simply being loaded from disk.
            %
            %
            % DESCRIPTION:
            %
            % This constructor method takes either a WhiskerSignalTrialArray as input or a list of file
            % name prefixes of .mat files that contain individual WhiskerTrialLite objects.
            %
            %
            %
            if nargin==0
                return
            end
            
            p = inputParser;
            p.addOptional('w_or_d', '', @(x) ischar(x) || isa(x,'Whisker.WhiskerSignalTrialArray'));
            
            % Arguments for use when w_or_d is a WhiskerSignalTrialArray:
            p.addParamValue('polyRoiInPix', NaN);
            p.addParamValue('r_in_mm', 1, @(x) isnumeric(x) && numel(x)==1);
            p.addParamValue('calc_forces', false, @islogical);
            % If calc_forces is false the following are ignored:
            p.addParamValue('whisker_radius_at_base', 33.5, @isnumeric);
            p.addParamValue('whisker_length', 16, @isnumeric);
            p.addParamValue('youngs_modulus', 5e9, @isnumeric);
            p.addParamValue('baseline_time_end', 0.1, @isnumeric);
            
            % Arguments for use when w_or_d is a directory path:
            p.addParamValue('include_files', {}, @(x) all(cellfun(@ischar,x)));
            p.addParamValue('ignore_files', {}, @(x) all(cellfun(@ischar,x)));
            
            p.parse(w_or_d,varargin{:});
            
            if nargin==0
                return
            end
            
            disp 'List of all arguments:'
            disp(p.Results)
            
            if strcmp(p.Results.w_or_d,'')
                error('First argument not valid or not given.')
            end
            
            if isa(p.Results.w_or_d,'Whisker.WhiskerSignalTrialArray')
                doConstructor_WhiskerSignalTrialArray;
            else % directory path given as first argument
                doConstructor_directory;
            end
            
            function doConstructor_WhiskerSignalTrialArray % subfunction
                w = p.Results.w_or_d;
                obj.mouseName = w.mouseName;
                obj.sessionName = w.sessionName;
                
                ntrials = length(w);
                
                obj.trials = cell(1, ntrials);
                for k=1:ntrials
                    disp(['Processing trial ' int2str(k) ' of ' int2str(ntrials)])
                    obj.trials{k} = Whisker.WhiskerTrialLite(w.trials{k},varargin{:});
%                     obj.trials{k} = Whisker.WhiskerTrialLite(w.trials{k},'r_in_mm',...
%                         p.Results.r_in_mm,'calc_forces',p.Results.calc_forces,...
%                         'whisker_radius_at_base',p.Results.whisker_radius_at_base,...
%                         'whisker_length',p.Results.whisker_length,'youngs_modulus',p.Results.youngs_modulus,...
%                         'baseline_time_end',p.Results.baseline_time_end);
                end
            end
            
            function doConstructor_directory % subfunction
                d = p.Results.w_or_d;
                
                if ~strcmp(d(end), filesep)
                    d = [d filesep];
                end
                
                currentDir = pwd;
                cd(d)
                
                fnall = arrayfun(@(x) x.name(1:(end-7)), dir([d '*_WL.mat']),'UniformOutput',false);
                
                if ~isempty(p.Results.include_files) % Make sure files are found. If not, ignored.
                    ind = ismember(p.Results.include_files,fnall);
                    fnall = p.Results.include_files(ind);
                    if sum(ind) ~= numel(ind)
                        disp('The following files in ''include_files'' were not found in directory ''d'' and will be skipped:')
                        disp(p.Results.include_files(~ind))
                    end
                end
                
                if ~isempty(p.Results.ignore_files)
                    ind = ~ismember(fnall,p.Results.ignore_files);
                    fnall = fnall(ind);
                end
                
                inBoth = intersect(p.Results.include_files,p.Results.ignore_files);
                if ~isempty(inBoth)
                    disp('The following files were given in BOTH ''include_files'' and ''ignore files'' and will be ignored:')
                    disp(inBoth)
                end
                
                nfiles = length(fnall);
                
                obj.trials = cell(1, nfiles);
                if ~isempty(fnall)
                    for k=1:nfiles
                        fn = fnall{k};
                        disp(['Loading ''_WL.mat'' file '  fn ', ' int2str(k) ' of ' int2str(nfiles)])
                        
                        load([fn '_WL.mat'],'wl');
                        if ~isa(wl,'Whisker.WhiskerTrialLite')
                            error(['File ' fn '_WL.mat did not contain a WhiskerTrialLite object.'])
                        end
                        obj.trials{k} = wl;
                    end
                    
                    obj.mouseName = obj.trials{1}.mouseName; % Assume all files are from same mouse.
                    obj.sessionName = obj.trials{1}.sessionName; % Assume all files are from same session.
                end
                
                cd(currentDir)
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

