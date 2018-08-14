classdef WhiskerTrialLite_2pad < handle
    %
    % Modified from
    % 
    %   WhiskerTrialLite < handle
    %
    % and also combining
    %   WhiskerTrialLietI < WhiskerTrialLite
    %
    % Contains just final timeseries measurements from whiskers (kappa, theta, etc)
    % and meta-information in order to match with other data (ephys, imaging, etc).
    %
    %
    % DHO, 3/10.
    %
    properties
        trialNum = [];
        trialType = '';
        whiskerNames = {};  % whiskerNames and trajectoryIDs must be of same length
        trajectoryIDs = []; % with matching elements.
        framePeriodInSec = 1/311; % 311 Hz
        pxPerMm = []; %  Inherited from WhiskerSignalTrial.
        mouseName = '';
        sessionName = '';
        trackerFileName = '';
        useFlag = 1;
        
        rInMm = [];
        frontRInMm = [];
        
        time = {};
        deltaKappa = {};  
        thetaAtBase = {}; 
        follicleCoordsX = {}; 
        follicleCoordsY = {}; 
        
        % The following require pole tracking:
        thetaAtContact = {};
        distanceToPole = [];
        Fnorm = {};
        Faxial = {};
        Flateral = {};
        M0 = {}; % Moment at the follicle.
        meanKappa = {};
        distanceToPoleCenter = []; % dummy variable for now. It is populated during force calculation 2018/07/26 JK
        
        barPos = []; %  Inherited from WhiskerSignalTrial. [frameNum XPosition YPosition]
        barPosOffset = []; % Inherited from WhiskerSignalTrial. [x y], either 1X2 or nframesX2
        barRadius = []; % Inherited from WhiskerSignalTrial.  In pixels. Must be radius of bar tracked by the bar tracker.
        
        M0I = {};
        contactInds = {};

        servoAngle = []; % Inherited from WhiskerSignalTrial.
        apPosition = []; % Inherited from WhiskerSignalTrial.
        radialDistance = [] ; % Inherited from WhiskerSignalTrial.
        whiskerEdgeCoord = []; % Inherited from WhiskerSignalTrial.        
        nof = []; % Number of Frames. Inherited from WhiskerSignalTrial.
        poleUpFrames = []; % Inherited from WhiskerSignalTrial. first timepoint is 1, not 0. 2017/04/13 JK
        poleMovingFrames = [];
        mirrorAngle = 0;
        
        touchHP = []; % in 2d
        protractionHP = []; % in 2d
        retractionHP = []; % in 2d
        touchBoundaryThicknessInPix = [];
        touchBoundaryBufferInPix = [];
        touchPsi1 = [];
        touchPsi2 = [];
        stdHistogramThreshold = [];
        whiskingAmpThreshold = [];
        distanceHistogramBin = [];
        maxPointsNearHyperplane = [];
        touchKappaSTDthreshold = [];
        
        allTouchFrames = [];
        retractionTouchFrames = [];
        retractionTFchunks = {};
        protractionTouchFrames = [];
        protractionTFchunks = {};
        protractionThreshold = [];
        retractionThreshold = [];
        obviousNoTouchFrames = []; % just for check
        fwkappamean = []; % just for check
        fwkappastd = []; % just for check
        protractionDistance = []; % just for check
        retractionDistance = []; % just for check
        protractionFrames = [];
        retractionFrames = [];
        
        
        retractionTouchFramesPre = [];
        retractionTFchunksPre = {};
        protractionTouchFramesPre = [];
        protractionTFchunksPre = {};
        
        intersect_2d = []; % for test, maybe just for now 2018/07/31 JK
        prothresholdMethod = 0; % 1 if > mean + std, 2 if max > 10, 3 if kappa std > threshold, 0 otherwise
        rethresholdMethod = 0; % 1 if > mean + std, 2 if max > 10, 3 if kappa std > threshold, 0 otherwise        
    end
    
    properties (Dependent = true)
        theta % For compatability, make this 'alias' to refer to thetaAtBase
        kappa % For compatability, make this 'alias' to refer to deltaKappa
        M0Combined
    end
    
        
    methods (Access = public)
        function obj = WhiskerTrialLite_2pad(ws,varargin)
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
            p.addRequired('ws', @(x) isa(x,'Whisker.WhiskerSignalTrial_2pad'));
            p.addParameter('calc_forces', false, @islogical);
       
            p.addParameter('whisker_radius_at_base', 33.5, @isnumeric);
            p.addParameter('whisker_length', 16, @isnumeric);
            p.addParameter('youngs_modulus', 5e9, @isnumeric);
            p.addParameter('baseline_time_or_kappa_value', [0 0.1], @isnumeric);
            p.addParameter('proximity_threshold', -1, @isnumeric);
            
            p.addParameter('behavior',[], @(x) isa(x,'Solo.BehavTrial2padArray'));
            
            p.addParameter('trial_type',{}, @ischar);
            p.addParameter('kappaTouchThreshold',[],@(x) isnumeric(x) && numel(x)==2); % 2 values for top-view and front-view kappa
            p.addParameter('durationThreshold',[], @isnumeric);
            p.addParameter('mirrorAngle', 0, @isnumeric); % averaged from all the trials in the session
            p.addParameter('rInMm',{}, @isnumeric);
            p.addParameter('touchHP',[], @(x) isnumeric(x) && size(x,1) == 3);
            p.addParameter('touchHPpeaks',[], @(x) isnumeric(x) && numel(x) ==2);
            p.addParameter('touchPsi1', [], @isnumeric);
            p.addParameter('touchPsi2', [], @isnumeric);
            p.addParameter('touchKappaSTDthreshold', 2, @isnumeric);
            p.addParameter('fwkappamean', 0, @isnumeric);
            p.addParameter('fwkappastd', 0, @isnumeric);
            p.addParameter('stdHistogramThreshold', 3, @isnumeric);
            p.addParameter('whiskingAmpThreshold', 2.5, @isnumeric);
            p.addParameter('distanceHistogramBin', 0.2, @isnumeric);
            p.addParameter('touchBoundaryThickness', 0.25, @isnumeric);
            p.addParameter('touchBoundaryBuffer', 0.1, @isnumeric);
            p.addParameter('maxPointsNearHyperplane', 10, @isnumeric);
            
            p.parse(ws,varargin{:});
            
            obj.trialNum = ws.trialNum;
            obj.trialType = ws.trialType;
            obj.whiskerNames = ws.whiskerNames;
            obj.trajectoryIDs = ws.trajectoryIDs;
            obj.framePeriodInSec = ws.framePeriodInSec;
            obj.pxPerMm = ws.pxPerMm;
            obj.mouseName = ws.mouseName;
            obj.sessionName = ws.sessionName;
            obj.trackerFileName = ws.trackerFileName;
            obj.follicleCoordsX = ws.follicleCoordsX;
            obj.follicleCoordsY = ws.follicleCoordsY;  
            obj.useFlag = ws.useFlag;
            obj.rInMm = p.Results.rInMm;
            
            ntraj = length(obj.trajectoryIDs);
                        
            obj.time = cell(1,ntraj);
            obj.deltaKappa = cell(1,ntraj);  
            obj.thetaAtBase = cell(1,ntraj);     
            obj.thetaAtContact = cell(1,ntraj);            
            obj.Fnorm = cell(1,ntraj);
            obj.Faxial = cell(1,ntraj);
            obj.Flateral = cell(1,ntraj);
            obj.M0 = cell(1,ntraj); % Moment at the follicle.
            obj.meanKappa = cell(1,ntraj); % Mean kappa over the ROI.
            obj.mirrorAngle = p.Results.mirrorAngle;
          
            obj.barPos = ws.barPos; %  Inherited from WhiskerSignalTrial. [frameNum XPosition YPosition]
            obj.barPosOffset = ws.barPosOffset; % Inherited from WhiskerSignalTrial. [x y], either 1X2 or nframesX2
            obj.barRadius = ws.barRadius; % Inherited from WhiskerSignalTrial.  In pixels. Must be radius of bar tracked by the bar tracker.
            
            obj.servoAngle = ws.angle;
            obj.apPosition = ws.apPosition;
            obj.radialDistance = ws.radialDistance;
            obj.trialType = ws.trialType;

            obj.touchPsi1 = p.Results.touchPsi1;
            obj.touchPsi2 = p.Results.touchPsi2;
            obj.touchKappaSTDthreshold = p.Results.touchKappaSTDthreshold;
            obj.fwkappamean = p.Results.fwkappamean;
            obj.fwkappastd = p.Results.fwkappastd;
            obj.touchBoundaryThicknessInPix = obj.pxPerMm * p.Results.touchBoundaryThickness; % in pixels for detecting touch frames (later to be used with kappa change)
            obj.touchBoundaryBufferInPix = obj.pxPerMm * p.Results.touchBoundaryBuffer; % in pixels for touch distance threshold
            obj.stdHistogramThreshold = p.Results.stdHistogramThreshold;
            obj.whiskingAmpThreshold = p.Results.whiskingAmpThreshold;
            obj.distanceHistogramBin = p.Results.distanceHistogramBin; % in pixels
            obj.maxPointsNearHyperplane = p.Results.maxPointsNearHyperplane;
            
            obj.whiskerEdgeCoord = ws.whiskerEdgeCoord;
            obj.poleUpFrames = ws.poleUpFrames;
            obj.poleMovingFrames = ws.poleMovingFrames;
            obj.distanceToPole = ws.dist2pole;

            obj.nof = ws.nof;
            obj.frontRInMm = ws.get_frontRInMm(obj.rInMm);

%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             obj.time = (0:obj.nof-1)*obj.framePeriodInSec;
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change this to using 'timestamp'
%             
            for k=1:ntraj
                tid = obj.trajectoryIDs(k);
                
                obj.time{k} = ws.get_time(tid);
                
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
                   [obj.M0{k},obj.Faxial{k},~,obj.deltaKappa{k},obj.Fnorm{k},...
                       obj.thetaAtBase{k},obj.thetaAtContact{k},obj.distanceToPoleCenter{k}, obj.meanKappa{k}, obj.Flateral{k}] = ...
                       ws.calc_M0_Faxial(tid,p.Results.r_in_mm,p.Results.whisker_radius_at_base, p.Results.whisker_length,...
                       p.Results.youngs_modulus,p.Results.baseline_time_or_kappa_value,p.Results.proximity_threshold);                    
                else
                    % Should consolidate into single function to optimize the following: 
                    if k == 1
                        [dk,~,~,~] = ws.get_kappa_at_roi_point(tid,p.Results.rInMm);
                        [tab,~] = ws.get_theta_at_base(tid);                        
                        inds = round(obj.time{k}/obj.framePeriodInSec) + 1;
                        obj.deltaKappa{k} = nan(obj.nof,1);
                        obj.deltaKappa{k}(inds) = dk;
                        obj.thetaAtBase{k} = nan(obj.nof,1);
                        obj.thetaAtBase{k}(inds) = tab + p.Results.mirrorAngle;
                    else % k == 2
                        [dk,~,~,~] = ws.get_kappa_at_roi_point(tid,obj.frontRInMm);
                        [tab,~] = ws.get_theta_at_base(tid);
                        inds = round(obj.time{k}/obj.framePeriodInSec) + 1;
                        obj.deltaKappa{k} = nan(obj.nof,1);
                        obj.deltaKappa{k}(inds) = dk;
                        obj.thetaAtBase{k} = nan(obj.nof,1);
                        obj.thetaAtBase{k}(inds) = tab;
                    end                    
                end
            end

            %% Temporary remedy 2017/05/30
%             thp1 = obj.thPolygon(1:end/2,:);
%             thp1 = [2*thp1(1,:) - thp1(end,:); thp1; 2*thp1(end,:) - thp1(1,:)];
%             thp2 = obj.thPolygon(end/2+1:end,:);
%             thp2 = [2*thp2(1,:) - thp2(end,:); thp2; 2*thp2(end,:) - thp2(1,:)];
%             obj.thPolygon = [thp1; thp2];            
            %% Finding touch frames
            if (contains(ws.sessionName, 'S') || contains(ws.sessionName, 'pre')) && ~strcmp(ws.trialType, 'oo')
                
                % first, sort out any obvious non-touch frames, where whisker does not overlap at all with top-view pole, except for 90 degrees 
                % (assigning these frames into touch frames happens more frequently in lower degrees (45~75)
                % Only for pole-up frames. Touch detection during pole-moving frames is highly unreliable
                
                obj.obviousNoTouchFrames = zeros(obj.nof,1,'logical');
                if obj.servoAngle ~= 90
                    height = size(ws.binvavg,1);
                    width = size(ws.binvavg,2);
                    topViewPole = zeros([height, width],'logical');
                    targetH = 0.6; targetW = 0.5; targetArea = zeros(size(ws.binvavg), 'logical'); targetArea(1:round(height*targetH), round((1-targetW)*width):end) = deal(1);
                    targetInd = find(targetArea);
                    bw = bwconncomp(ws.binvavg);
                    bwInd = find(cellfun(@(x) length(intersect(x,targetInd)),bw.PixelIdxList));
                    if isempty(bwInd)
                        error(['no top-view pole detected at mouse ', obj.mouseName, ' session ', obj.sessionName, ' trial # ', num2str(obj.trialNum)]);
                    elseif length(bwInd) > 1
                        pixLength = cellfun(@(x) length(x), bw.PixelIdxList);
                        [~,tempInd] = max(pixLength(bwInd));
                        bwInd = bwInd(tempInd);
                    end
                    topViewPole(bw.PixelIdxList{bwInd}) = deal(1);
                    
                    for tempi = 1 : length(obj.poleUpFrames)
                        trackerFrameTop = find(ws.trackerFrames{1} == obj.poleUpFrames(tempi),1);
                        if ~isempty(trackerFrameTop)
                            x = polyval(ws.polyFits{1}{1}(trackerFrameTop,:),linspace(0,1.3)); % stretch the whisker fitting outwardly 30 % for cases where whisker tracing is cut off because of the pole
                            y = polyval(ws.polyFits{1}{2}(trackerFrameTop,:),linspace(0,1.3));
                            xyi = find(y >= 1 & y <= height & x >= 1 & x <= width); % take only those within image dimension
                            x = x(xyi); y = y(xyi);
                            whiskerBW = zeros([height, width],'logical');
                            whiskerBW(sub2ind(size(whiskerBW),round(y),round(x))) = deal(1);
                            whiskerPoleIntersection = intersect(find(whiskerBW), find(topViewPole));
                            if isempty(whiskerPoleIntersection)
                                obj.obviousNoTouchFrames(obj.poleUpFrames(tempi)) = 1;
                            end
                        end
                    end
                end
                obj.obviousNoTouchFrames = find(obj.obviousNoTouchFrames);
                
                % Now find touch frames
                noNaNInd = intersect(find(~isnan(sum(ws.whiskerEdgeCoord,2))), ws.poleUpFrames);
                intersect_3d_total = [ws.whiskerEdgeCoord(noNaNInd,1), ws.whiskerEdgeCoord(noNaNInd,2), ws.apPosition(noNaNInd)];

                intersect_4d = [intersect_3d_total, ones(size(intersect_3d_total,1),1)]';
                if obj.touchPsi1 <= 90
                    projMat = viewmtx(obj.touchPsi1, 90 - obj.touchPsi2);
                else
                    projMat = viewmtx(obj.touchPsi1, -90 + obj.touchPsi2);
                end
                intersect_2d = projMat*intersect_4d;
                if size(intersect_2d,1) < size(intersect_2d,2)
                    intersect_2d = intersect_2d';
                end
                intersect_2d = intersect_2d(:,1:2);
                obj.intersect_2d = intersect_2d;
                
                touchHP4d = [p.Results.touchHP; ones(1,size(p.Results.touchHP,2))];
                touchHP = projMat*touchHP4d;
                obj.touchHP = unique(round(touchHP(1:2,:)'*100)/100,'rows');
                proHP4d = [p.Results.touchHP(1,:) + p.Results.touchHPpeaks(2); p.Results.touchHP(2:3,:); ones(1,size(p.Results.touchHP,2))];
                proHP = projMat*proHP4d;
                obj.protractionHP = unique(round(proHP(1:2,:)'*100)/100, 'rows');
                retHP4d = [p.Results.touchHP(1,:) + p.Results.touchHPpeaks(1); p.Results.touchHP(2:3,:); ones(1,size(p.Results.touchHP,2))];
                retHP = projMat*retHP4d;
                obj.retractionHP = unique(round(retHP(1:2,:)'*100)/100, 'rows');
                
                protractionDistance = obj.distance_and_side_from_line(intersect_2d, obj.protractionHP);
                tempProtFrames = find(protractionDistance > 0);
                retractionDistance = obj.distance_and_side_from_line(intersect_2d, obj.retractionHP);
                tempRetFrames = find(retractionDistance < 0);
                inFrames = setdiff(1:length(noNaNInd), union(tempProtFrames, tempRetFrames));
                inChunks = obj.get_chunks(inFrames);
                if ~isempty(inChunks)
%                     if inChunks{1}(1) == 1
%                         if abs(obj.distance_and_side_from_line(intersect_2d(inChunks{1}(end),:), obj.protractionHP)) ...
%                                 <= abs(obj.distance_and_side_from_line(intersect_2d(inChunks{1}(end),:), obj.retractionHP)) % (end) of whiskerEdgeCoord should be more stable than (1)
%                             tempProtFrames = [tempProtFrames; inChunks{1}];
%                         else
%                             tempRetFrames = [tempRetFrames; inChunks{1}];
%                         end
%                     else
%                         if ismember(inChunks{1}(1)-1, tempProtFrames)
%                             tempProtFrames = [tempProtFrames; inChunks{1}];
%                         elseif ismember(inChunks{1}(1)-1, tempRetFrames)
%                             tempRetFrames = [tempRetFrames; inChunks{1}];
%                         else
%                             error(['Cannot assign to a side of the pole in trial #', num2str(obj.trialNum)])
%                         end
%                     end
%                     for i = 2 : length(inChunks)
%                         if ismember(inChunks{i}(1)-1, tempProtFrames)
%                             tempProtFrames = [tempProtFrames; inChunks{i}];
%                         elseif ismember(inChunks{i}(1)-1, tempRetFrames)
%                             tempRetFrames = [tempRetFrames; inChunks{i}];
%                         else
%                             error(['Cannot assign to a side of the pole in trial #', num2str(obj.trialNum)])
%                         end
%                     end
                    prodist = cellfun(@(x) mean(abs(protractionDistance(x))), inChunks);
                    retdist = cellfun(@(x) mean(abs(retractionDistance(x))), inChunks);
                    for i = 1 : length(inChunks)
                        if ismember(inChunks{i}(1)-1, tempProtFrames)
                            tempProtFrames = [tempProtFrames; inChunks{i}];
                        elseif ismember(inChunks{i}(1)-1, tempRetFrames)
                            tempRetFrames = [tempRetFrames; inChunks{i}];
                        elseif prodist(i) > retdist(i)
                            tempRetFrames = [tempRetFrames; inChunks{i}];
                        else
                            tempProtFrames = [tempProtFrames; inChunks{i}];
                        end
                    end
                end
                tempProtFrames = sort(tempProtFrames);
                tempRetFrames = sort(tempRetFrames);
                obj.protractionFrames = noNaNInd(tempProtFrames);
                protractionDistance = protractionDistance(tempProtFrames);
                obj.retractionFrames = noNaNInd(tempRetFrames);
                retractionDistance = retractionDistance(tempRetFrames);

                obj.protractionDistance = protractionDistance;                
                obj.retractionDistance = retractionDistance;
                
%                 % (1) protraction touch
%                 if ~isempty( protractionDistance <= obj.touchBoundaryThickness) % meaning when there is possible touch frames
%                     threshold = obj.touchBoundaryBuffer; % minimum threshold value (correct when the hyperplane is perfect
%                     % (1) - 1. Number of points at the distance with the densest population has more than 1/3 of all possible
%                     % points between 0 and "touch boundary thickness" pixels from the boundary, when divided into 0.25 pixels
%                     closeFrames = find(protractionDistance <= obj.touchBoundaryThickness);
%                     tempDist = round(protractionDistance(closeFrames)*4);
%                     modeVal = mode(tempDist);
%                     if length(find(tempDist == modeVal)) > length(closeFrames)/3
%                         threshold = modeVal/4 + obj.touchBoundaryBuffer;
%                     else
%                         % (1) - 2. Closest points have kappa larger or shorter than the mean by a certain threshold, 
%                         % defined by std of free whisking frames
%                         protractionDistance( (protractionDistance < 0) ) = deal(0);
%                         minDistInds = find( round(protractionDistance) == min(round(protractionDistance)) );                    
%                         candid = find(obj.deltaKappa{1}(protractionFrames(minDistInds)) < obj.fwkappamean- obj.touchKappaSTDthreshold * obj.fwkappastd & ...
%                             obj.deltaKappa{1}(protractionFrames(minDistInds)) > obj.fwkappamean + obj.touchKappaSTDthreshold * obj.fwkappastd); % meaning there IS (ARE) protraction touch frames. There can be multiple frames with min value
%                         if ~isempty(candid)
%                             threshold = mode(round(protractionDistance(minDistInds(candid)))) + obj.touchBoundaryBuffer;                            
%                         end
%                     end
%                     obj.protractionTouchFrames = protractionFrames(protractionDistance <= threshold);
%                 end                
%                 obj.protractionTouchFrames = setdiff(obj.protractionTouchFrames, obj.obviousNoTouchFrames);
%                 
%                 % (2) retraction touch
%                 if ~isempty( retractionDistance >= -obj.touchBoundaryThickness) % meaning when there is possible touch frames
%                     threshold = -obj.touchBoundaryBuffer; % minimum threshold value (correct when the hyperplane is perfect
%                     % (2) - 1. Number of points at the distance with the densest population has more than 1/3 of all possible
%                     % points between 0 and "touch boundary thickness" pixels from the boundary, when divided into 0.25 pixels
%                     closeFrames = find(retractionDistance >= -obj.touchBoundaryThickness);
%                     tempDist = round(retractionDistance(closeFrames)*4);
%                     modeVal = mode(tempDist);
%                     if length(find(tempDist == modeVal)) > length(closeFrames)/3
%                         threshold = modeVal/4 - obj.touchBoundaryBuffer;
%                     else
%                         % (2) - 2. Closest points have kappa larger or shorter than the mean by a certain threshold, 
%                         % defined by std of free whisking frames
%                         retractionDistance( (retractionDistance > 0) ) = deal(0);
%                         maxDistInds = find( round(retractionDistance) == max(round(retractionDistance)) ); % max value is negative or 0                    
%                         candid = find(obj.deltaKappa{1}(retractionFrames(maxDistInds)) > obj.fwkappamean + obj.touchKappaSTDthreshold * obj.fwkappastd & ...
%                             obj.deltaKappa{1}(retractionFrames(maxDistInds)) < obj.fwkappamean - obj.touchKappaSTDthreshold * obj.fwkappastd); % meaning there IS (ARE) protraction touch frames. There can be multiple frames with min value
%                         if ~isempty(candid)
%                             threshold = mode(round(retractionDistance(maxDistInds(candid)))) - obj.touchBoundaryBuffer;                            
%                         end
%                     end
%                     obj.retractionTouchFrames = retractionFrames(retractionDistance >= threshold);
%                 end
%                 obj.retractionTouchFrames = setdiff(obj.retractionTouchFrames, obj.obviousNoTouchFrames);
                
                % a new way of calculating touch distance: within 5 pixels from the touch hp, find the peak distance with consecutive frames, 
                % and if the occurance at this distance is significantly higher than all the other distances, set this distance as the touch distance. 
                % (1) protraction touch
                if ~isempty( protractionDistance <= obj.touchBoundaryThicknessInPix) % meaning when there is possible touch frames                    
                    % (1) - 1. Number of points at the distance with the densest population has more than "n" X std of all possible
                    % points between 0 and "touch boundary thickness" pixels from the boundary, when divided into "bin" pixels
                    % In case of protraction, consider only during whisking (amplitude > "amplitude threshold")
                    closeFrames = find(protractionDistance <= obj.touchBoundaryThicknessInPix);                    
                    closeAdjFrames = closeFrames(find([0;diff(closeFrames)] == 1));
                    if ~isempty(closeAdjFrames)
                        [~, amplitude, ~, ~, ~, ~, phase, ~] = jkWhiskerDecomposition(obj.thetaAtBase{1});
                        whiskingInds = [1;find([0;diff(phase)]<0); length(phase)+1];
                        if ~isempty(whiskingInds)
                            whiskingFrames = [];
                            for i = 1 : length(whiskingInds)-1
                                if max(amplitude(whiskingInds(i):whiskingInds(i+1)-1)) > obj.whiskingAmpThreshold
                                    whiskingFrames = [whiskingFrames, whiskingInds(i):whiskingInds(i+1)-1];
                                end
                            end
                            closeProtractionDist = protractionDistance(intersect(closeAdjFrames, whiskingFrames));
                            if length(closeProtractionDist) > 1 && max(closeProtractionDist) - min(closeProtractionDist) > obj.distanceHistogramBin
                                closeProtractionDist = round(closeProtractionDist / obj.distanceHistogramBin) * obj.distanceHistogramBin;                                
                                [N, edges] = histcounts(closeProtractionDist, [min(closeProtractionDist) : obj.distanceHistogramBin : max(closeProtractionDist)]);
                                if max(N) > mean(N) + std(N) * obj.stdHistogramThreshold
                                    [~, maxind] = max(N);
                                    obj.prothresholdMethod = 1;
                                    obj.protractionThreshold = ( edges(maxind) + edges(maxind+1) ) / 2 + obj.touchBoundaryBufferInPix;
                                elseif max(N) > obj.maxPointsNearHyperplane
                                    [~, maxind] = max(N);
                                    obj.prothresholdMethod = 2;
                                    obj.protractionThreshold = ( edges(maxind) + edges(maxind+1) ) / 2 + obj.touchBoundaryBufferInPix;                                    
                                end
                            end
                        end
                        if isempty(obj.protractionThreshold) % in case where the above method could not be used.
                            % (1) - 2. Closest points have kappa larger or shorter than the mean by a certain threshold, 
                            % defined by std of free whisking frames
                            proDist = protractionDistance;
                            proDist( (proDist < 0) ) = deal(0);
                            proDist = round(proDist / obj.distanceHistogramBin) * obj.distanceHistogramBin;
                            minDistInds = find( proDist == min(proDist) );
                            candid = find(obj.deltaKappa{1}(obj.protractionFrames(minDistInds)) < obj.fwkappamean - obj.touchKappaSTDthreshold * obj.fwkappastd & ...
                                obj.deltaKappa{1}(obj.protractionFrames(minDistInds)) > obj.fwkappamean + obj.touchKappaSTDthreshold * obj.fwkappastd); % meaning there IS (ARE) protraction touch frames. There can be multiple frames with min value
                            if ~isempty(candid)
                                obj.prothresholdMethod = 3;
                                obj.protractionThreshold = mode(proDist(minDistInds(candid))) + obj.touchBoundaryBufferInPix;                            
                            end
                        end
                    end
                end
                if ~isempty(obj.protractionThreshold)
                    obj.protractionTouchFramesPre = obj.protractionFrames(obj.protractionDistance <= obj.protractionThreshold);
                    obj.protractionTouchFramesPre = setdiff(obj.protractionTouchFramesPre, obj.obviousNoTouchFrames);
                    obj.protractionTFchunksPre = obj.get_chunks(obj.protractionTouchFramesPre);
                end
                
                % (2) retraction touch
                if ~isempty( retractionDistance >= -obj.touchBoundaryThicknessInPix) % meaning when there is possible touch frames                    
                    % (2) - 1. Number of points at the distance with the densest population has more than "n" X std of all possible
                    % points between 0 and "touch boundary thickness" pixels from the boundary, when divided into "bin" pixels
                    closeFrames = find(retractionDistance >= -obj.touchBoundaryThicknessInPix);                    
                    closeAdjFrames = closeFrames(find([0;diff(closeFrames)] == 1));
                    if ~isempty(closeAdjFrames)
                        closeRetractionDist = retractionDistance(closeAdjFrames);
                        if length(closeRetractionDist) > 1 && max(closeRetractionDist) - min(closeRetractionDist) > obj.distanceHistogramBin 
                            closeRetractionDist = round(closeRetractionDist / obj.distanceHistogramBin) * obj.distanceHistogramBin;
                            [N, edges] = histcounts(closeRetractionDist, [min(closeRetractionDist) : obj.distanceHistogramBin : max(closeRetractionDist)]);
                            if max(N) > mean(N) + std(N) * obj.stdHistogramThreshold
                                obj.rethresholdMethod = 1;
                                [~, maxind] = max(N);
                                obj.retractionThreshold = ( edges(maxind) + edges(maxind+1) ) / 2 - obj.touchBoundaryBufferInPix;
                            elseif max(N) > obj.maxPointsNearHyperplane
                                obj.rethresholdMethod = 2;
                                [~, maxind] = max(N);
                                obj.retractionThreshold = ( edges(maxind) + edges(maxind+1) ) / 2 - obj.touchBoundaryBufferInPix;
                            end
                        end
                        if isempty(obj.retractionThreshold) % in case where the above method could not be used.
                            % (2) - 2. Closest points have kappa larger or shorter than the mean by a certain threshold, 
                            % defined by std of free whisking frames
                            retDist = retractionDistance;
                            retDist( (retDist > 0) ) = deal(0);
                            retDist = round(retDist / obj.distanceHistogramBin) * obj.distanceHistogramBin;
                            minDistInds = find( retDist == max(retDist) );
                            candid = find(obj.deltaKappa{1}(obj.retractionFrames(minDistInds)) < obj.fwkappamean - obj.touchKappaSTDthreshold * obj.fwkappastd & ...
                                obj.deltaKappa{1}(obj.retractionFrames(minDistInds)) > obj.fwkappamean + obj.touchKappaSTDthreshold * obj.fwkappastd); % meaning there IS (ARE) protraction touch frames. There can be multiple frames with min value
                            if ~isempty(candid)
                                obj.rethresholdMethod = 3;
                                obj.retractionThreshold = mode(retDist(minDistInds(candid))) - obj.touchBoundaryBufferInPix;                            
                            end
                        end
                    end
                end
                if ~isempty(obj.retractionThreshold)
                    obj.retractionTouchFramesPre = obj.retractionFrames(obj.retractionDistance >= obj.retractionThreshold);
                    obj.retractionTouchFramesPre = setdiff(obj.retractionTouchFramesPre, obj.obviousNoTouchFrames);
                    obj.retractionTFchunksPre = obj.get_chunks(obj.retractionTouchFramesPre);
                end
                
                % Combine both protraction and retraciton touches in pre-,
                % and re-chunk them to assign protraction and retraction
                % touches
                
                if ~isempty(obj.protractionTouchFramesPre) && ~isempty(obj.retractionTouchFramesPre)
                    protractionDistance = obj.distance_and_side_from_line(intersect_2d, obj.protractionHP);
                    tempProtFrames = find(protractionDistance > obj.protractionThreshold);
                    retractionDistance = obj.distance_and_side_from_line(intersect_2d, obj.retractionHP);
                    tempRetFrames = find(retractionDistance < obj.retractionThreshold);
                    
                    obj.allTouchFrames = union(obj.protractionTouchFramesPre, obj.retractionTouchFramesPre); % union results in sorted order
                    [~, allTouchFrames] = ismember(obj.allTouchFrames, noNaNInd); % change back to noNaNInd indexing for further processing
                    if ~isempty(find(allTouchFrames == 0,1)) % something's wrong
                        error(['allTouchFrames not included in noNaNInd in trial #', num2str(wl.trialNum)])
                    end
                    touchFramesChunks = obj.get_chunks(allTouchFrames);
                    prodist = cellfun(@(x) mean(abs(protractionDistance(x))), touchFramesChunks);
                    retdist = cellfun(@(x) mean(abs(retractionDistance(x))), touchFramesChunks);
                    proTF = [];
                    retTF = [];
                    for i = 1 : length(touchFramesChunks)
                        if ismember(touchFramesChunks{i}(1)-1, tempProtFrames)
                            proTF = [proTF; touchFramesChunks{i}];
                            tempProtFrames = [tempProtFrames; touchFramesChunks{i}];
                        elseif ismember(touchFramesChunks{i}(1)-1, tempRetFrames)
                            retTF = [retTF; touchFramesChunks{i}];
                            tempRetFrames = [tempRetFrames; touchFramesChunks{i}];
                        elseif prodist(i) > retdist(i)
                            retTF = [retTF; touchFramesChunks{i}];
                            tempRetFrames = [tempRetFrames; touchFramesChunks{i}];
                        else
                            proTF = [proTF; touchFramesChunks{i}];
                            tempProtFrames = [tempProtFrames; touchFramesChunks{i}];
                        end
                    end
                    obj.protractionFrames = noNaNInd(unique(tempProtFrames));
                    obj.retractionFrames = noNaNInd(unique(tempRetFrames));
                    obj.protractionDistance = protractionDistance(tempProtFrames);
                    obj.retractionDistance = retractionDistance(tempRetFrames);
                    if ~isempty(proTF)
                        obj.protractionTouchFrames = noNaNInd(sort(proTF));
                        obj.protractionTFchunks = obj.get_chunks(obj.protractionTouchFrames);
                    end
                    if ~isempty(retTF)
                        obj.retractionTouchFrames = noNaNInd(sort(retTF));
                        obj.retractionTFchunks = obj.get_chunks(obj.retractionTouchFrames);
                    end
                else
                    obj.protractionTouchFrames = obj.protractionTouchFramesPre;
                    obj.protractionTFchunks = obj.protractionTFchunksPre;
                    obj.retractionTouchFrames = obj.retractionTouchFramesPre;
                    obj.retractionTFchunks = obj.retractionTFchunksPre;
                end
            end
                

            %
            %
            % for debugging
%             figure, 
%             plot(intersect_2d(1, dist2HP >= obj.touchBoundaryThickness), intersect_2d(2, dist2HP >= obj.touchBoundaryThickness), 'k.'), hold on
%             plot(obj.thPolygon(:,1), obj.thPolygon(:,2), 'r-')
%             plot(intersect_2d(1, dist2HP < obj.touchBoundaryThickness ), intersect_2d(2, dist2HP < obj.touchBoundaryThickness ), 'b.')
            %
            %
            %           
        end
        
        function distance = distance_and_side_from_line(obj, points, line)
        % points right side of the line is assigned to be positive.
        % distanceNside returns distance from the line and + if the points are on the right side (protraction space), 0 on the line, - otherwise (retraction space). 
            [~, maxind] = max(line(:,1));
            v2 = line(maxind,:);
            [~, minind] = min(line(:,1));
            v1temp = line(minind,:);
            if v2(1) - v1temp(1) == 0
                error('Hyperplane 2d slope is 90 degrees')
            else
                hpSlope = (v2(2) - v1temp(2)) / (v2(1) - v1temp(1));
            end
            [~, minind] = min(points(:,1));
            v1 = [points(minind,1), v2(2) - (v2(1) - points(minind,1)) * hpSlope];

            v1_ = repmat([v1,0], [size(points,1),1]);
            v2_ = repmat([v2,0], [size(points,1),1]);
            points = [points, zeros(size(points,1),1)];
            a = v1_ - v2_;
            b = points - v2_;
%                 dist2touchHP = p_poly_dist(intersect_2d(1,:), intersect_2d(2,:), obj.touchHP(:,1)', obj.touchHP(:,2)');
            distance = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));
            if size(distance,2) > size(distance,1)
                distance = distance';
            end
            intersectSlopes = (points(:,2) - v1(2)) ./ (points(:,1) - v1(1));
            negind = find(intersectSlopes > hpSlope);
            distance(negind) = -distance(negind);
        end
        
        function chunks = get_chunks(obj, frames)
            if isempty(frames)
                chunks = {};
            else
                if size(frames,1) > size(frames,2)
                    frames = frames';
                end
                chunkPoints = [1, find(diff(frames)>1) + 1, length(frames)+1]; % +1 because of diff. first 1 for the first chunk. So this is actually start points of each chunk.
                chunks = cell(1,length(chunkPoints)-1); % 
                for i = 1 : length(chunks)
                    chunks{i} = [frames(chunkPoints(i) : chunkPoints(i+1)-1)]';
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
        
        function value = get.M0Combined(obj)
            if isempty(obj.barPos)
                value = [];
            else
               value = [1];
            end
        end
        
    end
       
end
