classdef WhiskerFinal_2pad < handle
    %
    %
    properties
        % basic information, from WL
        trialNum = [];
        trialType = '';
        framePeriodInSec = []; 
        pxPerMm = [];
        mouseName = '';
        sessionName = '';
        trackerFileName = '';
        
        poleAngle = [];
        polePosition = [];
        poleDistance = [] ;
        
        nof = [];
        time = []; % starting from 0, for the entire nof
        
        poleUpFrames = []; % First timepoint is 1, not 0. 2017/04/13 JK
        poleMovingFrames = [];
        
        % from WL (touch)
        protractionTFchunks = {}; % in frames, starting with 1
        protractionTFchunksByWhisking = {};
        retractionTFchunks = {};
        retractionTFchunksByWhisking = {};
        protractionTouchDuration = [];
        protractionTouchDurationByWhisking = [];
        retractionTouchDuration = [];
        retractionTouchDurationByWhisking = [];
        protractionSlide = {};
        protractionSlideByWhisking = {};
        retractionSlide = {};
        retractionSlideByWhisking = {};
        
        prothresholdMethod = []; % 1 if > mean + std, 2 if max > 10, 3 if kappa std > threshold, 0 otherwise
        rethresholdMethod = []; % 1 if > mean + std, 2 if max > 10, 3 if kappa std > threshold, 0 otherwise
        
        % from 3D (kinematics), recalculated to fill NaNs for not-tracked frames
        theta = []; % horizontal angle (top angle)
        phi = []; % elevation angle (front angle)
        kappaH = []; % horizaontal kappa
        kappaV = []; % vertical kappa
%         zeta = []; % roll angle, calculated by tangent line from the mask
        
        rInMm = [];
        arcLength = [];
    end
    
    properties (Dependent = true)        
        

    end
    
        
    methods (Access = public)
        function obj = WhiskerFinal_2pad(wl,w3)
            if nargin==0
                return
            end
            
            % basic information, from WL
            obj.trialNum = wl.trialNum;
            obj.trialType = wl.trialType;
            obj.framePeriodInSec = wl.framePeriodInSec;
            obj.pxPerMm = wl.pxPerMm;
            obj.mouseName = wl.mouseName;
            obj.sessionName = wl.sessionName;
            obj.trackerFileName = wl.trackerFileName;

            obj.poleAngle = wl.servoAngle;
            obj.polePosition = wl.apPosition;
            obj.poleDistance = wl.radialDistance;
            
            obj.nof = wl.nof;
            obj.time = (0:wl.nof-1) .* obj.framePeriodInSec;
            
            obj.poleUpFrames = wl.poleUpFrames;
            obj.poleMovingFrames = wl.poleMovingFrames;            
            
            % from WL (touch)
            obj.protractionTFchunks = wl.protractionTFchunks;
            obj.retractionTFchunks = wl.retractionTFchunks;
            if ~isempty(obj.protractionTFchunks)
                obj.protractionTouchDuration = cellfun(@(x) obj.time(x(end)) - obj.time(x(1)) + obj.framePeriodInSec, obj.protractionTFchunks);
            end
            if ~isempty(obj.retractionTFchunks)
                obj.retractionTouchDuration = cellfun(@(x) obj.time(x(end)) - obj.time(x(1)) + obj.framePeriodInSec, obj.retractionTFchunks);
            end
            
            if ~isempty(obj.protractionTFchunksByWhisking)
                obj.protractionTouchDuration = cellfun(@(x) obj.time(x(end)) - obj.time(x(1)) + obj.framePeriodInSec, obj.protractionTFchunksByWhisking);
            end
            if ~isempty(obj.retractionTFchunksByWhisking)
                obj.retractionTouchDuration = cellfun(@(x) obj.time(x(end)) - obj.time(x(1)) + obj.framePeriodInSec, obj.retractionTFchunksByWhisking);
            end
            
            obj.protractionSlide = wl.protractionSlide;
            obj.retractionSlide = wl.retractionSlide;
            obj.protractionSlideByWhisking = wl.protractionSlideByWhisking;
            obj.retractionSlideByWhisking = wl.retractionSlideByWhisking;
            
            obj.prothresholdMethod = wl.prothresholdMethod;
            obj.rethresholdMethod = wl.rethresholdMethod;
            
            % from 3D (kinematics), recalculated to fill NaNs for not-tracked frames
            obj.theta = nan(obj.nof,1);
            obj.phi = nan(obj.nof,1);
            obj.kappaH = nan(obj.nof,1);
            obj.kappaV = nan(obj.nof,1);
            [~,timeMatchingInds] = ismember(round(w3.time,4), round(obj.time,4));
            obj.theta(timeMatchingInds) = w3.theta;
            obj.phi(timeMatchingInds) = w3.phi;
            obj.kappaH(timeMatchingInds) = w3.kappaH;
            obj.kappaV(timeMatchingInds) = w3.kappaV;
            
            obj.rInMm = w3.rInMm;
            obj.arcLength = w3.lengthAlongWhisker;
            
            obj.protractionTFchunksByWhisking = wl.protractionTFchunksByWhisking;
            obj.retractionTFchunksByWhisking = wl.retractionTFchunksByWhisking;
        end
        
        
    end
    
    methods % Dependent property methods; cannot have attributes.

        
    end
       
end
