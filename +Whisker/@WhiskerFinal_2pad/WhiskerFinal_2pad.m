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
        retractionTFchunks = {};
        protractionTouchDuration = [];
        retractionTouchDuration = [];
        protractionSlide = {};
        retractionSlide = {};
        
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
        protractionTFchunksByWhisking
        protractionTouchDurationByWhisking
        protractionSlideByWhisking
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
                obj.protractionTouchDuration = cellfun(@(x) length(x), obj.protractionTFchunks);
            end
            if ~isempty(obj.retractionTFchunks)
                obj.retractionTouchDuration = cellfun(@(x) length(x), obj.retractionTFchunks);
            end
            obj.protractionSlide = wl.protractionSlide;
            obj.retractionSlide = wl.retractionSlide;
            
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
        end
        
        
    end
    
    methods % Dependent property methods; cannot have attributes.
        function value = get.protractionTFchunksByWhisking(obj)
            if isempty(obj.protractionTFchunks)
                value = {};
            else
                [onsetFrames, ~, ~, ~, peakFrames] = jkWhiskerOnsetNAmplitude(obj.theta, 2);
                wholeChunk = cell2mat(obj.protractionTFchunks');
                protractionFrames = cell(1,length(find(isfinite(peakFrames))));
                for i = 1 : length(find(isfinite(peakFrames)))
                    protractionFrames{i} = onsetFrames(i):peakFrames(i);
                end
                wholeProtractionFrames = cell2mat(protractionFrames);
                frames = intersect(wholeChunk, wholeProtractionFrames);
                if isempty(frames)
                    value = {};
                else
                    if size(frames,1) > size(frames,2)
                        frames = frames';
                    end
                    chunkPoints = [1, find(diff(frames)>1) + 1, length(frames)+1]; % +1 because of diff. first 1 for the first chunk. So this is actually start points of each chunk.
                    value = cell(1,length(chunkPoints)-1); % 
                    for i = 1 : length(value)
                        value{i} = [frames(chunkPoints(i) : chunkPoints(i+1)-1)]';
                    end
                end
            end
        end
        
        function value = get.protractionTouchDurationByWhisking(obj)
            if ~isempty(obj.protractionTFchunks)
                value = cellfun(@(x) length(x), obj.protractionTFchunksByWhisking);
            end
        end
        
        function value = get.protractionSlideByWhisking(obj)
            if ~isempty(obj.protractionTFchunksByWhisking)
                wholeChunks = cell2mat(obj.protractionTFchunks');
                protractionChunks = cell2mat(obj.protractionTFchunksByWhisking');
                chunkInds = [0, cumsum(cellfun(@(x) length(x), obj.protractionTFchunksByWhisking))];
                [~, matchingInd] = ismember(protractionChunks, wholeChunks);
                wholeSlide = cell2mat(obj.protractionSlide');                
                protractionSlide = wholeSlide(matchingInd);
                value = cell(1,length(obj.protractionTFchunksByWhisking));
                for i = 1 : length(value)
                    value{i} = protractionSlide(chunkInds(i)+1:chunkInds(i+1));
                end
            end
        end
        
    end
       
end
