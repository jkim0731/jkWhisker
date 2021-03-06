classdef WhiskerSignalTrial_2pad < Whisker.WhiskerSignalTrial
   
    % Infer whisker-pole contact frames (index and time) from angle &
    % radial distance tasks (2pad, by JK). The calculation is based on the
    % intersection points between whisker and the pole edge, from both the
    % front and top view. This is based on the assumption that whisker-pole
    % contact points are determined uniquely when combining both of the
    % views. 
    %
    % 2017/04/10 JK   
    
    % Added automatic detection of pole_available_frames 2018/03/06 JK
    
    properties
% uncomment this when I need to use show_nan_polyfits_with_mask        
        trackerData = {}; 
        whiskerPadOrigin = [];
        trackerFrames = {};
        whiskerPoleIntersection = {}; 
        whiskerEdgeCoord = [];
        imagePixelDimsXY = []; % [NumberOfXPixels NumberOfYPixels]
        poleAxesUp = cell(1,2); % {1} for top-view, {2} for front-view
        poleAxesMoving = {}; % axes for poles during moving. cell(nframes,2). only from poleMovingFrames 
        poleUpFrames = [];
        poleMovingFrames = [];
        angle = [];
        dist2pole = [];
        topPix = []; % frame-by-frame bottom-right pixel value in width (x-axis of the video) of the top-view pole. NaN if the pole is out of sight. 
        apUpPosition = []; % ap position during pole up
        binvavg = [];
        
        apPosition = []; % given anterior-posterior motor position during pole up frames, and estimated positions during pole moving frames.        
        radialDistance = []; 
        % to be calculated after WST is made once, collecting all data from all the trials of the same angle in the session
        mirrorAngle = []; % 
        % to be calculated after WST is made once, collecting all data from all the trials of the same angle in the session
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Follows normal image coordinates!! Different from whisker 
        % polynomial or polyfit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    properties (Dependent = true)
%         baseCoordinateTopview
%         lengthAlongWhiskerTopview
    end
    
    methods (Access = public)
        function obj = WhiskerSignalTrial_2pad(w, poleAxisTop, varargin) % w is whisker trial 2pad obj

            obj = obj@Whisker.WhiskerSignalTrial(w,varargin{:});                   
           
            obj.trackerData = w.trackerData;
            obj.whiskerPadOrigin = w.whiskerPadOrigin;
            obj.trackerFrames = w.trackerFrames;
            obj.imagePixelDimsXY = w.imagePixelDimsXY;
            obj.poleUpFrames = w.poleUpFrames;
            obj.poleMovingFrames = w.poleMovingFrames;
            obj.poleAxesUp = w.poleAxesUp;
            obj.poleAxesUp{1} = poleAxisTop;
            obj.poleAxesMoving = w.poleAxesMoving;
            for i = 1 : length(obj.poleMovingFrames)
                obj.poleAxesMoving{i,1} = poleAxisTop;
            end
            obj.angle = w.angle;
            obj.dist2pole = w.dist2pole;
            obj.topPix = w.topPix;
            obj.apUpPosition = w.apUpPosition;
            obj.radialDistance = w.radialDistance;
            obj.binvavg = w.binvavg;            
            obj.whiskerPoleIntersection = cell(obj.nof,length(obj.trajectoryIDs));
            obj.whiskerEdgeCoord = zeros(obj.nof,length(obj.trajectoryIDs));
            if ~strcmp(obj.trialType,'oo') % when it's not a catch trial 
                obj.find_whisker_pole_intersection;
            end
        end
        
        function obj = find_whisker_pole_intersection(obj)
%             Get whisker-pole intersection points.
%             If the tracked whisker is too short, then extend the whisker 
%               (extrapolate) from the tip, with the theta at tip preserved.
%             If the whisker still does not meet the pole edge, then return
%               [].
%             
%             2017/04/10 JK
% 
% 2018/07/09 JK
% Intersection is calculated from raw tracker data and pole edge, not from fitted data points since whisker tracking after the pole edge is not reliable
% So, moved from WST to WT
            if length(obj.trajectoryIDs) ~= 2
                error('Number of whisker should be 2')
            end            
            if isempty(obj.poleAxesUp) || isempty(obj.poleAxesMoving)
                return
            elseif isempty(obj.poleAxesUp{1}) || isempty(obj.poleAxesMoving{1}) || isempty(obj.poleUpFrames) || isempty(obj.poleMovingFrames)
                return
            else
                for i = 1 : 2
                    for k = 1 : obj.nof
                        if ~isempty(find(obj.trackerFrames{i} == k-1,1)) && (ismember(k,obj.poleUpFrames) || ismember(k,obj.poleMovingFrames))
                            frame_ind = find(obj.trackerFrames{i} == k-1,1);
                            xall = obj.trackerData{i}{frame_ind}{4};
                            xall = xall';
                            yall = obj.trackerData{i}{frame_ind}{5};
                            yall = yall';
                            C = [xall+1;yall+1]; % [y;x] order is changed to [x;y] for consistency 2018/06/13 JK
                            % whisker tracker's coordinate is based on zero (According to WhiskerTrial comments)

                            if ismember(k,obj.poleMovingFrames)
                                currAxis = obj.poleAxesMoving{obj.poleMovingFrames==k,i};
                            elseif ismember(k,obj.poleUpFrames)
                                currAxis = obj.poleAxesUp{i};
                            else
                                obj.whiskerPoleIntersection{k,i} = [];
                                obj.whiskerEdgeCoord(k,i) = NaN;
                                continue
                            end

                            temp = Whisker.InterX(currAxis,C); % Whisker.InterX only gets inputs as column pairs of points (x = C(1,:), y = C(2,:))
                            if ~isempty(temp)
                                temp = temp'; % row vector
                                if size(temp,1) > 1
                                    temp = sortrows(temp,-2); % sort temp descending order of the second column, which is 1st dim (or ty). Changed from 1st column to 2nd 2018/06/13 JK                                    
                                    temp = temp(1,:); % select the largest value (lowest in the video)
                                end
                                obj.whiskerPoleIntersection{k,i} = temp; % This is not changed, (for analysis of JK025~041), since the important one is just whiskerEdgeCoord, calculated symmetrically. 2018/06/13 JK
                                obj.whiskerEdgeCoord(k,i) = sqrt(sum((temp'-currAxis(:,1)).^2)); % the distances from each axis origin
                            else  % extrapolate the whisker and find the intersection with pole edge
                                q = linspace(0,1.3); % 30% stretching out. If longer is needed, it would not fit well anyway.
                                x = polyval(obj.polyFits{i}{1}(frame_ind,:), q);
                                y = polyval(obj.polyFits{i}{2}(frame_ind,:), q);
                                if size(x,1) > size(x,2)
                                    x = x'; y = y';
                                end
                                C = [x+1; y+1];
                                temp = Whisker.InterX(currAxis,C);
                                if ~isempty(temp)
                                    temp = temp'; % row vector
                                    if size(temp,1) > 1
                                        temp = sortrows(temp,-2); % sort temp descending order of the second column, which is 1st dim (or ty). Changed from 1st column to 2nd 2018/06/13 JK                                    
                                        temp = temp(1,:); % select the largest value (lowest in the video)
                                    end
                                    obj.whiskerPoleIntersection{k,i} = temp; % This is not changed, (for analysis of JK025~041), since the important one is just whiskerEdgeCoord, calculated symmetrically. 2018/06/13 JK
                                    obj.whiskerEdgeCoord(k,i) = sqrt(sum((temp'-currAxis(:,1)).^2)); % the distances from each axis origin
                                else
                                    obj.whiskerPoleIntersection{k,i} = [];
                                    obj.whiskerEdgeCoord(k,i) = NaN;
                                end
                            end
                        else
                            obj.whiskerPoleIntersection{k,i} = [];
                            obj.whiskerEdgeCoord(k,i) = NaN;
                        end
                    end
                end
            end
            
        end

        function frontRInMm = get_frontRInMm(obj, rInMm)
            frontRInMm = rInMm;
        end  
        
    end
    
    methods 
%         function value = get.baseCoordinateTopview(obj)
%             x = polyval(obj.polyFitsMask{1}{1}, linspace(0,1));
%             y = polyval(obj.polyFitsMask{1}{2}, linspace(0,1));
%             mask = [x;y];
%             value = zeros(length(obj.trackerData{1}),2);
%             for i = 1 : length(obj.trackerData{1})
%                 xall = obj.trackerData{1}{i}{4};
%                 xall = xall';
%                 yall = obj.trackerData{1}{i}{5};
%                 yall = yall';
%                 whisker = [xall+1;yall+1];
%                 temp = Whisker.InterX(whisker, mask);
%                 if ~isempty(temp)
%                     value(i,:) = temp;
%                 else
%                     value(i,:) = whisker(1,:);
%                 end
%             end
%         end
        
%         function value = get.lengthAlongWhiskerTopview(obj)
%             value = nan(length(obj.baseCoordinateTopview),1);
%             ind = find(isfinite(sum(obj.baseCoordinateTopview,2)));
%             for i = 1 : length(value)
%                 if ~isempty(obj.whiskerPoleIntersection{i,1}) && ismember(i, ind)
%                     whisker = [obj.trackerData{1}{i}{4}+1, obj.trackerData{1}{i}{5}+1];
%                     dist2base = sum((whisker - obj.baseCoordinateTopview(i,:)).^2);
%                     baseInd = find(dist2base == min(dist2base));
%                     dist2intersect = sum((whisker - obj.whiskerPoleIntersection{i,1}).^2);
%                     intersectInd = find(dist2intersect == min(dist2intersect));
%                     arcLength = [0; cumsum(sqrt(diff(whisker(:,1)).^2 + diff(whisker(:,2)).^2))];
%                     value(i) = arcLength(intersectInd) - arcLength(baseInd);
%                 end
%             end
%             value = [];
%         end
% It takes really long... do this when making w3
        
    end
       
end