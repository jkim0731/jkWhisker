classdef WhiskerTrial_2pad < Whisker.WhiskerTrial
   
    % Get pole edge and axis information, and estimated pole position (anterior-posterio), all frame-by-frame.
    % Pole available frames and pole moving frames.
    % 2018/04/12 JK   
    
    properties
        poleUpFrames = [];
        poleMovingFrames = [];
        poleAxesUp = cell(1,2); % {1} for top-view, {2} for front-view
        poleAxesMoving = {}; % axes for poles during moving. cell(nframes,2). only from poleMovingFrames
        topPix = []; % frame-by-frame bottom-right pixel value in width (x-axis of the video) of the top-view pole. NaN if the pole is out of sight. 
        angle = [];
        apUpPosition = [];
        radialDistance = [];
        nof = []; % number of frames (video, not tracker)
        dist2pole = [];
        binvavg = []; % binary average image of the video during pole up frames
        
        whiskerPoleIntersection = {}; 
        whiskerEdgeCoord = [];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Follows normal image coordinates!! Different from whisker 
        % polynomial or polyfit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    properties (Dependent = true)        

    end
    
    methods (Access = public)
        function obj = WhiskerTrial_2pad(tracker_file_name, trial_num, trajectory_nums, varargin)
            
            p = inputParser;

            p.addRequired('tracker_file_name', @ischar);
            p.addRequired('trial_num', @isnumeric);
            p.addRequired('trajectory_nums', @isnumeric);
            p.addParameter('mouseName', '', @ischar);
            p.addParameter('sessionName', '', @ischar);
            p.addParameter('trialType', '', @ischar);
            p.addParameter('angle', [], @isnumeric);
            p.addParameter('apUpPosition', [], @isnumeric);
            p.addParameter('radialDistance', [], @isnumeric);

            p.parse(tracker_file_name, trial_num, trajectory_nums, varargin{:});
            
            obj = obj@Whisker.WhiskerTrial(p.Results.tracker_file_name, p.Results.trial_num, p.Results.trajectory_nums, 'mouseName', p.Results.mouseName, 'sessionName', p.Results.sessionName, 'trialType', p.Results.trialType);
            obj.angle = p.Results.angle;
            obj.apUpPosition = p.Results.apUpPosition;
            obj.radialDistance = p.Results.radialDistance;
            if ~strcmp(p.Results.trialType, 'oo') && ~contains(p.Results.sessionName, 'piezo') && ~contains(p.Results.sessionName, 'spont') && ~isempty(obj.angle)
                [obj.nof, obj.poleUpFrames, obj.poleMovingFrames, obj.poleAxesUp, obj.poleAxesMoving, obj.topPix, obj.barPos, obj.binvavg] = Whisker.pole_edge_detection(obj.trackerFileName, obj.angle, obj.barRadius);
            end
            if ~isempty(obj.barPos) % only for 90 degrees
                obj.dist2pole = obj.distance_to_pole; 
            end
            ntraj = length(obj.trajectoryIDs);
            obj.whiskerPoleIntersection = cell(obj.nof,ntraj);
            obj.whiskerEdgeCoord = zeros(obj.nof,ntraj);
            if ~strcmp(obj.trialType,'oo') % when it's not a catch trial 
                obj.find_whisker_pole_intersection;
            end
            
        end
        
                
        function dist = distance_to_pole(obj) % only for 90 degrees
            dist = NaN(obj.nof,1);
            for i = 1 : size(obj.barPos,1) % frames with bar position
                frameNum = obj.barPos(i,1);
                if ~isempty(find(obj.trackerFrames{1} == frameNum-1,1))
                    frame_ind = find(obj.trackerFrames{1} == frameNum-1,1);
                    f = obj.trackerData{1}{frame_ind};
                    if numel(f{4}) > 1 % Tracker can sometimes (rarely) leave frame entries for a trajectory in whiskers file that have no pixels.
                        x = f{4};
                        if strcmp(obj.trackerFileFormat,'whisker0')
                            x = (x(1):x(2))';
                        end
                        y = f{5};
                        if size(x,1) > size(x,2)
                            x = x';
                        end
                        if size(y,1) > size(y,2)
                            y = y';
                        end                                           
                        distance = (obj.barPos(i,2) - x).^2 + (obj.barPos(i,3) - y).^2;
                        % distance2pole = distance from bar center to the
                        % line between 2 closest points
                        [~, distInd] = sort(distance);
                        closest = [x(distInd(1)), y(distInd(1)), 0]; % z axis value set to 0 for cross product calculation
                        nextclosest = [x(distInd(2)), y(distInd(2)), 0];
                        barcenter = [obj.barPos(i,2), obj.barPos(i,3), 0];
                        point2line = norm(cross(closest-nextclosest,barcenter - closest)) / norm(closest-nextclosest);
                        dist(obj.barPos(i,1)) = max(point2line - obj.barRadius,0);
                    end
                end
            end
        end
        
        function obj = find_whisker_pole_intersection(obj)
%             Get whisker-pole intersection points.
%             If the tracked whisker is too short, then extend the whisker 
%               (extrapolate) from the tip, with the theta at tip preserved.
%             If the whisker still does not meet the pole edge, then return
%               [].
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Follows MATLAB image coordinates!! Different from whisker polynomial or
% polyfit (x and y axis, [0,0] vs [1,1] as the origin (left top corner)

% Currently only for face at the bottom and rightward whisking position.
% 2017/04/10 JK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%                             q = linspace(0,1);                            
%                             fittedX = obj.polyFits{i}{1}(frame_ind,:);
%                             fittedY = obj.polyFits{i}{2}(frame_ind,:);
%                             xall = polyval(fittedX,q);
%                             yall = polyval(fittedY,q);
%                             
                            xall = obj.trackerData{i}{frame_ind}{4};
                            xall = xall';
                            yall = obj.trackerData{i}{frame_ind}{5};
                            yall = yall';
                            
                            C = [xall+1;yall+1]; % [y;x] order is changed to [x;y] for consistency 2018/06/13 JK
    %                         C = [ty'+1;tx'+1]; % converting whisker tracker points into normal MATLAB coordinates

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
                                % increase the whisker outward for 30% based on the polynomial fitting
                                polyDegree = 5;
                                wpo = obj.whiskerPadOrigin;
                                if sqrt(sum((wpo-[xall(end) yall(end)]).^2)) < sqrt(sum((wpo-[xall(1) yall(1)]).^2))
                                    % c(q_max) is closest to whisker pad origin, so reverse the (x,y) sequence
                                    xall = xall(end:-1:1);
                                    yall = yall(end:-1:1);
                                end
                                coeffX = Whisker.polyfit(linspace(0,1,length(xall)), xall, polyDegree);
                                coeffY = Whisker.polyfit(linspace(0,1,length(yall)), yall, polyDegree);
                                q = linspace(0,1.3); 
                                xall = polyval(coeffX,q);
                                yall = polyval(coeffY,q);
                                C = [xall+1;yall+1];
                                temp = Whisker.InterX(currAxis,C);
    %                             if ty(1) < ty(end) % follicle at the beginning of the vector (column)
    %                                 ty = flip(ty);
    %                                 tx = flip(tx);
    %                             end
    %                             p = polyfix(ty(end-3:end-1),tx(end-3:end-1),1,ty(end-3),tx(end-3)); % I need p(1) only.
    %                             tip = [ty(end-3)+1, tx(end-3)+1];
    %                             ext_tip = [tip(1)-20, tip(2)-p(1)*20];
    %                             L = [tip', ext_tip'];
    %                             temp = Whisker.InterX(obj.pole_axes{i},L);                            
                                if ~isempty(temp)
                                    temp = temp'; % row vector
                                    if size(temp,1) > 1
                                        temp = sortrows(temp,-2); % sort temp descending order of the first column, which is 1st dim (or ty) (This is the comment before the change 2018/06/13)
                                        temp = temp(1,:); % select the largest value (lowest in the video)
                                    end                            
                                    obj.whiskerPoleIntersection{k,i} = temp; 
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
        
    end
    
    methods (Access = private)
       

    end
       
end