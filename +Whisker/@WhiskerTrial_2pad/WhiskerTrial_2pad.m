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
        nof = []; % number of frames (video, not tracker)
        dist2pole = [];
        binvavg = []; % binary average image of the video during pole up frames
        
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

            p.parse(tracker_file_name, trial_num, trajectory_nums, varargin{:});
            
            obj = obj@Whisker.WhiskerTrial(p.Results.tracker_file_name, p.Results.trial_num, p.Results.trajectory_nums, 'mouseName', p.Results.mouseName, 'sessionName', p.Results.sessionName, 'trialType', p.Results.trialType);
            obj.angle = p.Results.angle;
            obj.apUpPosition = p.Results.apUpPosition;
            [obj.nof, obj.poleUpFrames, obj.poleMovingFrames, obj.poleAxesUp, obj.poleAxesMoving, obj.topPix, obj.barPos, obj.binvavg] = Whisker.pole_edge_detection(obj.trackerFileName, obj.angle, obj.barRadius);
            if ~isempty(obj.barPos)
                obj.dist2pole = obj.distance_to_pole; 
            end
        end
        
                
        function dist = distance_to_pole(obj)
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
        
    end
    
    methods (Access = private)
       

    end
       
end