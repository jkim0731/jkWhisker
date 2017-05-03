classdef WhiskerSignalTrial_2pad < Whisker.WhiskerSignalTrialI
   
    % Infer whisker-pole contact frames (index and time) from angle &
    % radial distance tasks (2pad, by JK). The calculation is based on the
    % intersection points between whisker and the pole edge, from both the
    % front and top view. This is based on the assumption that whisker-pole
    % contact points are determined uniquely when combining both of the
    % views. 
    %
    % 2017/04/10 JK
    
    properties
        contact_ind = {};
        contact_time = {};        
        trackerData = {};
        whisker_pole_intersection = {}; 
        whisker_edge_coord = [];
        imagePixelDimsXY = [400 250]; % [NumberOfXPixels NumberOfYPixels]
        pole_edge = cell(1,2); % edge detection of the pole
        pole_axes = cell(1,2); % axes for edges
        vavg = []; % average pic
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Follows normal image coordinates!! Different from whisker 
        % polynomial or polyfit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    properties (Dependent = true)        

    end
    
    methods (Access = public)
        function obj = WhiskerSignalTrial_2pad(w, varargin)

            obj = obj@Whisker.WhiskerSignalTrialI(w,varargin{:});
            
            ntraj = length(obj.trajectoryIDs);
            nframes = size(obj.polyFits{1}{1},1);
            obj.contact_ind = cell(nframes,ntraj);
            obj.contact_time = cell(nframes,ntraj);
            obj.trackerData = w.trackerData;
            obj.whisker_pole_intersection = cell(nframes,ntraj);
            obj.whisker_edge_coord = zeros(nframes,ntraj);
            obj.imagePixelDimsXY = w.imagePixelDimsXY;
            [obj.pole_edge, obj.pole_axes, obj.vavg] = Whisker.pole_edge_detection(obj.trackerFileName);
            obj.find_whisker_pole_intersection;            
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

            npoints = 100;
            if length(obj.trajectoryIDs) ~= 2
                error('Number of whisker should be 2')
            end            
            nframes = size(obj.polyFits{1}{1},1);            
            if isempty(obj.polyFits)
                error('obj.polyFits is empty.')
            end
            if isempty(obj.pole_edge{1}) || isempty(obj.pole_edge{2})
                [obj.pole_edge, obj.pole_axes, obj.vavg] = Whisker.pole_edge_detection(obj.trackerFileName);
            end
            
            for i = 1 : 2
%                 fittedX = obj.polyFits{i}{1};
%                 fittedY = obj.polyFits{i}{2};
%                 q = linspace(0,1,npoints);
                
                for k=1:nframes
                    tx = obj.trackerData{i}{k}{4};
                    ty = obj.trackerData{i}{k}{5};
%                     px = fittedX(k,:);
%                     py = fittedY(k,:);
                    C = [ty'+1;tx'+1]; % converting whisker tracker points into normal MATLAB coordinates
%                     C = [polyval(py,q)+1; polyval(px,q)+1]; % converting whisker tracker points into normal MATLAB coordinates
                    
%                     try
                        temp = Whisker.InterX(obj.pole_axes{i},C); % Whisker.InterX only gets inputs as column pairs of points (x = C(1,:), y = C(2,:))                        
                        if ~isempty(temp)
                            temp = temp'; % row vector
                            if size(temp,1) > 1
                                temp = sortrows(temp,-1); % sort temp descending order of the first column, which is 1st dim (or ty)
                                temp = temp(1,:); % select the largest value (lowest in the video)
                            end
                            obj.whisker_pole_intersection{k,i} = temp; 
                            obj.whisker_edge_coord(k,i) = sqrt(sum((temp'-obj.pole_axes{i}(:,1)).^2)); % the distances from each axis origin
                        else  % extrapolate the whisker and find the intersection with pole edge
                            % polyfit from the last 4 points, linear fitting, and then drawing from the 5th last point
                            % Extend for 20 pixels y direction (negative 1st dim in MATLAB convention). Should be enough. If that's not enough, it means extrapolation itself cannot be correct anyway.
                            if ty(1) < ty(end) % follicle at the beginning of the vector (column)
                                ty = flip(ty);
                                tx = flip(tx);
                            end
                            p = polyfix(ty(end-3:end-1),tx(end-3:end-1),1,ty(end-3),tx(end-3)); % I need p(1) only.
                            tip = [ty(end-3)+1, tx(end-3)+1];
                            ext_tip = [tip(1)-20, tip(2)-p(1)*20];
%                             tip = C(:,end)'; % row vector
%                             tip_1back = C(:,end-1)'; % row vector
%                             tip_5back = C(:,end-5)';
%                             delta = tip_5back - tip_1back; % row vector; the angle at the tip, -5 to -1 end, since the end is erroneous
%                             if delta(1) > 0 
%                                 steps = (tip(1)-1)/delta(1);
%                                 x_intersect = tip(2) - steps * delta(2);
%                                 ext_tip = [1, x_intersect]; % extrapolated tip. InterX works for negative values too.
%                             elseif delta(1) == 0
%                                 ext_tip = [1, tip(2)];
%                             else % very unlikely that this will happen, but just in case...
%                                 steps = (tip(1) - obj.imagePixelDimsXY(2))/delta(1);
%                                 x_intersect = tip(2) - steps * delta(2);
%                                 ext_tip = [obj.imagePixelDimsXY(2), x_intersect];
%                             end
                            L = [tip', ext_tip'];
                            temp = Whisker.InterX(obj.pole_axes{i},L);                            
                            if ~isempty(temp)
                                temp = temp'; % row vector
                                if size(temp,1) > 1
                                    temp = sortrows(temp,-1); % sort temp descending order of the first column, which is 1st dim (or ty)
                                    temp = temp(1,:); % select the largest value (lowest in the video)
                                end                            
                                obj.whisker_pole_intersection{k,i} = temp; 
                                obj.whisker_edge_coord(k,i) = sqrt(sum((temp'-obj.pole_axes{i}(:,1)).^2)); % the distances from each axis origin
                            else
                                obj.whisker_pole_intersection{k,i} = [];
                                obj.whisker_edge_coord(k,i) = NaN;
                            end
                        end
%                     catch
%                         obj.whisker_pole_intersection{k,i} = [];
%                         obj.whisker_edge_coord{k,i} = [];
%                     end
                end
            end
        end
        
    end
    
    methods (Access = private)
       

    end
       
end