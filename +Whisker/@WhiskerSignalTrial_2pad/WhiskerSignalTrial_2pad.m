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
%         trackerData = {}; 
        trackerFrames = {};
        whiskerPoleIntersection = {}; 
        whiskerEdgeCoord = [];
        imagePixelDimsXY = []; % [NumberOfXPixels NumberOfYPixels]
        poleAxesUp = cell(1,2); % {1} for top-view, {2} for front-view
        poleAxesMoving = {}; % axes for poles during moving. cell(nframes,2). only from poleMovingFrames 
%         vavg = []; % average pic
        poleUpFrames = [];
        poleMovingFrames = [];
        angle = [];
        dist2pole = [];
        topPix = []; % frame-by-frame bottom-right pixel value in width (x-axis of the video) of the top-view pole. NaN if the pole is out of sight. 
        apUpPosition = []; % ap position during pole up
        binvavg = [];
        
        apPosition = []; % given anterior-posterior motor position during pole up frames, and estimated positions during pole moving frames.        
        % to be calculated after WST is made once, collecting all data from all the trials of the same angle in the session
        mirrorAngle = []; % 
        % to be calculated after WST is made once, collecting all data from all the trials of the same angle in the session
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Follows normal image coordinates!! Different from whisker 
        % polynomial or polyfit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    properties (Dependent = true)        

    end
    
    methods (Access = public)
        function obj = WhiskerSignalTrial_2pad(w, varargin)

            obj = obj@Whisker.WhiskerSignalTrial(w,varargin{:});
            
            ntraj = length(obj.trajectoryIDs);       
% uncomment this when I need to use show_nan_polyfits_with_mask            
%             obj.trackerData = w.trackerData;
            obj.whiskerPoleIntersection = cell(obj.nof,ntraj);
            obj.whiskerEdgeCoord = zeros(obj.nof,ntraj);
            
            obj.trackerFrames = w.trackerFrames;
            obj.imagePixelDimsXY = w.imagePixelDimsXY;
            obj.poleAxesUp = w.poleAxesUp;
            obj.poleAxesMoving = w.poleAxesMoving;
            obj.poleUpFrames = w.poleUpFrames;
            obj.poleMovingFrames = w.poleMovingFrames;
            obj.angle = w.angle;
            obj.dist2pole = w.dist2pole;
            obj.topPix = w.topPix;
            obj.apUpPosition = w.apUpPosition;
            obj.binvavg = w.binvavg;
            
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

            if length(obj.trajectoryIDs) ~= 2
                error('Number of whisker should be 2')
            end            
            if isempty(obj.polyFits)
                error('obj.polyFits is empty.')
            end
            if isempty(obj.poleAxesUp) || isempty(obj.poleAxesMoving)
                return
            elseif isempty(obj.poleAxesUp{1}) || isempty(obj.poleAxesMoving{1}) || isempty(obj.poleUpFrames) || isempty(obj.poleMovingFrames)
                return
            else
                for i = 1 : 2
                    for k = 1 : obj.nof
                        if ~isempty(find(obj.trackerFrames{i} == k-1,1)) && (ismember(k,obj.poleUpFrames) || ismember(k,obj.poleMovingFrames))
                            q = linspace(0,1);
                            frame_ind = find(obj.trackerFrames{i} == k-1,1);
                            fittedX = obj.polyFits{i}{1}(frame_ind,:);
                            fittedY = obj.polyFits{i}{2}(frame_ind,:);
                            xall = polyval(fittedX,q);
                            yall = polyval(fittedY,q);
                            C = [yall+1;xall+1];
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
                                    temp = sortrows(temp,-1); % sort temp descending order of the first column, which is 1st dim (or ty)
                                    temp = temp(1,:); % select the largest value (lowest in the video)
                                end
                                obj.whiskerPoleIntersection{k,i} = temp; 
                                obj.whiskerEdgeCoord(k,i) = sqrt(sum((temp'-currAxis(:,1)).^2)); % the distances from each axis origin
                            else  % extrapolate the whisker and find the intersection with pole edge
                                % increase the whisker outward for 30% based on the polynomial fitting
                                q = linspace(0,1); 
                                xall = polyval(fittedX,q);
                                yall = polyval(fittedY,q);
                                C = [yall+1;xall+1];
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
                                        temp = sortrows(temp,-1); % sort temp descending order of the first column, which is 1st dim (or ty)
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

        function frontRInMm = get_frontRInMm(obj, rInMm)
            [~,~,y,x,~] = obj.get_theta_kappa_at_point(0,rInMm); % x and y in (1,nframes)
            [~,~,y0,x0,~] = obj.get_theta_kappa_at_point(0,0);
            ratio = zeros(1,length(x));
            for i = 1 : length(x)
                xtip = x0(i) + (x(i)-x0(i))*10;
                ytip = y0(i) + (y(i)-y0(i))*10;                
                C1 = [xtip, x0(i); ytip, y0(i)];
                C2 = [obj.poleAxesUp{1}(2,:);obj.poleAxesUp{1}(1,:)];
                
                P = Whisker.InterX(C1,C2); % Find points where whisker and mask curves intersect. Slower but more
                                           % accurate version that isn't limited in resolution by the number of
                                           % points whisker and mask are evaluated at.
                if isempty(P) % stretch the pole axis by 3X
                    C2 = [obj.poleAxesUp{1}(2,1) - (obj.poleAxesUp{1}(2,end) - obj.poleAxesUp{1}(2,1)), obj.poleAxesUp{1}(2,end) + (obj.poleAxesUp{1}(2,end) - obj.poleAxesUp{1}(2,1));
                        obj.poleAxesUp{1}(1,1) - (obj.poleAxesUp{1}(1,end) - obj.poleAxesUp{1}(1,1)), obj.poleAxesUp{1}(1,end) + (obj.poleAxesUp{1}(1,end) - obj.poleAxesUp{1}(1,1))];
                    P = Whisker.InterX(C1,C2);
                end
                if isempty(P)
                    ratio(i) = nan;
                else
                    if size(P,2) > 1   % Don't need for faster version, which handles this.
                        disp('Found more than 1 intersection of whisker and pole edge (top-view for calculating the ratio); using only first.')
                        P = P(:,1);
                    end
                    ratio(i) = sqrt((P(1) - x(i))^2 + (P(2) - y(i))^2) / sqrt((P(1) - x0(i))^2 + (P(2) - y0(i))^2); 
                end
            end
            
            [~,~,y0,x0,~] = obj.get_theta_kappa_at_point(1,0);
            frontRInMm = zeros(length(obj.time{2}),1);
            noNaNInd = ismember(obj.time{2},obj.time{1});
            for i = 1 : length(noNaNInd)
                if noNaNInd == 0
                    frontRInMm(i) = NaN;
                else                    
                    C1 = [x0(i), x0(i); obj.imagePixelDimsXY(2), obj.imagePixelDimsXY(1)];
                    C2 = [obj.poleAxesUp{2}(2,:);obj.poleAxesUp{2}(1,:)];
                    P = Whisker.InterX(C1,C2);
                    if isempty(P) % stretch the pole axis by 3X
                        C2 = [obj.poleAxesUp{1}(2,1) - (obj.poleAxesUp{1}(2,end) - obj.poleAxesUp{1}(2,1)), obj.poleAxesUp{1}(2,end) + (obj.poleAxesUp{1}(2,end) - obj.poleAxesUp{1}(2,1));
                            obj.poleAxesUp{1}(1,1) - (obj.poleAxesUp{1}(1,end) - obj.poleAxesUp{1}(1,1)), obj.poleAxesUp{1}(1,end) + (obj.poleAxesUp{1}(1,end) - obj.poleAxesUp{1}(1,1))];
                        P = Whisker.InterX(C1,C2);
                    end
                    topInd = find(obj.time{1} == obj.time{2}(i),1,'first');
                    if ~isempty(P) && ~isempty(topInd)
                        if size(P,2) > 1
                            disp('Found more than 1 intersection of whisker and pole edge (front-view); using only first.')
                            P = P(:,1);
                        end
                        rinmmEdge = [obj.poleAxesUp{2}(2,:); obj.poleAxesUp{2}(1,:) + ratio(topInd) * (y0(i) - P(2))];
                        q = linspace(0,1);
                        frontx = polyval(obj.polyFits{2}{1}(i,:), q);
                        fronty = polyval(obj.polyFits{2}{2}(i,:), q);
                        C1 = [frontx; fronty];
                        P2 = Whisker.InterX(C1, rinmmEdge);
                        if isempty(P2)                            
                            rinmmEdge = [obj.poleAxesUp{2}(2,1) - (obj.poleAxesUp{2}(2,end) - obj.poleAxesUp{2}(2,1)), obj.poleAxesUp{2}(2,end) + (obj.poleAxesUp{2}(2,end) - obj.poleAxesUp{2}(2,1));
                                        obj.poleAxesUp{2}(1,1) - (obj.poleAxesUp{2}(1,end) - obj.poleAxesUp{2}(1,1)), obj.poleAxesUp{2}(1,end) + (obj.poleAxesUp{2}(1,end) - obj.poleAxesUp{2}(1,1))];
                            P2 = Whisker.InterX(rinmmEdge, [frontx;fronty]);
                        end
                        if ~isempty(P2)
                            if size(P,2) > 1
                                disp('Found more than 1 intersection of whisker and fictional pole edge (front-view); using only first.')
                                P2 = P2(:,1);
                            end
                            
                            pxDot = polyder(frontx);
                            pyDot = polyder(fronty);

                            xDot = polyval(pxDot,q);
                            yDot = polyval(pyDot,q);

                            dq = [0 diff(q)];
                
                            % Arc length as a function of q, after integration below:
                            R = cumsum(sqrt(xDot.^2 + yDot.^2) .* dq); % arc length segments, in pixels, times dq.

                            C = C1 - repmat([x0(i);y0(i)],[1 size(C1,2)]);
                            err = sqrt(C(1,:).^2 + C(2,:).^2);
                            ind = err==min(err);
                            
                            R = R - R(ind);
                            
                            C = C1 - repmat(P2, [1 size(C1,2)]);
                            err2 = sqrt(C(1,:).^2 + C(2,:).^2);
                            ind2 = err2==min(err2);
                            frontRInMm(i) = R(ind2)/obj.pxPerMm;
                        end                        
                    end
                end
            end
        end
        
        
%         function show_nan_polyfits_with_mask(obj, varargin)

%         Don't need it for now. If I need it, I have to populate
%         trackerData from the properties.

%             % to graphically check what went wrong for NaN frames of
%             % polyFits. Currently for 2pad top and front view
%             % 2018/02/27 JK
%             num_whisker = length(obj.polyFits);
%             if length(obj.polyFitsMask) ~= num_whisker
%                 error('Number of masks does not match with number of whiskers.')
%             end
%             errorframes = [];
%             for i = 1 : num_whisker
%                 errorframes = union(errorframes, find(isnan(obj.polyFits{i}{1}(:,1))));
%             end            
%             q = linspace(0,1);            
%             figure
%             for k = 1 : length(errorframes)
%                 for i = 1 : num_whisker    
%                     hold all
%                     plot(obj.trackerData{i}{errorframes(k)}{4}, obj.trackerData{i}{errorframes(k)}{5}, 'k-')
%                     plot(polyval(obj.polyFitsMask{i}{1},q), polyval(obj.polyFitsMask{i}{2},q), 'r-')           
%                     title(['Frame # ', num2str(errorframes(k))])                    
%                 end
%                 w = 0;
%                 while w == 0
%                     w = waitforbuttonpress;
%                 end
%                 clf
%             end
%             
%         end
        
    end
    
    methods (Access = private)
       

    end
       
end