classdef Whisker3D_2pad < handle
    %
    % From WST file, reconstruct 3D whisker shape in each frame.
    % Using tracker data, not fitted data.
    % For just one whisker, imaged from top and front view simultaneously.
    % Also calculates kappa and angle at base of top-view and front-view.
    %
    % 2018/08/14 JK
    %
    properties
        trialNum = [];
        trialType = '';
        framePeriodInSec = 1/311; % 311 Hz
        pxPerMm = []; %  Inherited from WhiskerSignalTrial.
        mouseName = '';
        sessionName = '';
        trackerFileName = '';
        rInMm = [];
        
        servoAngle = []; % Inherited from WhiskerSignalTrial.
        apPosition = []; % Inherited from WhiskerSignalTrial.
        radialDistance = [] ; % Inherited from WhiskerSignalTrial.
        nof = []; % Number of Frames. Inherited from WhiskerSignalTrial.
        poleUpFrames = []; % Inherited from WhiskerSignalTrial. first timepoint is 1, not 0. 2017/04/13 JK
        poleMovingFrames = [];
        mirrorAngle = 0;
        cameraAngle = 0; % for error until JK056 (4.2 degrees)
        
        time = [];
        kappaH = []; % horizaontal kappa
        theta = []; % horizontal angle (top angle)
        kappaV = []; % vertical kappa
        phi = []; % elevation angle (front angle)
        zeta = []; % roll angle, calculated by tangent line from the mask
        
        trackerData = {}; % {n}(:,1) for anterior-posterior axis, {n}(:,2) for radial axis, and {n}(:,3) for vertical axis. Starts from the mask.        
        fit3Data = {}; % same as in trackerData, except that it's for polynomial fitting (using polyfitn by John D'Errico (https://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn)
        follicle = []; % (:,1) for anterior-posterior axis, (:,2) for radial axis, and (:,3) for vertical axis.
        
    end
    
    properties (Dependent = true)

    end
    
        
    methods (Access = public)
        function obj = Whisker3D_2pad(ws, varargin)
            if nargin==0
                return
            end
            
            p = inputParser;
            p.addRequired('ws', @(x) isa(x,'Whisker.WhiskerSignalTrial_2pad'));
            p.addParameter('rInMm', 3, @(x) isnumeric(x) && numel(x)==1 );
       
            p.parse(ws,varargin{:});
            
            obj.rInMm = p.Results.rInMm;
            
            obj.trialNum = ws.trialNum;
            obj.trialType = ws.trialType;
            obj.framePeriodInSec = ws.framePeriodInSec;
            obj.pxPerMm = ws.pxPerMm;
            obj.mouseName = ws.mouseName;
            obj.sessionName = ws.sessionName;
            obj.trackerFileName = ws.trackerFileName;            
            
            obj.servoAngle = ws.angle;
            obj.apPosition = ws.apPosition;
            obj.radialDistance = ws.radialDistance;
            obj.nof = ws.nof;
            obj.poleUpFrames = ws.poleUpFrames;
            obj.poleMovingFrames = ws.poleMovingFrames;            
            obj.mirrorAngle = ws.mirrorAngle;
            R = [cosd(-obj.mirrorAngle) -sind(-obj.mirrorAngle); sind(-obj.mirrorAngle) cosd(-obj.mirrorAngle)]; % rotation matrix in top view
            
            % Compensating for camera angle error before 2018/11/13
            if strcmp(obj.mouseName(1:2), 'JK') && str2double(obj.mouseName(3:end)) < 60
                obj.cameraAngle = 4.2;
            end
            
            obj.time = intersect(ws.time{1}, ws.time{2});
            
            % calculating 3D tracker Data
            [~,tdtopind] = ismember(obj.time, ws.time{1});
            [~,tdfrontind] = ismember(obj.time, ws.time{2});
            obj.trackerData = cell(length(tdtopind),1);
            
            wpo = ws.whiskerPadOrigin;
            vwidth = ws.imagePixelDimsXY(1);
            for i = 1 : length(tdtopind)
                if tdtopind(i) && tdfrontind(i)
                    x = ws.trackerData{1}{tdtopind(i)}{4};
                    y = ws.trackerData{1}{tdtopind(i)}{5};
                    z = ws.trackerData{2}{tdfrontind(i)}{4};
                    w = ws.trackerData{2}{tdfrontind(i)}{5};
                    whiskerTop = [x'; y'];
                    maskTop = [polyval(ws.polyFitsMask{1}{1},linspace(-0.3,1.3)); polyval(ws.polyFitsMask{1}{2},linspace(-0.3,1.3))];
                    whiskerFront = [z'; w'];
                    maskFront = [polyval(ws.polyFitsMask{2}{1},linspace(-0.3,1.3)); polyval(ws.polyFitsMask{2}{2},linspace(-0.3,1.3))];

                    Ptop = Whisker.InterX(whiskerTop, maskTop);
                    Pfront = Whisker.InterX(whiskerFront, maskFront);

                    tempData = NaN(length(x),3);
                    if ~isempty(Ptop) && ~isempty(Pfront) % only consider where tracking data intersects with the mask in both views                    
                        if sqrt(sum((wpo-[x(end) y(end)]).^2)) < sqrt(sum((wpo-[x(1) y(1)]).^2))
                            % c(q_max) is closest to whisker pad origin, so reverse the (x,y) sequence
                            x = x(end:-1:1);
                            y = y(end:-1:1);
                        end
                        for j = 1 : length(x)
                            if y(j) <= Ptop(2)                            
                                dist = Ptop(2) - y(j); % lateral distance from the mask, calculated from the top-view
                                line = [1, vwidth; Pfront(2) - dist, Pfront(2) - dist]; % corresponding line for front-view                            
                                P = Whisker.InterX(line, whiskerFront);
                                if size(P,2) == 1 % else, there is no intersection or more than 1 intersection, which in that case cannot correctly reconstruct 3D shape                                
                                    % projection to the axis orthogonal to the body axis. Height (j,3) does not have to change
                                    v = R * [x(j); y(j)];
                                    tempData(j,1) = v(1);
                                    tempData(j,2) = v(2);
                                    tempData(j,3) = P(1);
                                end
                            end
                        end
                        obj.trackerData{i} = tempData(isfinite(sum(tempData,2)),:);
                    end
                end
            end            
            ind = find(cellfun(@(x) length(x), obj.trackerData));
            obj.follicle = zeros(length(ind), 3);
            if length(ind) < length(obj.time)
                obj.time = obj.time(ind);
                tempData = obj.trackerData;
                obj.trackerData = cell(length(ind),1);
                for i = 1 : length(ind)
                    obj.trackerData{i} = tempData{ind(i)};
                    obj.follicle(i,:) = obj.trackerData{i}(1,:);
                end
            end
            
            % Calculating whisker kinematics
            % 1. find the calculation point for kappas 
            obj.kappaH = nan(length(ind),1);
            obj.kappaV = nan(length(ind),1);
            obj.theta = nan(length(ind),1);
            obj.phi = nan(length(ind),1);
            
            for fi = 1 : length(ind)
                preDist = (obj.rInMm-1) * obj.pxPerMm;
                postDist = (obj.rInMm+1) * obj.pxPerMm;
                whiskerPixLengths = zeros(length(obj.trackerData{fi}),1);
                for i = 2 : length(whiskerPixLengths)
                    whiskerPixLengths(i) = sqrt(sum((obj.trackerData{fi}(i,:) - obj.trackerData{fi}(i-1,:)).^2));
                end            
                prePoint = find(cumsum(whiskerPixLengths) > preDist, 1, 'first');
                postPoint = find(cumsum(whiskerPixLengths) > postDist, 1, 'first');

                if ~isempty(prePoint) && ~isempty(postPoint) && postPoint - prePoint > 3
                    x = obj.trackerData{fi}(prePoint:postPoint,1);
                    y = obj.trackerData{fi}(prePoint:postPoint,2);
                    z = obj.trackerData{fi}(prePoint:postPoint,3);

                    q = linspace(0,1, postPoint-prePoint+1);
                    px = polyfit(q',x,2);
                    py = polyfit(q',y,2);
                    pz = polyfit(q',z,2);

                    % horizontal & vertical kappa
                    pxDot = polyder(px);
                    pxDoubleDot = polyder(pxDot);

                    pyDot = polyder(py);
                    pyDoubleDot = polyder(pyDot);

                    pzDot = polyder(pz);
                    pzDoubleDot = polyder(pzDot);

                    xDot = polyval(pxDot,q);
                    xDoubleDot = polyval(pxDoubleDot,q);

                    yDot = polyval(pyDot,q);
                    yDoubleDot = polyval(pyDoubleDot,q);

                    zDot = polyval(pzDot,q);
                    zDoubleDot = polyval(pzDoubleDot,q);

                    kappasH = (xDot.*yDoubleDot - yDot.*xDoubleDot) ./ ((xDot.^2 + yDot.^2).^(3/2)) * obj.pxPerMm; % SIGNED CURVATURE, in 1/mm.
                    kappasV = (zDot.*yDoubleDot - yDot.*zDoubleDot) ./ ((zDot.^2 + yDot.^2).^(3/2)) * obj.pxPerMm; % SIGNED CURVATURE, in 1/mm.
                    midpoint = round((postPoint - prePoint)/2);
                    obj.kappaH(fi) = kappasH(midpoint);
                    obj.kappaV(fi) = kappasV(midpoint);
                end
            end
            
            obj.theta = ws.get_theta_at_base(0) - obj.mirrorAngle;
            obj.phi = ws.get_theta_at_base(1) - obj.cameraAngle;
            
            % 3D fitting of the tracked Data
            obj.fit3Data = cell(length(obj.trackerData),1);
%             for i = 1 : length(obj.trackerData)
%                 obj.fit3Data{i} = polyfitn(obj.trackerData{1});
%             end
            
        end
        
        function show_all_3D(obj)
            figure, hold on,
            for i = 1 : length(obj.trackerData)
                plot3(obj.trackerData{i}(:,1), obj.trackerData{i}(:,2), obj.trackerData{i}(:,3), '-')
            end
        end
        
        function show_onebyone_3D(obj)
            figure, 
            i = 1; plot3(obj.trackerData{i}(:,1), obj.trackerData{i}(:,2), obj.trackerData{i}(:,3), 'k-'); hold on, 
            plot3(obj.trackerData{i}(1,1), obj.trackerData{i}(1,2), obj.trackerData{i}(1,3), 'r.', 'markersize', 20); 
            hold off, 
            view(3);
            [az, el] = view;
            xl = [min(cellfun(@(x) min(x(:,1)), obj.trackerData)), max(cellfun(@(x) max(x(:,1)), obj.trackerData))];
            yl = [min(cellfun(@(x) min(x(:,2)), obj.trackerData)), max(cellfun(@(x) max(x(:,2)), obj.trackerData))]; 
            zl = [min(cellfun(@(x) min(x(:,3)), obj.trackerData)), max(cellfun(@(x) max(x(:,3)), obj.trackerData))];
            while( i > 0 )                
                plot3(obj.trackerData{i}(:,1), obj.trackerData{i}(:,2), obj.trackerData{i}(:,3), 'k-', 'linewidth', 2), hold on
                if i < 51
                    j = 1;
                else
                    j = i - 50;
                end
                for k = j : i-1
                    plot3(obj.trackerData{k}(:,1), obj.trackerData{k}(:,2), obj.trackerData{k}(:,3), 'color', [0.7 0.7 0.7])
                end
                if i > length(obj.trackerData) - 50
                    j = length(obj.trackerData);
                else
                    j = i + 50;
                end
                for k = i+1:j
                    plot3(obj.trackerData{k}(:,1), obj.trackerData{k}(:,2), obj.trackerData{k}(:,3), 'color', [0.7 0.7 0.7])                    
                end                
                plot3(obj.trackerData{i}(:,1), obj.trackerData{i}(:,2), obj.trackerData{i}(:,3), 'k-', 'linewidth', 4),
                plot3(obj.trackerData{i}(1,1), obj.trackerData{i}(1,2), obj.trackerData{i}(1,3), 'r.', 'markersize', 25),
                hold off
                title({'Navigate using keyboard'; ['Trial # ', obj.trackerFileName]; ['Frame # ', num2str(round(obj.time(i)/obj.framePeriodInSec))]});
                xlabel('Rostro-caudal'), ylabel('Medio-lateral'), zlabel('Dorso-ventral')                
                set(gca, 'linewidth', 3, 'fontweight', 'bold', 'fontsize', 15, 'style', 'equal')
                view(az,el), xlim(xl), ylim(yl), zlim(zl);
                [i, az, el] = obj.keyboard_navigation_3d(i, length(obj.trackerData));
            end
        end
        
        function [next, az, el] = keyboard_navigation_3d(obj, curr, maxnum)
            w = 0;
            while w < 1
                w = waitforbuttonpress;
            end
            if w
                value = double(get(gcf,'CurrentCharacter'));
                [az, el] = view;
                switch value
                    case 28 % <-
                        if curr == 1
                            next = maxnum;
                        else
                            next = curr - 1;
                        end

                    case 29 % ->
                        if curr == maxnum
                            next = 1;
                        else
                            next = curr + 1;
                        end

                    case 30 % up arrow
                        if curr > maxnum - 10
                            next = curr + 10 - maxnum;
                        else
                            next = curr + 10;
                        end

                    case 31 % down arrow
                        if curr < 10
                            next = maxnum - 10 + curr;
                        else
                            next = curr - 10;
                        end

                    case 93 % ]
                        if curr > maxnum - 100
                            next = curr + 100 - maxnum;
                        else
                            next = curr + 100;
                        end

                    case 91 % [
                        if curr < 100
                            next = maxnum - 100 + curr;
                        else
                            next = curr - 100;
                        end

                    case 27 % when pressing esc                
                        next = 0;
                end
            end
        end
    end
end
        