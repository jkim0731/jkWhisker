classdef Whisker3D_2pad < handle
    %
    % From WST file, reconstruct 3D whisker shape in each frame.
    % Using tracker data, not fitted data.
    % For just one whisker, imaged from top and front view simultaneously.
    % Also calculates kappa and angle at base of top-view and front-view.
    %
    % 2018/08/14 JK
    %
    % 1. 3D mask, dome-like shape instead of palisades -> mask3d
    % 2. Kv calculated at the projection to the angle tangent to the whisker
    % at the calculated point.
    % 3. Remove double mask-crossing fragments.
    % 2019/08/31 JK
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
        kappaV = []; % vertical kappa, projected to the angle 
        phi = []; % elevation angle (front angle)
        zeta = []; % roll angle, calculated by tangent line from the mask
        alpha = []; % the angle of the tangent line of the whisker segment used for kappa calculation
        
        trackerData = {}; % {n}(:,1) for anterior-posterior axis, {n}(:,2) for radial axis, and {n}(:,3) for vertical axis. Starts from the follicle (closest point to the face)        
        mask3d = [];
        fit3Data = {}; % same as in trackerData, except that it's for polynomial fitting (using polyfitn by John D'Errico (https://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn)
        fitorder = 5;
        base = []; % Base of the whisker. It is the point where whisker crosses with the mask. (:,1) for anterior-posterior axis, (:,2) for radial axis, and (:,3) for vertical axis.
        baseInd = []; % index of the base in each trackerData frame
%         fitBase = []; % base position with fitted 3D data
        intersectPoint = []; % index of the whisker-pole intersection
        prePoint = []; % for visual confirmation
        postPoint = [];

        lengthAlongWhisker = [];
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
            rotMatMirror = [cosd(obj.mirrorAngle) -sind(obj.mirrorAngle); sind(obj.mirrorAngle) cosd(obj.mirrorAngle)]; % rotation matrix in top view
            rotMatMirrorInv = [cosd(-obj.mirrorAngle) -sind(-obj.mirrorAngle); sind(-obj.mirrorAngle) cosd(-obj.mirrorAngle)]; % rotation matrix in top view
            % Caution: the rotation is based on 'plot function coordinate
            % system', which in 'imshow coordiate' should be in different
            % direction, hence omitting the negative sign 2018/12/10 JK
            
            vwidth = ws.imagePixelDimsXY(1);
            vheight = ws.imagePixelDimsXY(2);
            rotRefP = [vwidth/2; vheight/2];
            
            % Compensating for camera angle error before 2018/11/13
            if strcmp(obj.mouseName(1:2), 'JK') && str2double(obj.mouseName(3:end)) < 60
                obj.cameraAngle = 4.2;
            end
            rotMatCamera = [cosd(-obj.cameraAngle) -sind(-obj.cameraAngle); sind(-obj.cameraAngle) cosd(-obj.cameraAngle)];
            obj.time = intersect(ws.time{1}, ws.time{2});
            
            %% calculating 3D mask
            % basic settings
            tempMaskZ = polyval(ws.polyFitsMask{1}{1}, linspace(-0.3, 1.3));
            zInterval = max(tempMaskZ) - min(tempMaskZ); % in pixels
            numPoints = round(zInterval/ws.pxPerMm*33)+1; % total num points in defining the 3d mask is going to be numPoints^2
            % 30 um interval.

            % calculate a reference line at the top-down view (phi = 0).
            maskTopYtemp = polyval(ws.polyFitsMask{1}{2}, linspace(-0.3, 1.3, numPoints));
            maskTopXtemp = polyval(ws.polyFitsMask{1}{1}, linspace(-0.3, 1.3, numPoints));
            maskTopYapexOG = min(maskTopYtemp);
            
            % find the apex point of front-view mask (in y direction)
            maskFrontYtemp = polyval(ws.polyFitsMask{2}{2}, linspace(-0.3, 1.3, numPoints));
            maskFrontZtemp = polyval(ws.polyFitsMask{2}{1}, linspace(-0.3, 1.3, numPoints));
            maskFront = (rotMatCamera * [maskFrontZtemp - rotRefP(1); maskFrontYtemp - rotRefP(2)] + rotRefP)';
            maskFrontZ = maskFront(:,1);
            maskFrontY = maskFront(:,2);
            apexFrontY = min(maskFrontY);
            
            % along the z axis, move the reference line following the relationship in
            % front mask (maskFrontY and maskFrontz)
            obj.mask3d = zeros(numPoints^2,3);
            for i = 1 : numPoints
                obj.mask3d((i-1)*numPoints+1 : i*numPoints,1:2) = (rotMatMirror*([maskTopXtemp; maskTopYtemp + (maskFrontY(i) - apexFrontY)] - rotRefP) + rotRefP)';
                obj.mask3d((i-1)*numPoints+1 : i*numPoints,3) = maskFrontZ(i);
            end

            %% calculating 3D tracker Data
            [~,tdtopind] = ismember(obj.time, ws.time{1});
            [~,tdfrontind] = ismember(obj.time, ws.time{2});
            obj.trackerData = cell(length(tdtopind),1);
            obj.base = zeros(length(obj.trackerData), 3);
            obj.baseInd = zeros(length(obj.trackerData), 1);
            obj.intersectPoint = zeros(length(obj.trackerData), 3);
            obj.fit3Data = cell(length(obj.trackerData),1);
%             obj.fitBase = zeros(length(obj.trackerData), 3);
            wpo = ws.whiskerPadOrigin;
            
            for i = 1 : length(tdtopind)
                if ws.time{1}(tdtopind(i)) == ws.time{2}(tdfrontind(i)) % just in case
                    z = polyval(ws.polyFits{2}{1}(tdfrontind(i),:), linspace(0,1));
                    w = polyval(ws.polyFits{2}{2}(tdfrontind(i),:), linspace(0,1));
                    if sqrt(sum((wpo-[z(end) w(end)]).^2)) < sqrt(sum((wpo-[z(1) w(1)]).^2))
                        % c(q_max) is closest to whisker pad origin, so reverse the (z,w) sequence
                        z = z(end:-1:1);
                        w = w(end:-1:1);
                    end
                    whiskerFrontOG = [z; w];
                    PfrontOG = Whisker.InterX(whiskerFrontOG, [maskFrontZtemp; maskFrontYtemp]);
                    
                    if isempty(PfrontOG)
                        z = polyval(ws.polyFits{2}{1}(tdfrontind(i),:), linspace(-0.1, 1.1));
                        w = polyval(ws.polyFits{2}{2}(tdfrontind(i),:), linspace(-0.1, 1.1));
                        if sqrt(sum((wpo-[z(end) w(end)]).^2)) < sqrt(sum((wpo-[z(1) w(1)]).^2))
                            % c(q_max) is closest to whisker pad origin, so reverse the (z,w) sequence
                            z = z(end:-1:1);
                            w = w(end:-1:1);
                        end
                        whiskerFrontOG = [z;w];
                        PfrontOG = Whisker.InterX(whiskerFrontOG, [maskFrontZtemp; maskFrontYtemp]);
                    end

                    if ~isempty(PfrontOG) % only consider where tracking data intersects with the mask in both views
                        PfrontOG = PfrontOG(:,end); % just in case where there is 2 intersection points.
                        Pfront = rotMatCamera*(PfrontOG - rotRefP) + rotRefP;
                        
                        % top-view mask is picked from 3D mask, based on the z-axis intersection value.
                        whiskerBaseZ = Pfront(1);
                        [~,maskZind] = min(abs(maskFrontZ - whiskerBaseZ));
                        maskZ = maskFrontZ(maskZind);
                        mask3dInd = find(obj.mask3d(:,3) == maskZ);
                        maskTop = obj.mask3d(mask3dInd,1:2);
                        maskTop = maskTop';
                        
                        x = polyval(ws.polyFits{1}{1}(tdtopind(i),:), linspace(0, 1));
                        y = polyval(ws.polyFits{1}{2}(tdtopind(i),:), linspace(0, 1));
                        if sqrt(sum((wpo-[x(end) y(end)]).^2)) < sqrt(sum((wpo-[x(1) y(1)]).^2))
                            % c(q_max) is closest to whisker pad origin, so reverse the (x,y) sequence
                            x = x(end:-1:1);
                            y = y(end:-1:1);
                        end
                        whiskerTop = rotMatMirror * ([x; y] - rotRefP) + rotRefP;
                        Ptop = Whisker.InterX(whiskerTop, maskTop);
                    
                        if isempty(Ptop)
                            x = polyval(ws.polyFits{1}{1}(tdtopind(i),:), linspace(-0.1, 1.1))';
                            y = polyval(ws.polyFits{1}{2}(tdtopind(i),:), linspace(-0.1, 1.1))';
                            if sqrt(sum((wpo-[x(end) y(end)]).^2)) < sqrt(sum((wpo-[x(1) y(1)]).^2))
                                % c(q_max) is closest to whisker pad origin, so reverse the (x,y) sequence
                                x = x(end:-1:1);
                                y = y(end:-1:1);
                            end
                            whiskerTop = rotMatMirror * ([x; y] - rotRefP) + rotRefP;
                            Ptop = Whisker.InterX(whiskerTop, maskTop);
                        end
                        if ~isempty(Ptop) % only consider where tracking data intersects with the mask in both views
                            Ptop = Ptop(:,end);
                            tempData = NaN(length(x),3);
                    
                            tempData(:,1:2) = whiskerTop'; % rotate in regard to the base (intersection between whisker and mask)
                            for j = 1 : length(x)
                                lateralDist = y(j) - maskTopYapexOG; % lateral distance from the mask, calculated from the top-view, before rotation
                                line = [1, vwidth; apexFrontY + lateralDist, apexFrontY + lateralDist]; % corresponding line for front-view, after rotation
                                whiskerFront = rotMatCamera*(whiskerFrontOG-rotRefP) + rotRefP;
                                P = Whisker.InterX(line, whiskerFront);
%                                 if size(P,2) == 1 % else, there is no intersection or more than 1 intersection, which in that case cannot correctly reconstruct 3D shape
                                if size(P,2) > 0
                                    % projection to the axis orthogonal to the body axis. Height (j,3) does not have to change
                                    P = P(:,end);
                                    tempData(j,3) = P(1);
                                end
                            end
                            distFromBase = sum((Ptop - whiskerTop).^2);
                            [~,baseInd] = nanmin(distFromBase);
                            finiteInds = find(isfinite(sum(tempData,2)));
                            [~, obj.baseInd(i)] = min(abs(finiteInds - baseInd));
                            tempData = tempData(finiteInds,:);
                            obj.base(i,:) = tempData(obj.baseInd(i),:);
                            s = cumsum(sqrt([0; diff(tempData(:,1))].^2 + [0; diff(tempData(:,2))].^2) + [0; diff(tempData(:,3))].^2);
                            q = s ./ max(s);

                            if max(s) > obj.rInMm * obj.pxPerMm % 3D tracking should be at least rInMm long
                                obj.trackerData{i} = tempData;
                                px = polyfit(q',tempData(:,1)',obj.fitorder);
                                py = polyfit(q',tempData(:,2)',obj.fitorder);
                                pz = polyfit(q',tempData(:,3)',obj.fitorder);
                                obj.fit3Data{i} = [px', py', pz'];
                            end
                            frameNum = round(ws.time{1}(tdtopind(i))/ws.framePeriodInSec + 1);
                            if ~isempty(ws.whiskerPoleIntersection{frameNum,1}) && ~isempty(ws.whiskerPoleIntersection{frameNum,2})
                                
                                obj.intersectPoint(i,1:2) = (rotMatMirror * (ws.whiskerPoleIntersection{frameNum,1}' - rotRefP) + rotRefP)';
                                % y axis is inverted for ws.whiskerPolerIntersection
                                tempFrontIntersect = rotMatCamera * (ws.whiskerPoleIntersection{frameNum,2}' - rotRefP) + rotRefP;
                                obj.intersectPoint(i,3) = tempFrontIntersect(1);
                            end
                        end
                    end
                else
                    error('Frame mismatch between top and front view') % just in case    
                end
            end            
            ind = find(cellfun(@(x) length(x), obj.trackerData));
            
            if length(ind) < length(obj.time)
                obj.time = obj.time(ind);
                tempData = obj.trackerData;
                tempFit = obj.fit3Data;
                tempBase = obj.base;
                tempBaseInd = obj.baseInd;
                tempIntersectPoint = obj.intersectPoint;
                obj.trackerData = cell(length(ind),1);
                obj.fit3Data = cell(length(ind),1);
                obj.base = zeros(length(ind),3);
                obj.baseInd = zeros(length(ind),1);
                obj.intersectPoint = zeros(length(ind),3);
                for i = 1 : length(ind)
                    obj.trackerData{i} = tempData{ind(i)};
                    obj.fit3Data{i} = tempFit{ind(i)};
                    obj.base(i,:) = tempBase(ind(i),:);
                    obj.baseInd(i) = tempBaseInd(ind(i));
                    obj.intersectPoint(i,:) = tempIntersectPoint(ind(i),:);
                end
            end
            
            % Calculating whisker kinematics
            obj.kappaH = nan(length(ind),1);
            obj.kappaV = nan(length(ind),1);
            obj.theta = nan(length(ind),1);
            obj.phi = nan(length(ind),1);
            obj.prePoint = nan(length(ind),1);
            obj.postPoint = nan(length(ind),1);
            obj.alpha = nan(length(ind),1);
            
            for fi = 1 : length(ind)
                % 1. find the calculation point for kappas 
                preDist = (obj.rInMm-1) * obj.pxPerMm;
                postDist = (obj.rInMm+1) * obj.pxPerMm;
                whiskerPixLengths = zeros(length(obj.trackerData{fi}),1);
                for i = 2 : length(whiskerPixLengths)
                    whiskerPixLengths(i) = sqrt(sum((obj.trackerData{fi}(i,:) - obj.trackerData{fi}(i-1,:)).^2));
                end
                % length from the first pixel of tracker data to the intersection with the mask (i.e., "base")
                baseLength = sum(whiskerPixLengths(1:obj.baseInd(fi)));
                if ~isempty(find(cumsum(whiskerPixLengths) - baseLength >= preDist, 1, 'first')) && ~isempty(find(cumsum(whiskerPixLengths) - baseLength >= postDist, 1, 'first'))
                    obj.prePoint(fi) = find(cumsum(whiskerPixLengths) - baseLength >= preDist, 1, 'first');
                    obj.postPoint(fi) = find(cumsum(whiskerPixLengths) - baseLength >= postDist, 1, 'first');

                    if ~isempty(obj.prePoint(fi)) && ~isempty(obj.postPoint(fi)) && obj.postPoint(fi) - obj.prePoint(fi) > 3
                        x = obj.trackerData{fi}(obj.prePoint(fi):obj.postPoint(fi),1);
                        y = obj.trackerData{fi}(obj.prePoint(fi):obj.postPoint(fi),2);
                        z = obj.trackerData{fi}(obj.prePoint(fi):obj.postPoint(fi),3);

                        q = linspace(0,1, obj.postPoint(fi)-obj.prePoint(fi)+1);
                        px = polyfit(q',x,2);
                        py = polyfit(q',y,2);
                        
                        midpoint = round((obj.postPoint(fi) - obj.prePoint(fi))/2);
                        
                        % figure out the angle of tangent plane (alpha)
                        pxDot = polyder(px);
                        pyDot = polyder(py);
                        xDot = polyval(pxDot,q);
                        yDot = polyval(pyDot,q);
                        if strcmp(ws.faceSideInImage,'top') && strcmp(ws.protractionDirection,'rightward')
                            alpha = atand(xDot(midpoint) / yDot(midpoint));
                        elseif strcmp(ws.faceSideInImage,'top') && strcmp(ws.protractionDirection,'leftward')
                            alpha = -atand(xDot(midpoint) / yDot(midpoint));
                        elseif strcmp(ws.faceSideInImage,'left') && strcmp(ws.protractionDirection,'downward')
                            alpha = atand(yDot(midpoint) / xDot(midpoint));
                        elseif strcmp(ws.faceSideInImage,'left') && strcmp(ws.protractionDirection,'upward')
                            alpha = -atand(yDot(midpoint) / xDot(midpoint));
                        elseif strcmp(ws.faceSideInImage,'right') && strcmp(ws.protractionDirection,'upward')
                            alpha = atand(yDot(midpoint) / xDot(midpoint));
                        elseif strcmp(ws.faceSideInImage,'right') && strcmp(ws.protractionDirection,'downward')
                            alpha = -atand(yDot(midpoint) / xDot(midpoint));
                        elseif strcmp(ws.faceSideInImage,'bottom') && strcmp(ws.protractionDirection,'rightward')
                            alpha = -atand(xDot(midpoint) / yDot(midpoint));
                        elseif strcmp(ws.faceSideInImage,'bottom') && strcmp(ws.protractionDirection,'leftward')
                            alpha = atand(xDot(midpoint) / yDot(midpoint));
                        end
                        
                        obj.alpha(fi) = alpha;
                        
                        % Rotate the segment 
                        Rvertical = [cosd(-alpha) -sind(-alpha); sind(-alpha) cosd(-alpha)]; % rotation matrix in top view
                        refPoint = [x(midpoint); y(midpoint)];
                        projected = [(Rvertical * [x' - refPoint(1); y' - refPoint(2)] + refPoint)', z];
                        
                        % 2. horizontal & vertical kappa
                        px = polyfit(q',projected(:,1),2);
                        py = polyfit(q',projected(:,2),2);
                        pz = polyfit(q',projected(:,3),2);

                        % 2. horizontal & vertical kappa
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

                        obj.kappaH(fi) = kappasH(midpoint);
                        obj.kappaV(fi) = kappasV(midpoint);

                    end
                end
                
                % 3. theta and phi at the base                
                q = linspace(0,1);
                px = obj.fit3Data{fi}(:,1);
                py = obj.fit3Data{fi}(:,2);
                pz = obj.fit3Data{fi}(:,3);
                
                pxDot = polyder(px);                
                pyDot = polyder(py);
                pzDot = polyder(pz);
                
                xDot = polyval(pxDot,q);                
                yDot = polyval(pyDot,q);
                zDot = polyval(pzDot,q);
                
                dq = [0 diff(q)];
                
                % Arc length as a function of q, after integration below:
                R = cumsum(sqrt(xDot.^2 + yDot.^2 + zDot.^2) .* dq); % arc length segments, in pixels, times dq.
                rind = find(R >= baseLength, 1, 'first');
                
                % Angle (in degrees) as a function of q:
                % Protraction means theta is increasing.
                % Theta is 0 when perpendicular to the midline of the mouse.
                if strcmp(ws.faceSideInImage,'top') && strcmp(ws.protractionDirection,'rightward')
                    thetas = atand(xDot ./ yDot);
                    phis = atand(zDot ./ yDot);
                elseif strcmp(ws.faceSideInImage,'top') && strcmp(ws.protractionDirection,'leftward')
                    thetas = -atand(xDot ./ yDot);
                    phis = -atand(zDot ./ yDot);
                elseif strcmp(ws.faceSideInImage,'left') && strcmp(ws.protractionDirection,'downward')
                    thetas = atand(yDot ./ xDot);
                    phis = atand(yDot ./ zDot);
                elseif strcmp(ws.faceSideInImage,'left') && strcmp(ws.protractionDirection,'upward')
                    thetas = -atand(yDot ./ xDot);
                    phis = -atand(yDot ./ zDot);
                elseif strcmp(ws.faceSideInImage,'right') && strcmp(ws.protractionDirection,'upward')
                    thetas = atand(yDot ./ xDot); 
                    phis = atand(yDot ./ zDot);
                elseif strcmp(ws.faceSideInImage,'right') && strcmp(ws.protractionDirection,'downward')
                    thetas = -atand(yDot ./ xDot);
                    phis = -atand(yDot ./ zDot);
                elseif strcmp(ws.faceSideInImage,'bottom') && strcmp(ws.protractionDirection,'rightward')
                    thetas = -atand(xDot ./ yDot);
                    phis = -atand(zDot ./ yDot);
                elseif strcmp(ws.faceSideInImage,'bottom') && strcmp(ws.protractionDirection,'leftward')
                    thetas = atand(xDot ./ yDot);
                    phis = atand(zDot ./ yDot);
                else
                    error('Invalid value of property ''faceSideInImage'' or ''protractionDirection''')
                end
                obj.theta(fi) = thetas(rind);
                obj.phi(fi) = phis(rind);
            end
            
            obj.lengthAlongWhisker = obj.get_lengthAlongWhisker;
        end
        
%         function value = get_baseCoordinateTopview(obj, mask)
%             value = zeros(length(obj.trackerData),2);            
%             for i = 1 : length(obj.trackerData)
%                 xall = obj.trackerData{i}{4};
%                 xall = xall';
%                 yall = obj.trackerData{i}{5};
%                 yall = yall';
%                 whisker = [xall+1;yall+1];
%                 temp = Whisker.InterX(whisker, mask);
%                 if isempty(temp)
%                     value(i,:) = whisker(:,1)';
%                 else
%                     value(i,:) = temp(:,1);
%                 end
%             end
%         end
%         
%         function value = get_lengthAlongWhiskerTopview(obj, whiskerPoleIntersection)
%             value = nan(length(obj.baseCoordinateTopview),1);
%             ind = intersect(find(obj.baseCoordinateTopview(:,1)), find(isfinite(sum(obj.baseCoordinateTopview,2))));
%             for i = 1 : length(value)
%                 if ~isempty(whiskerPoleIntersection{i,1}) && ismember(i, ind)
%                     whisker = [obj.trackerData{i}{4}+1, obj.trackerData{i}{5}+1];
%                     dist2base = sum((whisker - obj.baseCoordinateTopview(i,:)).^2,2);
%                     baseInd = find(dist2base == min(dist2base));
%                     dist2intersect = sum((whisker - whiskerPoleIntersection{i,1}).^2,2);
%                     intersectInd = find(dist2intersect == min(dist2intersect));
%                     arcLength = [0; cumsum(sqrt(diff(whisker(:,1)).^2 + diff(whisker(:,2)).^2))];
%                     value(i) = abs(arcLength(intersectInd) - arcLength(baseInd));
%                 end
%             end
%         end
%         
%         function value = get_lengthAlongWhisker(obj)
%             value = nan(length(obj.trackerData),1);
%             inds = find(isfinite(obj.lengthAlongWhiskerTopview))';
%             for i = inds
%                 x = obj.trackerData{i}(obj.baseInd(i):end,1);
%                 y = obj.trackerData{i}(obj.baseInd(i):end,2);
%                 z = obj.trackerData{i}(obj.baseInd(i):end,3);
%                 arcLengths = cumsum((diff(x).^2 + diff(y).^2));
%                 arcLengths3d = cumsum(diff(x).^2 + diff(y).^2 + diff(z).^2);
%                 [~, ind] = min(abs(arcLengths - obj.lengthAlongWhiskerTopview(i)));
%                 value(i) = arcLengths3d(ind);
%             end
%         end
        function value = get_lengthAlongWhisker(obj)
            value = nan(length(obj.trackerData),1);
            inds = intersect(find(sum(obj.intersectPoint,2))', obj.poleUpFrames);
            if ~strcmp(obj.trialType, 'oo')
                for i = inds                
                    x = obj.trackerData{i}(:,1);
                    y = obj.trackerData{i}(:,2);
                    z = obj.trackerData{i}(:,3);
                    distToIntersection = (x-obj.intersectPoint(i,1)).^2 + (y-obj.intersectPoint(i,2)).^2 + (z-obj.intersectPoint(i,3)).^2;
                    [~, intersectInd] = min(distToIntersection);
                    x = obj.trackerData{i}(obj.baseInd(i):intersectInd, 1);
                    y = obj.trackerData{i}(obj.baseInd(i):intersectInd, 2);
                    z = obj.trackerData{i}(obj.baseInd(i):intersectInd, 3);
                    temp = cumsum(sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2));
                    value(i) = temp(end);
                end
            end
        end
        
        function show_one_3D(obj, frameNum, varargin)
%             q = linspace(0,1);
%             x = polyval(obj.fit3Data{frameNum}(:,1),q)';
%             y = polyval(obj.fit3Data{frameNum}(:,2),q)';
%             z = polyval(obj.fit3Data{frameNum}(:,3),q)';

            x = obj.trackerData{frameNum}(obj.baseInd(frameNum):end,1);
            y = obj.trackerData{frameNum}(obj.baseInd(frameNum):end,2);
            z = obj.trackerData{frameNum}(obj.baseInd(frameNum):end,3);

            plot3(x,y,z, 'k-')
            
            xt = obj.trackerData{frameNum}(obj.prePoint(frameNum):obj.postPoint(frameNum),1);
            yt = obj.trackerData{frameNum}(obj.prePoint(frameNum):obj.postPoint(frameNum),2);
            zt = obj.trackerData{frameNum}(obj.prePoint(frameNum):obj.postPoint(frameNum),3);
%             qt = linspace(0,1, obj.postPoint(frameNum)-obj.prePoint(frameNum)+1);
%             px = polyfit(qt',xt,2);
%             py = polyfit(qt',yt,2);
%             pz = polyfit(qt',zt,2);
%             plot3(polyval(px,qt), polyval(py,qt), polyval(pz,qt), 'm-');
            plot3(xt, yt, zt, 'm-');
            plot3(obj.base(frameNum,1), obj.base(frameNum,2), obj.base(frameNum,3), 'r.', 'markersize', 10)
            if nargin > 2 && varargin{1}
                whisker = [x,y,z];
                ind = find(sum((whisker-obj.intersectPoint(frameNum,:)).^2,2) == min(sum((whisker-obj.intersectPoint(frameNum,:)).^2,2)) );
                plot3(x(ind), y(ind), z(ind), 'b.', 'markersize', 10)
            end            
            axis equal
            if nargin > 3
                view(varargin{2})
            else
                view(3)
            end
            xlabel('A-P'), ylabel('M-L'), zlabel('D-V')
        end
        
        function show_all_3D(obj)
            figure, hold on,
            for i = 1 : length(obj.trackerData)
                plot3(obj.trackerData{i}(:,1), obj.trackerData{i}(:,2), obj.trackerData{i}(:,3), 'k-')
            end
            for i = 1 : length(obj.trackerData)
                if isfinite(obj.prePoint(i))
                    plot3(obj.trackerData{i}(obj.prePoint(i):obj.postPoint(i),1), obj.trackerData{i}(obj.prePoint(i):obj.postPoint(i),2), obj.trackerData{i}(obj.prePoint(i):obj.postPoint(i),3), 'm-')
                end
            end
            for i = 1 : length(obj.trackerData)
                plot3(obj.base(i,1), obj.base(i,2), obj.base(i,3), 'r.')
            end
            for i = 1 : length(obj.trackerData)
                if isfinite(obj.intersectPoint(i,1)) && obj.intersectPoint(i,1) > 0
                    plot3(obj.intersectPoint(i,1), obj.intersectPoint(i,2), obj.intersectPoint(i,3), 'b.')
                end
            end
%             plot3(obj.mask3d(:,1), obj.mask3d(:,2), obj.mask3d(:,3), '.', 'color', [0.7 0.7 0.7])
            axis equal
        end
        
        function show_all_3D_polyfit(obj)
            q = linspace(0,1);
            figure, hold on,
            for i = 1 : length(obj.trackerData)
                plot3(polyval(obj.fit3Data{i}(:,1),q), polyval(obj.fit3Data{i}(:,2),q), polyval(obj.fit3Data{i}(:,3),q), 'k-')
            end
            for i = 1 : length(obj.trackerData)
                x = obj.trackerData{i}(obj.prePoint(i):obj.postPoint(i),1);
                y = obj.trackerData{i}(obj.prePoint(i):obj.postPoint(i),2);
                z = obj.trackerData{i}(obj.prePoint(i):obj.postPoint(i),3);

                q = linspace(0,1, obj.postPoint(i)-obj.prePoint(i)+1);
                px = polyfit(q',x,2);
                py = polyfit(q',y,2);
                pz = polyfit(q',z,2);
                plot3(polyval(px,q), polyval(py,q), polyval(pz,q), 'm-');
            end
            for i = 1 : length(obj.trackerData)
                plot3(obj.base(i,1), obj.base(i,2), obj.base(i,3), 'r.')
            end                    
            axis equal
        end        
        
        function show_onebyone_3D_polyfit(obj)
            q = linspace(0,1);
            figure,
            i = 1;
            
            basex = obj.base(i,1) / obj.pxPerMm;
            basey = -obj.base(i,2) / obj.pxPerMm;
            basez = obj.base(i,3) / obj.pxPerMm;
            x = polyval(obj.fit3Data{i}(:,1),q) / obj.pxPerMm - basex;
            y = -polyval(obj.fit3Data{i}(:,2),q) / obj.pxPerMm - basey;
            z = polyval(obj.fit3Data{i}(:,3),q) / obj.pxPerMm - basez;
            plot3(x, y, z, 'k-'); hold on, 
            plot3(0, 0, 0, 'r.', 'markersize', 20); 
            hold off, 
            view(3);
%             view([-90, 0]);
            [az, el] = view;
            xl = [min(cellfun(@(x) min(x(:,1)) / obj.pxPerMm - basex, obj.trackerData)), max(cellfun(@(x) max(x(:,1)) / obj.pxPerMm - basex, obj.trackerData))];
            yl = [min(cellfun(@(x) -max(x(:,2)) / obj.pxPerMm - basey, obj.trackerData)), max(cellfun(@(x) -min(x(:,2)) / obj.pxPerMm - basey, obj.trackerData))]; 
            zl = [min(cellfun(@(x) min(x(:,3)) / obj.pxPerMm - basez, obj.trackerData)), max(cellfun(@(x) max(x(:,3)) / obj.pxPerMm - basez, obj.trackerData))];
            while( i > 0 )
                basex = obj.base(i,1) / obj.pxPerMm;
                basey = -obj.base(i,2) / obj.pxPerMm;
                basez = obj.base(i,3) / obj.pxPerMm;
                x = polyval(obj.fit3Data{i}(:,1),q) / obj.pxPerMm - basex;
                y = -polyval(obj.fit3Data{i}(:,2),q) / obj.pxPerMm - basey;
                z = polyval(obj.fit3Data{i}(:,3),q) / obj.pxPerMm - basez;
                plot3(x, y, z, 'k-'); hold on, 
                if i < 51
                    j = 1;
                else
                    j = i - 50;
                end
                for k = j : i-1
                    x = polyval(obj.fit3Data{k}(:,1),q) / obj.pxPerMm - basex;
                    y = -polyval(obj.fit3Data{k}(:,2),q) / obj.pxPerMm - basey;
                    z = polyval(obj.fit3Data{k}(:,3),q) / obj.pxPerMm - basez;
                    plot3(x, y, z, 'k-', 'color', [0.7 0.7 0.7])                    
                end
                if i > length(obj.trackerData) - 50
                    j = length(obj.trackerData);
                else
                    j = i + 50;
                end
                for k = i+1:j
                    x = polyval(obj.fit3Data{k}(:,1),q) / obj.pxPerMm - basex;
                    y = -polyval(obj.fit3Data{k}(:,2),q) / obj.pxPerMm - basey;
                    z = polyval(obj.fit3Data{k}(:,3),q) / obj.pxPerMm - basez;
                    plot3(x, y, z, 'k-', 'color', [0.7 0.7 0.7])
                end
                x = polyval(obj.fit3Data{i}(:,1),q) / obj.pxPerMm - basex;
                y = -polyval(obj.fit3Data{i}(:,2),q) / obj.pxPerMm - basey;
                z = polyval(obj.fit3Data{i}(:,3),q) / obj.pxPerMm - basez;
                plot3(x, y, z, 'k-', 'linewidth', 2); hold on, 
                plot3(0, 0, 0, 'r.', 'markersize', 25);
                hold off
                title({'Navigate using keyboard'; ['Trial # ', obj.trackerFileName]; ['Frame # ', num2str(round(obj.time(i)/obj.framePeriodInSec))]});
                xlabel('AP (mm)'), ylabel('ML (mm)'), zlabel('DV (mm)')                
%                 set(gca, 'linewidth', 3, 'fontweight', 'bold', 'fontsize', 15)
                axis equal
                view(az,el), xlim(xl), ylim(yl), zlim(zl);
                [i, az, el] = obj.keyboard_navigation_3d(i, length(obj.trackerData));
            end
        end
        
        function show_onebyone_3D_polyfit_with_touch(obj)
            q = linspace(0,1);
            figure,
            i = 1;
            
            basex = obj.base(i,1) / obj.pxPerMm;
            basey = -obj.base(i,2) / obj.pxPerMm;
            basez = obj.base(i,3) / obj.pxPerMm;
            x = polyval(obj.fit3Data{i}(:,1),q) / obj.pxPerMm - basex;
            y = -polyval(obj.fit3Data{i}(:,2),q) / obj.pxPerMm - basey;
            z = polyval(obj.fit3Data{i}(:,3),q) / obj.pxPerMm - basez;
            plot3(x, y, z, 'color'); hold on, 
            plot3(0, 0, 0, 'r.', 'markersize', 20); 
            hold off, 
            view(3);
%             view([-90, 0]);
            [az, el] = view;
            xl = [min(cellfun(@(x) min(x(:,1)) / obj.pxPerMm - basex, obj.trackerData)), max(cellfun(@(x) max(x(:,1)) / obj.pxPerMm - basex, obj.trackerData))];
            yl = [min(cellfun(@(x) -max(x(:,2)) / obj.pxPerMm - basey, obj.trackerData)), max(cellfun(@(x) -min(x(:,2)) / obj.pxPerMm - basey, obj.trackerData))]; 
            zl = [min(cellfun(@(x) min(x(:,3)) / obj.pxPerMm - basez, obj.trackerData)), max(cellfun(@(x) max(x(:,3)) / obj.pxPerMm - basez, obj.trackerData))];
            while( i > 0 )
                basex = obj.base(i,1) / obj.pxPerMm;
                basey = -obj.base(i,2) / obj.pxPerMm;
                basez = obj.base(i,3) / obj.pxPerMm;
                x = polyval(obj.fit3Data{i}(:,1),q) / obj.pxPerMm - basex;
                y = -polyval(obj.fit3Data{i}(:,2),q) / obj.pxPerMm - basey;
                z = polyval(obj.fit3Data{i}(:,3),q) / obj.pxPerMm - basez;
                plot3(x, y, z, 'k-'); hold on, 
                if i < 51
                    j = 1;
                else
                    j = i - 50;
                end
                for k = j : i-1
                    x = polyval(obj.fit3Data{k}(:,1),q) / obj.pxPerMm - basex;
                    y = -polyval(obj.fit3Data{k}(:,2),q) / obj.pxPerMm - basey;
                    z = polyval(obj.fit3Data{k}(:,3),q) / obj.pxPerMm - basez;
                    plot3(x, y, z, 'k-', 'color', [0.7 0.7 0.7])                    
                end
                if i > length(obj.trackerData) - 50
                    j = length(obj.trackerData);
                else
                    j = i + 50;
                end
                for k = i+1:j
                    x = polyval(obj.fit3Data{k}(:,1),q) / obj.pxPerMm - basex;
                    y = -polyval(obj.fit3Data{k}(:,2),q) / obj.pxPerMm - basey;
                    z = polyval(obj.fit3Data{k}(:,3),q) / obj.pxPerMm - basez;
                    plot3(x, y, z, 'k-', 'color', [0.7 0.7 0.7])
                end
                x = polyval(obj.fit3Data{i}(:,1),q) / obj.pxPerMm - basex;
                y = -polyval(obj.fit3Data{i}(:,2),q) / obj.pxPerMm - basey;
                z = polyval(obj.fit3Data{i}(:,3),q) / obj.pxPerMm - basez;
                plot3(x, y, z, 'k-', 'linewidth', 2); hold on, 
                plot3(0, 0, 0, 'r.', 'markersize', 25);
                hold off
                title({'Navigate using keyboard'; ['Trial # ', obj.trackerFileName]; ['Frame # ', num2str(round(obj.time(i)/obj.framePeriodInSec))]});
                xlabel('AP (mm)'), ylabel('ML (mm)'), zlabel('DV (mm)')                
%                 set(gca, 'linewidth', 3, 'fontweight', 'bold', 'fontsize', 15)
                axis equal
                view(az,el), xlim(xl), ylim(yl), zlim(zl);
                [i, az, el] = obj.keyboard_navigation_3d(i, length(obj.trackerData));
            end
        end
        
        function show_onebyone_3D_print(obj)
            
            
            
%             viewVal = [30, -30];
            q = linspace(0,1);
            figure,
            subplot(221)            
            i = 1;
            
            basePoint = obj.base(i,:);
            
            x = polyval(obj.fit3Data{i}(:,1),q);
            y = polyval(obj.fit3Data{i}(:,2),q);
            z = polyval(obj.fit3Data{i}(:,3),q);
            
            plot3(x, y, z, 'k-')
            ax1 = gca;
            hold(ax1,'on')
            plot3(ax1, (obj.base(i,1) - basePoint(1)) / obj.pxPerMm, -(obj.base(i,2) - basePoint(2)) / obj.pxPerMm, (obj.base(i,3) - basePoint(3)) / obj.pxPerMm, 'r.', 'markersize', 20); 
            hold(ax1,'off')
            
%             view(ax1, viewVal)

            subplot(223)            
            i = 1;
            x = polyval(obj.fit3Data{i}(:,1),q);
            y = polyval(obj.fit3Data{i}(:,2),q);
            z = polyval(obj.fit3Data{i}(:,3),q);
            plot3(x, y, z, 'k-'); 
            ax2 = gca;
            hold(ax2,'on')
            plot3(ax2, (obj.base(i,1) - basePoint(1)) / obj.pxPerMm, -(obj.base(i,2) - basePoint(2)) / obj.pxPerMm, (obj.base(i,3) - basePoint(3)) / obj.pxPerMm, 'r.', 'markersize', 20);
            hold(ax2,'off')       
            view(ax2, [-90, 0])
            
            subplot(222)            
            i = 1;
            x = polyval(obj.fit3Data{i}(:,1),q);
            y = polyval(obj.fit3Data{i}(:,2),q);
            z = polyval(obj.fit3Data{i}(:,3),q);
            plot3(x, y, z, 'k-')
            ax3 = gca;
            hold(ax3,'on')
            plot3(ax3, (obj.base(i,1) - basePoint(1)) / obj.pxPerMm, -(obj.base(i,2) - basePoint(2)) / obj.pxPerMm, (obj.base(i,3) - basePoint(3)) / obj.pxPerMm, 'r.', 'markersize', 20);
            hold(ax3,'off')            
            view(ax3, [0, -90])
            
            while( i > 0 )                
%                 plot3(ax1, (polyval(obj.fit3Data{i}(:,1),q) - basePoint(1)) / obj.pxPerMm, -(polyval(obj.fit3Data{i}(:,2),q) - basePoint(2)) / obj.pxPerMm, ...
%                     (polyval(obj.fit3Data{i}(:,3),q) - basePoint(3)) / obj.pxPerMm, 'k-', 'linewidth', 2)
                [az, el] = view(ax1);
                basePoint = obj.base(i,:);
                
                x = (polyval(obj.fit3Data{i}(:,1),q) - basePoint(1)) / obj.pxPerMm;
                y = -(polyval(obj.fit3Data{i}(:,2),q) - basePoint(2)) / obj.pxPerMm;
                z = (polyval(obj.fit3Data{i}(:,3),q) - basePoint(3)) / obj.pxPerMm;
                
                xlimVal = [min(x)-1 max(x)+1];
                ylimVal = [min(y)-1 max(y)+1];
                zlimVal = [min(z)-1 max(z)+1];
                
                basex = (obj.base(i,1) - basePoint(1)) / obj.pxPerMm;
                basey = -(obj.base(i,2) - basePoint(2)) / obj.pxPerMm;
                basez = (obj.base(i,3) - basePoint(3)) / obj.pxPerMm;
                
                plot3(ax1, x, y, z, 'k-', 'linewidth', 2);
                hold(ax1,'on')
                plot3(ax1, basex, basey, basez, 'r.', 'markersize', 25);
                minInd = zeros(2,1);
                [~, minInd(1)] = min( abs(  sqrt( (x-basex).^2 + (y-basey).^2 + (z-basez).^2 )  -  (obj.rInMm-1)) );
                [~, minInd(2)] = min( abs(  sqrt( (x-basex).^2 + (y-basey).^2 + (z-basez).^2 )  -  (obj.rInMm+1)) );                
                plot3(ax1, x(minInd(1):minInd(2)), y(minInd(1):minInd(2)), z(minInd(1):minInd(2)), 'm-', 'linewidth', 2);
                hold(ax1,'off')
                title(ax1, {['Trial # ', obj.trackerFileName]; ['Frame # ', num2str(round(obj.time(i)/obj.framePeriodInSec))]});
                xlabel(ax1, 'AP (mm)'), ylabel(ax1, 'ML (mm)'), zlabel(ax1, 'DV (mm)')
                axis(ax1, 'equal')
                xlim(ax1, xlimVal), ylim(ax1, ylimVal), zlim(ax1, zlimVal)
%                 view(ax1, viewVal)
                view(ax1, [az, el])
                
                plot3(ax2, x, y, z, 'k-', 'linewidth', 2)
                hold(ax2,'on')                
                plot3(ax2, basex, basey, basez, 'r.', 'markersize', 25);
                plot3(ax2, x(minInd(1):minInd(2)), y(minInd(1):minInd(2)), z(minInd(1):minInd(2)), 'm-', 'linewidth', 2);
                hold(ax2,'off')
                title(ax2, 'Front view');
                xlabel(ax2, 'AP (mm)'), ylabel(ax2, 'ML (mm)'), zlabel(ax2, 'DV (mm)')
                axis(ax2, 'equal')
                xlim(ax2, xlimVal), ylim(ax2, ylimVal), zlim(ax2, zlimVal)
                view(ax2, [-90, 0])

                plot3(ax3, x, y, z, 'k-', 'linewidth', 2)
                hold(ax3,'on')                
                plot3(ax3, basex, basey, basez, 'r.', 'markersize', 25);
                plot3(ax3, x(minInd(1):minInd(2)), y(minInd(1):minInd(2)), z(minInd(1):minInd(2)), 'm-', 'linewidth', 2);
                hold(ax3,'off')
                title(ax3, 'Top-down view');
                xlabel(ax3, 'AP (mm)'), ylabel(ax3, 'ML (mm)'), zlabel(ax3, 'DV (mm)')
                axis(ax3, 'equal')
                xlim(ax3, xlimVal), ylim(ax3, ylimVal), zlim(ax3, zlimVal)
                view(ax3, [0, -90])
                
                i = obj.keyboard_navigation_3d(i, length(obj.trackerData));
            end
        end        
        
        function show_onebyone_3D_video(obj)
            viewVal = [-30, 30];
            q = linspace(0,1);
            figure,
            subplot(221)            
            i = 1;
            x = polyval(obj.fit3Data{i}(:,1),q);
            y = polyval(obj.fit3Data{i}(:,2),q);
            z = polyval(obj.fit3Data{i}(:,3),q);
            xl = [min(cellfun(@(x) min(x(:,1)), obj.trackerData)), max(cellfun(@(x) max(x(:,1)), obj.trackerData))];
            yl = [min(cellfun(@(x) min(x(:,2)), obj.trackerData)), max(cellfun(@(x) max(x(:,2)), obj.trackerData))]; 
            zl = [min(cellfun(@(x) min(x(:,3)), obj.trackerData)), max(cellfun(@(x) max(x(:,3)), obj.trackerData))];
            
            plot3(x, y, z, 'k-')
            ax1 = gca;
            hold(ax1,'on')
            plot3(ax1, obj.base(i,1), obj.base(i,2), obj.base(i,3), 'r.', 'markersize', 20); 
            hold(ax1,'off')
            view(ax1, viewVal),  xlim(ax1, xl), ylim(ax1, yl), zlim(ax1, zl);

            subplot(223)            
            i = 1;
            x = polyval(obj.fit3Data{i}(:,1),q);
            y = polyval(obj.fit3Data{i}(:,2),q);
            z = polyval(obj.fit3Data{i}(:,3),q);
            plot3(x, y, z, 'k-'); 
            ax2 = gca;
            hold(ax2,'on')
            plot3(ax2, obj.base(i,1), obj.base(i,2), obj.base(i,3), 'r.', 'markersize', 20); 
            hold(ax2,'off')       
            view(ax2, [-90, 0]),  xlim(ax2, xl), ylim(ax2, yl), zlim(ax2, zl);
            
            
            subplot(222)            
            i = 1;
            x = polyval(obj.fit3Data{i}(:,1),q);
            y = polyval(obj.fit3Data{i}(:,2),q);
            z = polyval(obj.fit3Data{i}(:,3),q);
            plot3(x, y, z, 'k-')
            ax3 = gca;
            hold(ax3,'on')
            plot3(ax3, obj.base(i,1), obj.base(i,2), obj.base(i,3), 'r.', 'markersize', 20); 
            hold(ax1,'off')            
            view(ax3, [0, 90]),  xlim(ax3, xl), ylim(ax3, yl), zlim(ax3, zl);
            


            while( i > 0 )                
                plot3(ax1, polyval(obj.fit3Data{i}(:,1),q), polyval(obj.fit3Data{i}(:,2),q), polyval(obj.fit3Data{i}(:,3),q), 'k-', 'linewidth', 2)
                x = polyval(obj.fit3Data{i}(:,1),q);
                y = polyval(obj.fit3Data{i}(:,2),q);
                z = polyval(obj.fit3Data{i}(:,3),q);
                basex = obj.base(i,1);
                basey = obj.base(i,2);
                basez = obj.base(i,3);
                plot3(ax1, x, y, z, 'k-', 'linewidth', 2);
                hold(ax1,'on')
                plot3(ax1, obj.base(i,1), obj.base(i,2), obj.base(i,3), 'r.', 'markersize', 25);
                minInd = zeros(2,1);
                [~, minInd(1)] = min( abs(  sqrt( (x-basex).^2 + (y-basey).^2 + (z-basez).^2 )  -  (obj.rInMm-1) * obj.pxPerMm  ) );
                [~, minInd(2)] = min( abs(  sqrt( (x-basex).^2 + (y-basey).^2 + (z-basez).^2 )  -  (obj.rInMm+1) * obj.pxPerMm  ) );                
                plot3(ax1, x(minInd(1):minInd(2)), y(minInd(1):minInd(2)), z(minInd(1):minInd(2)), 'm-', 'linewidth', 2);
                hold(ax1,'off')
                title(ax1, {['Trial # ', obj.trackerFileName]; ['Frame # ', num2str(round(obj.time(i)/obj.framePeriodInSec))]});
                xlabel(ax1, 'Rostro-caudal'), ylabel(ax1, 'Medio-lateral'), zlabel(ax1, 'Dorso-ventral')
                axis(ax1, 'equal')
                view(ax1, viewVal), xlim(ax1, xl), ylim(ax1, yl), zlim(ax1, zl);
                
                plot3(ax2, polyval(obj.fit3Data{i}(:,1),q), polyval(obj.fit3Data{i}(:,2),q), polyval(obj.fit3Data{i}(:,3),q), 'k-', 'linewidth', 2)
                hold(ax2,'on')
                plot3(ax2, x, y, z, 'k-', 'linewidth', 2)
                plot3(ax2, obj.base(i,1), obj.base(i,2), obj.base(i,3), 'r.', 'markersize', 25);
                plot3(ax2, x(minInd(1):minInd(2)), y(minInd(1):minInd(2)), z(minInd(1):minInd(2)), 'm-', 'linewidth', 2);
                hold(ax2,'off')
                title(ax2, 'Front view');
                xlabel(ax2, 'Rostro-caudal'), ylabel(ax2, 'Medio-lateral'), zlabel(ax2, 'Dorso-ventral')
                axis(ax2, 'equal')
                view(ax2, [-90, 0]), xlim(ax2, xl), ylim(ax2, yl), zlim(ax2, zl);

                plot3(ax3, polyval(obj.fit3Data{i}(:,1),q), polyval(obj.fit3Data{i}(:,2),q), polyval(obj.fit3Data{i}(:,3),q), 'k-', 'linewidth', 2)
                hold(ax3,'on')
                plot3(ax3, x, y, z, 'k-', 'linewidth', 2)
                plot3(ax3, obj.base(i,1), obj.base(i,2), obj.base(i,3), 'r.', 'markersize', 25);
                plot3(ax3, x(minInd(1):minInd(2)), y(minInd(1):minInd(2)), z(minInd(1):minInd(2)), 'm-', 'linewidth', 2);
                hold(ax3,'off')
                title(ax3, 'Top-down view');
                xlabel(ax3, 'Rostro-caudal'), ylabel(ax3, 'Medio-lateral'), zlabel(ax3, 'Dorso-ventral')
                axis(ax3, 'equal')
                view(ax3, [0, 90]), xlim(ax3, xl), ylim(ax3, yl), zlim(ax3, zl);
                
                i = obj.keyboard_navigation_3d(i, length(obj.trackerData));
            end
        end
        
        function show_onebyone_3D(obj)
            figure, 
            i = 1; plot3(obj.trackerData{i}(:,1), obj.trackerData{i}(:,2), obj.trackerData{i}(:,3), 'k-'); hold on, 
            plot3(obj.base(i,1), obj.base(i,2), obj.base(i,3), 'r.', 'markersize', 20); 
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
                plot3(obj.base(i,1), obj.base(i,2), obj.base(i,3), 'r.', 'markersize', 25),
                hold off
                title({'Navigate using keyboard'; ['Trial # ', obj.trackerFileName]; ['Frame # ', num2str(round(obj.time(i)/obj.framePeriodInSec))]});
                xlabel('Rostro-caudal'), ylabel('Medio-lateral'), zlabel('Dorso-ventral')                
                set(gca, 'fontsize', 14)
                axis equal
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
    methods

    end
end
        