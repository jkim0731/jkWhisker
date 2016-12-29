classdef WhiskerSignalTrial_NX < Whisker.WhiskerSignalTrial
    %
    %
    %   WhiskerSignalTrial_NX < Whisker.WhiskerSignalTrial
    %
    %
    %
    %
    properties
        bar_pos_trial = [];
        bar_time_win = [];
        SoloTrial = [];
        contacts = {};
        mKappaNearBar = {};
        mThetaNearBar = {};
        imagePixelDimsXY = [];
        distToBar = {};
        
%         touch_on = {};
%         touch_off = {};
        
    end
    
    properties (Dependent = true)
        contact_params % physical parameters around contact times.
        whisk_amplitude
%         dToBar % distance from whisker to bar, after trajectory extrapolation
        deltaKappa
    end
    
    methods (Access = public)
        function obj = WhiskerSignalTrial_NX(w, varargin)
            %
            if nargin==0
                ws_args = {};
            end
            
            if isa(w, 'Whisker.WhiskerTrial')
                ws_args = {w, varargin{1}};
            else
                ws_args = {};
            end
            obj = obj@Whisker.WhiskerSignalTrial(ws_args{:});
            
            if isa(w, 'Whisker.WhiskerSignalTrial')
                ws_prop = properties(w);
                for i=1:length(ws_prop),
                    obj.(ws_prop{i}) = w.(ws_prop{i});
                end
            end
            ntraj = length(obj.trajectoryIDs);
%             obj.deltaKappa = cell(1,ntraj);
            obj.mKappaNearBar = cell(1,ntraj);
            obj.mThetaNearBar = cell(1,ntraj);
            obj.contacts = cell(1,ntraj);
        end
        
        function bar_pos_trial = get_bar_pos_trial(obj,bar_pos_trial)
            if nargin < 2
                bar_pos_trial = uigetfile('*.mat','Select bar coordinates file');
            elseif ischar(bar_pos_trial)
                x = load(bar_pos_trial);
                xx = fieldnames(x);
                bar_pos_trial = x.(xx)(obj.trialNum);
            end
            
            obj.bar_pos_trial = bar_pos_trial;
        end
        
        function tip_coords = get_tip_coords(obj,tid,varargin)
            % varargin{1}, restrictTime, [startTimeInSec endTimeInSec]
            t = obj.get_time(tid);
            if ~isempty(varargin) & ~isempty(varargin{1})
                restrictTime = varargin{1};
            else
                restrictTime = [min(t) max(t)];
            end
            
            tip_coords = nan(length(t),2);
            
            frameInds = find(t >= restrictTime(1) & t <= restrictTime(2));
            
            ind = find(obj.trajectoryIDs==tid);
            
            fittedX = obj.polyFits{ind}{1};
            fittedY = obj.polyFits{ind}{2};
            
            nframes = length(frameInds);
            %             x = cell(1,nframes);
            %             y = cell(1,nframes);
            %
            q = linspace(0,1);
            
            for k=1:nframes
                f = frameInds(k);
                
                px = fittedX(f,:);
                py = fittedY(f,:);
                x = polyval(px,q);
                y = polyval(py,q);
                tip_coords(f,1) = x(end);
                tip_coords(f,2) = y(end);
            end
            
        end
        
        function dist = get_tip_bar_distance(obj, tid, varargin)
            % distance from traced whisker tip to bar center, no
            % extrapolation. 
            if ~isempty(obj.bar_pos_trial)
                bc = obj.bar_pos_trial;
            else
                bc = obj.barPosClean;
            end
                t = obj.get_time(tid);
            if ~isempty(varargin) & ~isempty(varargin{1})
                restrictTime = varargin{1};
            else
                restrictTime = [min(t) max(t)];
            end
            
            dist = nan(length(t),1);
            
            frameInds = find(t >= restrictTime(1) & t <= restrictTime(2));
            
            tip_coords = obj.get_tip_coords(tid,restrictTime);
            
            for i = 1:length(frameInds)
                dist(frameInds(i)) = norm(tip_coords(frameInds(i),:) - bc(i,[2 3]));
            end
        end
        
        function ang_diff = get_whisker_bar_angular_difference(obj,tid,varargin)
            wid = find(obj.trajectoryIDs == tid);
            if ~isempty(obj.bar_pos_trial)
                bar_coord = obj.bar_pos_trial;
            else
                bar_coord = obj.barPos(1000,[2 3]);
            end
            ang_bar = atand((obj.follicleCoordsX{wid}-bar_coord(1))/abs(bar_coord(2)-obj.follicleCoordsY{wid}));
            % if ang_diff < 0, whisker is at more protracted postion than
            % bar location. 
            ang_diff =  obj.theta{wid} - ang_bar;
        end
        
        function contacts = get_contacts(obj,tid,varargin)
            % Contact time detection based on whisker tip to bar distance
            % and curvature change
            % INPUT: varargin{1}, the time window where bar is presented.
            %        Does not need to be precise, e.g., [1.01 3]
            if isempty(varargin)
                bar_time = obj.bar_time_win;
            else
                bar_time = varargin{1};
            end
            
            if isempty(bar_time)
                bar_time = input('Input bar time window [start end]: ');
            end
            w_ind = find(obj.trajectoryIDs == tid);
            bar_frames = find(obj.time{w_ind}>bar_time(1) & obj.time{w_ind}<bar_time(2));
            
            %              d0 = obj.get_tip_bar_distance(tid);
            
            contacts = whisker_touch_detection2(obj, tid, bar_frames);
            obj.contacts{w_ind} = contacts;
        end
        
        function [mTheta, mKappa, rROI] = mean_theta_kappa_near_bar(obj,tid, len_avg, proximity_threshold)
            % based on
            % [rNearest,thetaNearest,kappaNearest,yNearest,xNearest,dist,t] = get_r_theta_kappa_nearest_bar(obj,tid,proximity_threshold)
            %
            %
            % INPUTS:
            %       tid:  Trajectory ID (as an integer) starting with 0.
            %       len_avg:  the length (in mm) on the whisker near the bar over
            %               which to compute the mean kappa and theta.
            %       proximity_threshold:  Optional argument giving distance from nearest point on whisker
            %               to bar center, in units of bar radius, beyond which
            %               the whisker will be extrapolated along the last theta in
            %               order to determine distance between whisker and bar.
            %
            % OUTPUTS:
            %       rROI: Arc-length (radial) distance of points in the region near center of the bar. Units of pixels.
            %       mTheta: mean Theta of points in the region near the bar.
            %       mKappa: mean Kappa at the points within the region near the bar.
            %
            if nargin < 3,
                len_avg = 2; % default, 2 mm.
            end
            nPix = round(len_avg * obj.pxPerMm); % convert to number of pixels
            w_ind = find(obj.trajectoryIDs == tid);
            t = obj.time{w_ind};
            f = t / obj.framePeriodInSec;
            nframes = length(f);
%             [rNearest,thetaNearest,kappaNearest,yNearest,xNearest,dist,t] = obj.get_r_theta_kappa_nearest_bar(tid);
            [R,THETA,KAPPA] = obj.arc_length_theta_and_kappa(tid);
            
            mTheta = zeros(1,nframes);
            mKappa = zeros(1,nframes);
            rROI = cell(1,nframes);
            for k=1:nframes
                inds = find(R{k} > R{k}(end)-nPix);
                rROI{k} = R{k}(inds);
                mTheta(k) = mean(THETA{k}(inds));
                mKappa(k) = mean(KAPPA{k}(inds));
            end
            obj.mKappaNearBar{w_ind} = mKappa;
            obj.mThetaNearBar{w_ind} = mTheta;
            
        end
        
        function out = get_position_extrema(obj, tid, varargin)
            w_ind = find(obj.trajectoryIDs == tid);
            th = obj.theta{w_ind};
            th_medfilt = medfilt1(th,3);
            if size(th_medfilt,1)==1
                th_medfilt = th_medfilt';
            end
            [lmaxval, lmaxind] = lmax_pw(th_medfilt, 3);
            [lminval, lminind] = lmin_pw(th_medfilt, 3);
            out.lmaxind = lmaxind;
            out.lminind = lminind;
        end
        
        function dToBar = get_distToBar(obj,tid)
            % Compute the shortest distance from whisker to bar center
            % Extrapolate whisker poly fit lines if necessary.
            
            dToBar = cell(1,length(obj.trajectoryIDs));
            if nargin<2
                tid = obj.trajectoryIDs;
            end
            for i = 1:length(tid);
                wNo = find(obj.trajectoryIDs == tid(i));
                t = obj.time{wNo};
                f = t / obj.framePeriodInSec;
                nframes = length(f);
                q = linspace(0,1); % This must be same as that used in obj.arc_length_theta_and_kappa().
                
                fittedX = obj.polyFits{wNo}{1};
                fittedY = obj.polyFits{wNo}{2};
                
                for k=1:nframes
                    %                 disp(['Computing nearest-point values for frame=' int2str(k)])
                    px = fittedX(k,:);
                    py = fittedY(k,:);
                    
                    x = polyval(px,q);
                    y = polyval(py,q);
                    
                    xb = obj.bar_pos_trial(1);
                    yb = obj.bar_pos_trial(2);
                    xf = obj.follicleCoordsX{wNo}(k); % follicle coordinates
                    yf = obj.follicleCoordsY{wNo}(k);
                    Lb = norm([xb-xf yb-yf]);
                    Ltip = norm([x(end)-xf  y(end)-yf]);
                    
                    deltaL = abs(Lb - Ltip); 
                    % extrapolate poly fit if whisker length is smaller
                    % than bar distance to follicle.
                    x_ex = []; y_ex=  [];
                    if Lb > Ltip
                        q_ex = [q linspace(1, 1+deltaL/Ltip, round(deltaL/Ltip))]; 
                        x_ex = polyval(px,q_ex);
                        y_ex = polyval(py,q_ex);
%                         
%                         I = mean(diff(y));
%                         yi = by-obj.barRadius*2: abs(I): min(y);
%                         tic;
%                         xi = interp1(y,x,yi,'linear','extrap');
%                         if xi(end)>obj.imagePixelDimsXY(1) || xi(end)<=1
%                             xi = []; yi = [];
%                         end
%                         tt = toc; 
%                         if tt>5, error('Extrapolating time out'); end
%                         x = [x xi(2:end)];
%                         y = [y yi(2:end)];
                    end
                    d = [];
                    for j=1:length(x)
                        d(j) = norm([x(j);y(j)] - [xb;yb]); % sqrt((x(j)-bx)^2 + (y(j)-by)^2);
                    end
                    dToBar{wNo}(k) = min(d);
                end
                obj.distToBar{wNo} = dToBar{wNo};
            end
        end 
        
        function plot_fitted_whisker_frameInds(obj, tid, frameInds, varargin)
            %
            %   plot_fitted_whisker_time_projection(obj, tid, varargin)
            %
            % tid: trajectory ID.  Only a single ID is allowed.
            %
            % varargin{1}: Plot color/symbol string.
            %
            % varargin{2}: Optional 2 x 1 vector giving starting and ending times (in seconds) to include in
            %               returned image, inclusive, starting with 0.
            %               Of format: [startTimeInSec endTimeInSec]. Can be empty ([]) to allow
            %               access to varargin{3}.
            %
            % varargin{3}: Cell array with two elements. First element is a 1x2 vector giving radial distance
            %              beginning and end points (inclusive) to draw in the color given in
            %              varargin{1}, in format [startRadialDistance stopRadialDistance]. Radial distance
            %              is whisker arc length moving outward from follicle and is in units of pixels.
            %              The second element of the cell array is the plot color/symbol string to use when
            %              plotting the radial segment of the whisker between startRadialDistance and
            %              stopRadialDistance.
            %
            %
            %
            if numel(tid) > 1
                error('Only a single trajectory ID is alowed.')
            end
            
            if nargin > 3
                plotString = varargin{1};
            else
                plotString = 'k-';
            end
            if nargin > 4
                if ~iscell(varargin{2})
                    error('Wrong format for varargin{3}.')
                end
                plotString2 = varargin{3}{2};
                rstart = varargin{3}{1}(1);
                rstop = varargin{3}{1}(2);
            end

            hold on; axis ij
            
            ind = obj.trajectoryIDs==tid;
            if max(ind) < 1
                error('Could not find specified trajectory ID.')
            end
            
            if isempty(obj.polyFits)
                error('obj.polyFits is empty.')
            end
            
            fittedX = obj.polyFits{ind}{1};
            fittedY = obj.polyFits{ind}{2};
            
            nframes = length(frameInds);
            x = cell(1,nframes);
            y = cell(1,nframes);
            
            q = linspace(0,1);
            
            for k=1:nframes
                f = frameInds(k);
                
                px = fittedX(f,:);
                py = fittedY(f,:);
                x{k} = polyval(px,q);
                y{k} = polyval(py,q);
            end
            
            if nargin > 4
                [R,THETA,KAPPA] = obj.arc_length_theta_and_kappa(tid);
                for k=1:nframes
                    xx = x{k};
                    yy = y{k};
                    f = frameInds(k);
                    radDist = R{f};
                    plot(xx,yy,plotString)
                    rind = radDist >= rstart & radDist <= rstop;
                    plot(xx(rind),yy(rind),plotString2)
                end
            else
                if length(frameInds) > 1
                    fr_color = color_mapping_nx(frameInds,[],jet);
                else
                    fr_color = 'r';
                end
                for k=1:nframes
                    plot(x{k},y{k},plotString,'Color',fr_color(k,:))
                end
            end
        end
        
    end
    
    methods % Dependent property methods; cannot have attributes.
        
        function value = get.contact_params(obj)
            nTraj = length(obj.trajectoryIDs);
            value = cell(1,nTraj);
            for i = 1 : nTraj % whiskers
                if i > length(obj.contacts)
                    continue
                else
                    tid = obj.trajectoryIDs(i);
                    deltaKappa = obj.deltaKappa{i};
                    theta = obj.theta{i};
                    theta_LM = obj.get_position_extrema(tid);
                    time_stamp = obj.get_time(tid);
                    vel = obj.get_velocity(tid);
                    
                    for ii = 1:length(obj.contacts{i}) % loop through touch windows for the current whisker
                        inds = obj.contacts{i}{ii}; % index of the first contact frame
                        value{i}{ii}.trialNo = obj.trialNum;
                        value{i}{ii}.frameInds = inds;
                        value{i}{ii}.ts = time_stamp(inds);
                        value{i}{ii}.vel = nanmean(vel(inds(1)-4 : inds(1))); % pre-contact velocity
                        kappaCont = deltaKappa(inds);
                        [maxDeltaKappa, I] = max(abs(kappaCont));
                        value{i}{ii}.maxDeltaKappa = maxDeltaKappa; % max curvature change in this touch window
                        value{i}{ii}.curv_dot = [0 diff(kappaCont)]./[0 diff(time_stamp(inds))];
                        lm_inds = sort([theta_LM.lmaxind theta_LM.lminind]);
                        lm_ind_0 = lm_inds(find(lm_inds < inds(1), 1, 'last'));
                        value{i}{ii}.theta_set = theta(lm_ind_0);
                        value{i}{ii}.theta1 = theta(inds(1)); % whisker position upon touch
                        value{i}{ii}.theta_slope = (value{i}{ii}.theta_set - value{i}{ii}.theta1)/(time_stamp(inds(1))-time_stamp(lm_ind_0));
                        
                    end
                end
            end
        end
        
        function value = get.whisk_amplitude(obj)
            value = cell(1,length(obj.trajectoryIDs));
            for i = 1:length(value)
                if i > length(obj.theta)
                    continue
                else
                    th = obj.theta{i};
                    [lmn, indmin] = lmin_pw(th',15);
                    [lmx, indmax] = lmax_pw(th',15);
                    xi = 1:length(th);
                    indmin = unique(indmin);
                    lmn = th(indmin);
                    indmax = unique(indmax);
                    lmx = th(indmax);
                    lmni = interp1(indmin,lmn,xi);
                    lmxi = interp1(indmax,lmx,xi);
                    value{i} = abs(lmni)-abs(lmxi);
                end
            end
        end
        
        function deltaKappa = get.deltaKappa(obj)
            % deltaKappa, is the kappa subtracted by the mean over a
            % baseline period define by a low velocity pre-bar period.
            % kappa is averaged over an arc length roi for each whiskers
            tid = obj.trajectoryIDs;
            
            if ~isempty(obj.bar_time_win)
                baseline_window = [0 obj.bar_time_win(1)-0.4];
            else
                baseline_window = [0 0.5]; % ad hoc, alwasy >0.5 sec deley before bar moves.
            end
            t = obj.get_time(tid);
            bs_win_ind = find(t< baseline_window(2) & t > baseline_window(1));
            
            w_ind = find(obj.trajectoryIDs == tid);
%             [R,theta,kappa] = obj.arc_length_theta_and_kappa(tid);
            
            % get one kappa value (averaged from an roi segment) for each frame 
%             for k = 1:length(R) % frames
%                 ind1 = find(R{k}< roi(2) & R{k}>roi(1));
%                 mkappa = mean(kappa{k}(ind1));
%             end
%             [t,theta,kappa] = obj.mean_theta_and_kappa(tid,radial_window_theta,radial_window_kappa)  
            % define the baseline
            % find low velocity period
            for i = 1:length(tid)
                kappa = obj.kappa{w_ind(i)};
                vel = obj.get_velocity(tid(i));
                ind = find(abs(vel)<prctile(abs(vel),20));
            % baseline kappa is averaged across frames with low velocity,
            % but within defined time window, e.g., pre-bar time.
                bs_kappa = nanmean(kappa(ind(ismember(ind,bs_win_ind))));
                deltaKappa{w_ind(i)} = (kappa - bs_kappa)*1000; % in 1/meter
            end
            
        end
        
            
    end
end
