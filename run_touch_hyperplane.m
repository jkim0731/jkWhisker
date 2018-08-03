function run_touch_hyperplane(mouseName, sessionName, bSession, whiskerBaseDir, ppm, varargin)
% For 2pad (2-port angle distance), WT, WST, and WL is generated during prebuild_JK.m

%% Setup whisker array builder 
useGPU = 0;
options.WindowStyle = 'normal';
optionsFin.WindowStyle = 'normal';
optionsFin.Units = 'normalized';
optionsFin.Position = [0.3 0.3 0.2 0.2];

try
    whisker_d = [whiskerBaseDir mouseName sessionName '\'];
    cd(whisker_d)
    wstFnList = dir('*_WST.mat');
    wstNums = zeros(length(wstFnList),1);
    for wi = 1 : length(wstFnList)
        wstNums(wi) = str2double(wstFnList(wi).name(1:end-8));
    end
catch
    error('No WST files')
end    

if nargin == 6
    skip = varargin{1};
else
    skip = 'skip';
end

if strcmp(skip,'skip') && ~isempty(dir([whisker_d, '*touch_hp.mat']))
    disp('Hyperplane already exists. Exit.')
    return
end

if strcmp(sessionName,'S99')
    sessionName = 'S17';
elseif strcmp(sessionName,'S91')
    sessionName = 'S01';
end

filelist=dir([whisker_d '*.measurements']);

dirTrialNums=zeros(1,size(filelist,1));
% trialNums=[];  % enter which trial nums to process 

% %%
% Assign the trial numbers to existing .measurements files in the directory
% NOTE : This assumes that the .measurements files have leading numbers
% corresponding to trial number in string positions 1:end-13 of the file
% name. These index numbers may need to be changed to match up to the
% numerical code of the trial number.  (2016/09/05 JK)

for i=1:length(filelist)
    dirTrialNums(i)=str2double(filelist(i).name(1:end-13)); % extract out the trial number from each measurements file present in directory
end
trialNums = sort(dirTrialNums);
trialNums = trialNums(~isnan(trialNums));
trialNums = intersect(trialNums,bSession.trialNums); % try only the ones with behavior trials

includef=cell(size(trialNums,1),1);
for i = 1: length(trialNums)
    includef{i} = num2str(trialNums(i));
end


% %% Make whisker-pole touch space for each type of trial, from 10 randomly selected trials (of each type)
% Currently, only dealing with 4 types of trials: 'rc', 'rf', 'lc', 'lf'
% Should make something different for straight pole touch in S00. 
% 2017/04/11 JK

% Modifying to deal with different types of session, such as DISCRETE
% target angles and distractors.
% 2018/02/26 JK

% Determining the type of task and number of trial types within this session
servo_values = cellfun(@(x) x.servoAngle, bSession.trials);
servo_values = unique(servo_values(2:end)); % excluding trial #1, which is a dummy trial
distance_values = cellfun(@(x) x.motorDistance, bSession.trials);
distance_values = unique(distance_values(2:end)); % excluding trial #1, which is a dummy trial   
distance_values((distance_values == 0)) = []; % excluding catch trials
servo_distance_pair = cell(length(servo_values),length(distance_values));
for i = 1 : length(servo_values)
    for j = 1 : length(distance_values)
        servo_distance_pair{i,j} = [servo_values(i), distance_values(j)];
    end
end

% Initialize
steps_hp = cell(length(servo_values),length(distance_values));
num_points_in_hp = cell(length(servo_values),length(distance_values));
psi1 = zeros(length(servo_values),length(distance_values));
psi2 = zeros(length(servo_values),length(distance_values));
touch_hp = cell(length(servo_values),length(distance_values)); % touch hyperplanes
hp_peaks = cell(length(servo_values),length(distance_values)); % touch hyperplane peak points. 2 points for each hyperplane
trial_nums = cell(length(servo_values),length(distance_values));
thPolygon = cell(length(servo_values),length(distance_values));

thflist = dir([whisker_d, '*touch_hp.mat']);
if ~isempty(thflist)
    load(thflist(1).name);
end

%%
for iservo = 1 : length(servo_values)
    for idist = 1 : length(distance_values)              
        tt_ind = intersect(find(cellfun(@(x) (x.servoAngle == servo_values(iservo)),bSession.trials)), find(cellfun(@(x) (x.motorDistance == distance_values(idist)),bSession.trials)));
        tt_wst_ind = find(cellfun(@(x) ismember(x.trialNum, wstNums), bSession.trials));                
        tt_ind = intersect(tt_ind, tt_wst_ind);
        temp_files = cell(length(tt_ind),1);
        trial_nums{iservo,idist} = zeros(length(tt_ind),1);
        poleTipCoords = zeros(length(tt_ind),2); % [x coordinates, y coordinates] of poleUpFrames (average them)
        apPosition = zeros(length(tt_ind),1); % ap motor position from the behavior file
        for j = 1 : length(tt_ind)
            temp_files{j} = num2str(bSession.trials{tt_ind(j)}.trialNum);
            trial_nums{iservo,idist}(j) = bSession.trials{tt_ind(j)}.trialNum;

            load([temp_files{j},'_WST.mat']) % loading ws
            poleTipCoords(j,:) = mean(ws.topPix(ws.poleUpFrames,:));
            apPosition(j) = bSession.trials{tt_ind(j)}.motorApPosition;
        end

        ws = Whisker.WhiskerSignalTrialArray_2pad(whisker_d,'include_files',temp_files);            
        done_flag = 1; 
        psi1_polygon_answer = 'Yes'; % for re-drawing of polygon for psi1
        psi2Flip = 0;
        while (done_flag)
            intersect_3d_total = [];
            for tnum = 1 : length(ws.trials)
                    topInd = find(~isnan(ws.trials{tnum}.whiskerEdgeCoord(:,1)));
                    frontInd = find(~isnan(ws.trials{tnum}.whiskerEdgeCoord(:,2)));
                    noNaNInd = intersect(ws.trials{tnum}.poleUpFrames, intersect(topInd, frontInd));
                    intersect_3d_total = [intersect_3d_total; ws.trials{tnum}.whiskerEdgeCoord(noNaNInd,1), ws.trials{tnum}.whiskerEdgeCoord(noNaNInd,2), ws.trials{tnum}.apPosition(noNaNInd)];
            end
            intersect_3d_total = unique(round(intersect_3d_total,2),'rows');
            psi1_answer = 'Yes'; % for overall psi1 detection            
            while(strcmp(psi1_answer,'Yes'))                    
                while(strcmp(psi1_polygon_answer,'No'))
                    h2 = figure('units','normalized','outerposition',[0 0 1 1]); plot(intersect_3d_total(:,1), intersect_3d_total(:,2), 'k.', 'MarkerSize', 0.1), hold on
                    pre_poly = []; % points of the polygon
                    i = 1;
                    temp_point = ginput(1);
                    while(~isempty(temp_point)) % finish drawing polygon by pressing "enter"
                        pre_poly = [pre_poly; temp_point];
                        plot(pre_poly(i,1), pre_poly(i,2), 'bo', 'MarkerSize', 3)
                        if i > 1 
                            plot(pre_poly(i-1:i,1), pre_poly(i-1:i,2), 'b-')
                        end
                        i = i + 1;
                        temp_point = ginput(1);
                    end
                    plot([pre_poly(end,1);pre_poly(1,1)], [pre_poly(end,2);pre_poly(1,2)], 'b-')
%                         questTitle='Polygon pre-selection'; 
%                         start(timer('StartDelay',1,'TimerFcn',@(o,e)set(findall(0,'Tag',questTitle),'WindowStyle','normal')));         
                    psi1_polygon_answer = MFquestdlg([0.5,0.3], 'Is the drawing correct?', 'Region pre-selection', 'Yes', 'No', 'Yes');
                    switch psi1_polygon_answer
                        case 'Yes'
                            close all
                            in = inpolygon(intersect_3d_total(:,1), intersect_3d_total(:,2), pre_poly(:,1), pre_poly(:,2));
                            intersect_3d_crop = intersect_3d_total(in,:);
                        case 'No'
                            close(h2)
                            h2 = figure('units','normalized','outerposition',[0 0 1 1]); plot(intersect_3d_total(:,1), intersect_3d_total(:,2), 'k.', 'MarkerSize', 0.1), hold on
                    end                
                end            
                if exist('intersect_3d_crop','var') 
                    intersect_3d = intersect_3d_crop;
                    clear intersect_3d_crop
                else
                    intersect_3d = intersect_3d_total;
                end
                h1 = figure('units','normalized','outerposition',[0 0 1 1]); plot3(intersect_3d(:,1), intersect_3d(:,2), intersect_3d(:,3), 'k.', 'MarkerSize', 0.1)
                title(['Angle = ', num2str(servo_values(iservo)), ', Distance = ', num2str(distance_values(idist))]), xlabel('Top-view intersection coord'), ylabel('Front-view intersection coord'), zlabel('Pole position')

                %% when interested in certain points in the figure                    
%                     % 'oo': 245   246   353   382   436   558   559   
%    
%                         zvalue = 68000;
%                         tnumHigher = intersect(tt_ind, find(cellfun(@(x) x.motorApPosition < zvalue, bSession.trials)))
%                         tnumLower = intersect(tt_ind, find(cellfun(@(x) x.motorApPosition > zvalue, bSession.trials)))
%                         zvalue = 42410;
%                         tnum = intersect(tt_ind, find(cellfun(@(x) abs(x.motorApPosition - zvalue) < 10, bSession.trials)))
%                         tnum_ws = find(cellfun(@(x) x.trialNum == tnum(1), ws.trials))
%                         ws.trials{tnum_ws}.trackerFileName
%                         figure, plot3(ws.trials{tnum_ws}.whiskerEdgeCoord(:,1), ws.trials{tnum_ws}.whiskerEdgeCoord(:,2), 1:length(ws.trials{tnum_ws}.whiskerEdgeCoord(:,1)))
%                         xlabel('Top-view intersection coord'), ylabel('Front-view intersection coord'), zlabel('Frame #'), title(num2str(tnum))
                %% Calculate psi1 % takes ~ 15 sec
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Manual selection
                ind_opt = 1; % optimal peak index. Starting from 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                disp('Calculating psi1')

%                         angle_steps_pre = -5 : 179;
                angle_steps_pre = 1 : 185;
%                     angle_steps_pre = servo_values(iservo)+80 : servo_values(iservo)+105; % 180414 JK. contraint for initial guess to reduce computation time. From previous results of JK025 and JK027 with -5:185.
%                     angle_steps_pre = servo_values(iservo)-100 : servo_values(iservo)-75; % 180414 JK. contraint for initial guess to reduce computation time. From previous results of JK025 and JK027 with -5:185.
% Weird things happening with restricting pre angle steps. When angle_steps_pre > 180. 
                stds_pre = zeros(length(angle_steps_pre),1);
                for i = 1 : length(angle_steps_pre)
                    A = viewmtx(angle_steps_pre(i),0);
                    x4d = [intersect_3d, ones(size(intersect_3d,1),1)]';
                    x2d = A*x4d;
                    x2d_pix = [floor(x2d(1,:));floor(x2d(2,:)/100)];
                    x2d_dim = [max(x2d_pix(2,:)) - min(x2d_pix(2,:)) + 1, max(x2d_pix(1,:)) - min(x2d_pix(1,:)) + 1];
                    x2d_proj = zeros(x2d_dim);
                    j_offset = min(x2d_pix(1,:)) - 1;
                    i_offset = min(x2d_pix(2,:)) - 1;
                    for j = 1 : size(x2d_pix,2)
                        x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) = x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) + 1;
                    end
                    stds_pre(i) = std(x2d_proj(x2d_proj(:)~=0));
                end
                [~, I] = findpeaks(smooth(smooth(stds_pre)),'sortstr','descend');    

                answer = 'No'; % first answer for psi1 (Is psi1 correct?)
                while(strcmp(answer,'No'))
                    max_psi1_pre = angle_steps_pre(I(ind_opt));
                    h2 = figure('units','normalized','outerposition',[0.5 0 0.5 1]); subplot(1,2,1), plot(1:length(stds_pre), smooth(smooth(stds_pre,5))), hold on, plot(I(ind_opt),stds_pre(I(ind_opt)),'ro')
                    subplot(1,2,2), plot(1:length(stds_pre), stds_pre)
                    max_std = 0;
                    angle_steps = max_psi1_pre-0.9:0.1:max_psi1_pre+0.9; % resolution reduced to 0.1 from 0.01, flanking boundaries to 1 from 5 degrees for the sake of computation time 180414 JK
                    stds = zeros(length(angle_steps),1);
                    x2d_final = zeros(size(x2d_proj));
                    for i = 1:length(angle_steps)
                        A = viewmtx(angle_steps(i),0);
                        x4d = [intersect_3d, ones(size(intersect_3d,1),1)]';
                        x2d = A*x4d;
                        x2d_pix = [floor(x2d(1,:));floor(x2d(2,:)/100)]; % shrinking down to 100 times for better visualization
                        x2d_dim = [max(x2d_pix(2,:)) - min(x2d_pix(2,:)) + 1, max(x2d_pix(1,:)) - min(x2d_pix(1,:)) + 1];
                        x2d_proj = zeros(x2d_dim);

                        j_offset = min(x2d_pix(1,:)) - 1;
                        i_offset = min(x2d_pix(2,:)) - 1;
                        for j = 1 : length(x2d_pix)
                            x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) = x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) + 1;
                        end
                        temp_std = std(x2d_proj(x2d_proj(:)~=0));
                        
                        stds(i) = temp_std;
                        if temp_std > max_std
                            max_std = temp_std;
                            psi1(iservo,idist) = angle_steps(i);
                            x2d_final = x2d_proj;
                        end
                    end
                    psi1(iservo,idist) = psi1(iservo,idist)-90;
                    subplot(1,2,2), plot(1:length(stds), stds)

                    A = viewmtx(psi1(iservo,idist)+90,0);
                    x4d = [intersect_3d, ones(size(intersect_3d,1),1)]';
                    x2d = A*x4d;
                    h3 = figure('WindowStyle','normal','units','normalized','outerposition',[0 0 1 1]); plot(x2d(1,:), x2d(2,:),'k.', 'MarkerSize',3)        
                    answer = MFquestdlg([0.5, 0.3], 'Is psi1 correct?', 'whisker-pole intersection coordinate scatter', 'Yes', 'No', 'Yes');
                    switch answer
                        case 'Yes'
                            close all
                            psi1_answer = 'No'; % get out of this large while loop
                        case 'No' 
                            close(h2), close(h3)
                            ind_opt = ind_opt + 1;
                            if length(I) >= ind_opt
                                continue
                            else       
                                psi1_answer = MFquestdlg([0.5, 0.3],'Do you want to set (change) polygon drawing?', 'psi1 polygon', 'Yes', 'No', 'Yes');
                                if strcmp(psi1_answer,'Yes')
                                    answer = 'Yes'; % to get out of innner while loop 
                                    psi1_polygon_answer = 'No';
                                else
                                    error(['No observable psi1 at ' sessionName ' of ' mouseName])
                                end
                            end
                    end
                end
            end
                %% Manual check of the view
            %     A = viewmtx(I(6)+5,0);
            %     x4d = [intersect_3d, ones(size(intersect_3d,1),1)]';
            %     x2d = A*x4d;
            %     h3 = figure; plot(x2d(1,:), x2d(2,:),'k.', 'MarkerSize',3)

            %% Psi2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Manual selection
            answer1 = 'No'; answer2 = 'Yes';
            disp('Calculating psi2')
            x2d_flip = flip(x2d_final,1);
            while(strcmp(answer1,'No')) % answer1 = Is psi2 correct?
                h2 = figure('units','normalized','outerposition',[0 0 1 1]); imagesc(x2d_flip), hold on       
                BW = zeros(size(x2d_flip)) + 1;
                while(strcmp(answer2,'No')) % Draw polygon to select regions for radon transform
                    psi2_poly = []; 
                    i = 1;
                    temp_point = ginput(1); 
                    while(~isempty(temp_point)) % finish drawing polygon by pressing "enter"
                        psi2_poly = [psi2_poly;temp_point];                    
                        plot(psi2_poly(i,1), psi2_poly(i,2), 'yo', 'MarkerSize', 3)
                        if i > 1 
                            plot(psi2_poly(i-1:i,1), psi2_poly(i-1:i,2), 'y-')
                        end
                        i = i + 1;
                        temp_point = ginput(1);
                    end
                    plot([psi2_poly(end,1);psi2_poly(1,1)], [psi2_poly(end,2);psi2_poly(1,2)], 'y-')
                    answer2 = MFquestdlg([0.5, 0.3], 'Is the polygon right?', 'Polygon pre-selection', 'Yes', 'No', 'Yes');
                    switch answer2
                        case 'Yes'
                            close all
                            [sub_1, sub_2] = ind2sub(size(x2d_flip),1:length(x2d_flip(:)));
                            in = inpolygon(sub_2, sub_1, round(psi2_poly(:,1)), round(psi2_poly(:,2)));
                            BW = zeros(size(x2d_flip));
                            BW(in) = 1;                    
                        case 'No'
                            close(h2)
                            h2 = figure('units','normalized','outerposition',[0 0 1 1]); imagesc(x2d_flip), hold on
                    end
                end  
                x2d_edge = x2d_flip .* BW;
        %%
                theta = 0:0.01:180;
                R = radon(x2d_edge, theta);
                [~, max_ind] = max(std(R));
                psi2(iservo,idist) = (max_ind-1)*0.01;
                psi2(iservo,idist) = atand(tand(psi2(iservo,idist))/100); % psi2 adjusted because it was calculated with pole position divided by 100

                linex = zeros(3,2); 
                x2d_flip_1strow = x2d_flip(1,:);
                [~, linex(1,1)] = max(x2d_flip_1strow); x2d_flip_1strow(max(1,linex(1,1)-5):min(size(x2d_flip,2),linex(1,1)+5)) = 0; 
                [~, linex(2,1)] = max(x2d_flip_1strow); x2d_flip_1strow(max(1,linex(2,1)-5):min(size(x2d_flip,2),linex(2,1)+5)) = 0; 
                [~, linex(3,1)] = max(x2d_flip_1strow); 
                for i = 1 : 3
                    linex(i,2) = linex(i,1) + tand((max_ind-1)*0.01)*size(x2d_flip,1);
                end
                liney = [1 size(x2d_flip,1);1 size(x2d_flip,1);1 size(x2d_flip,1)];
                figure('units','normalized','outerposition',[0 0 1 1]), 
                subplot(121), imagesc(R), xlabel('Angle = 0:0.01:180'), ylabel('Projected values'), axis square
                subplot(122); imagesc(x2d_flip, [0 1000]), axis square, hold on, 
                for i = 1 : 3
                    line(linex(i,:),liney(i,:), 'LineWidth', 3, 'Color', [1 1 1])
                end

                answer1 = MFquestdlg([0.5, 0.3], 'Is psi2 correct?', 'whisker-pole intersection coordinate side-view', 'Yes', 'No', 'Yes');
                switch answer1
                    case 'Yes'
                        close all
                    case 'No'                 
                        answer3 = MFquestdlg([0.5, 0.3], 'Do you want to draw a region again?', 'Re-drawing the region', 'Yes', 'No', 'Yes');
                        close all
                        switch answer3
                            case 'Yes'
                                answer2 = 'No';
                            case 'No'                                     
                                answer4 = MFquestdlg([0.5, 0.3], 'Do you want to manually draw psi2?', 'Manual psi2', 'Yes', 'No', 'Yes');
                                if strcmp(answer4, 'Yes')                                        
                                    while (true)                                            
                                        figure('units','normalized','outerposition',[0 0 1 1])
                                        imagesc(x2d_flip, [0 1000]), axis square, hold on,
                                        [x,y] = ginput(2);                                            
                                        plot(x, y, 'wo-', 'MarkerSize', 5, 'LineWidth', 5)                                            
                                        answer5 = MFquestdlg([0.5, 0.3], 'Is this correct?', 'Manual psi2', 'Yes', 'No', 'Yes');
                                        if strcmp(answer5, 'Yes')
                                            angle = atand(diff(x)/diff(y)); % the image is 90 degrees rotated.
                                            psi2(iservo, idist) = atand(tand(angle)/100);
                                            answer1 = 'Yes';
                                            break
                                        else
                                            answer6 = MFquestdlg([0.5, 0.3], 'Do you want to specify the angle?', 'psi2', 'Yes', 'No', 'Yes');
                                            if strcmp(answer6, 'Yes')
                                                psi2(iservo, idist) = str2double(inputdlg('psi2 angle', 'psi2', 1, {''}, options));
                                                answer1 = 'Yes';
                                                break
                                            else
                                                error(['No optimal psi2 found in ' sessionName ' of ' mouseName])
                                            end
                                        end
                                    end

                                else
                                    answer6 = MFquestdlg([0.5, 0.3], 'Do you want to specify the angle?', 'psi2', 'Yes', 'No', 'Yes');
                                    if strcmp(answer6, 'Yes')
                                        psi2(iservo,idist) = str2double(inputdlg('psi2 angle', 'psi2', 1, {''}, options));
                                        answer1 = 'Yes';
                                        break
                                    else
                                        error(['No optimal psi2 found in ' sessionName ' of ' mouseName])
                                    end
                                end
                        end
                end
                if psi1(iservo,idist) > 90
                    psi2(iservo, idist) = -psi2(iservo, idist);
                    psi2Flip = 1;
                end
            end

            close all
            %% Calculate touch hyperplanes
            answer7 = 'Yes';
            while strcmp(answer7, 'Yes')
                answer7 = 'No'; % stay in this while loop only when certain condition is met (psi2 flip for some 90 degrees)

                disp('Sweeping the hyperplane')
                maxdist = ceil(sqrt(max(intersect_3d_total(:,1).^2) + max(intersect_3d_total(:,2).^2))); xmin = -maxdist; xmax = maxdist;                        
                zmin_data = min(intersect_3d_total(:,3)); zmax_data = max(intersect_3d_total(:,3)); zdiff = zmax_data - zmin_data; zmax = zmax_data + zdiff; zmin = zmin_data - zdiff;
                ymin_data = min(intersect_3d_total(:,2)); ymax_data = max(intersect_3d_total(:,2));
                z = zmin:zmax;
                xyz = zeros((length(z))*(xmax-xmin+1),3);
                for i = xmin:xmax
                    xyz((i-xmin)*length(z)+1 : (i-xmin+1)*length(z),:) = [ones(length(z),1)*i + mean(intersect_3d_total(:,1)), zeros(length(z),1) + mean(intersect_3d_total(:,2)), z'];
                end

                hp_decision = 'No';
                while(strcmp(hp_decision,'No'))
%                             [xyz_psi1, ~, ~] = AxelRot(xyz',psi1(iservo,idist),[0 0 1], 0); % rotate psi1 degrees counterclockwise around z axis
                    [xyz_psi1, ~, ~] = AxelRot(xyz',psi1(iservo,idist),[0 0 1], [mean(intersect_3d_total(:,1)) mean(intersect_3d_total(:,2)) 0]); % rotate psi1 degrees counterclockwise around z axis

                    zcenter = floor(mean([zmax_data, zmin_data]));
%                             x0 = [0 0 zcenter];
                    x0 = [mean(intersect_3d_total(:,1)) mean(intersect_3d_total(:,2)) zcenter];                            
                    if psi1(iservo,idist) == 90
                        u = [0 1 0];
                    else
                        u = [1 tand(psi1(iservo,idist)) 0];
                    end
                    [xyz_psi2, ~, ~] = AxelRot(xyz_psi1, psi2(iservo,idist), u, x0); 

                    xyz_psi2(:,xyz_psi2(3,:) < zmin_data) = [];
                    xyz_psi2(:,xyz_psi2(3,:) > zmax_data) = [];
                    xyz_psi2(:,xyz_psi2(2,:) < ymin_data) = [];
                    xyz_psi2(:,xyz_psi2(2,:) > ymax_data) = [];
                    if isempty(xyz_psi2) 
                        [xyz_psi2, ~, ~] = AxelRot(xyz_psi1, psi2(iservo,idist), u, x0); 
                        xyz_psi2(:,xyz_psi2(2,:) < ymin_data) = [];
                        xyz_psi2(:,xyz_psi2(2,:) > ymax_data) = [];
                    end
                    %%
%                     figure, plot3(intersect_3d_total(:,1),intersect_3d_total(:,2), intersect_3d_total(:,3),'k.', 'MarkerSize',3), xlabel('top'), ylabel('front'), zlabel('pos'), hold on
%                     zmaxind = find(xyz_psi2(3,:) == zmax_data); zminind = find(xyz_psi2(3,:) == zmin_data);
%                     for zi = 1 : length(zmaxind)
%                         xind = find(abs(xyz_psi2(1,:) - xyz_psi2(1,zmaxind(zi))) < 1);
%                         [~, minzind] = min(xyz_psi2(3,xind));                        
%                         plot3(xyz_psi2(1,[zmaxind(zi), xind(minzind)]), xyz_psi2(2,[zmaxind(zi), xind(minzind)]), xyz_psi2(3,[zmaxind(zi), xind(minzind)]), 'r-')
%                     end

                    %% test
                    figure('units','normalized','outerposition',[0 0.3 0.5 0.5]), 
                    plot3(intersect_3d_total(:,1),intersect_3d_total(:,2), intersect_3d_total(:,3),'k.', 'MarkerSize',3), xlabel('top'), ylabel('front'), zlabel('pos'), hold on                    
                    plot3(xyz_psi2(1,:), xyz_psi2(2,:), xyz_psi2(3,:), 'r.')
                    %% ~ 0.5 min (depending on the length of "steps" and the size of xyz_psi2)
                    intersect_pix = round(intersect_3d_total);

                    if isempty(steps_hp{iservo, idist})
                        steps_hp{iservo, idist} = -ppm*4 : -ppm*2 + 20;
                        steps_hp{iservo, idist} = round(steps_hp{iservo, idist});
                    end
                    num_points = zeros(length(steps_hp{iservo, idist}),1);
                    parfor i = 1:length(steps_hp{iservo, idist}) % this is time consuming...
                        hp = round(xyz_psi2);
                        hp(1,:) = hp(1,:)+ steps_hp{iservo, idist}(i);    
                        num_points(i) = sum(ismember(intersect_pix, hp','rows'));
                    end

                    h1 = figure('WindowStyle','normal','units','normalized','outerposition',[0.5 0.3 0.5 0.5]); plot(steps_hp{iservo, idist},num_points(:), 'k-', 'LineWidth', 3), xlabel('translocation (pix)'), ylabel('# intersection')

                    answer = MFquestdlg([0.5, 0.3], 'Does the result look correct?', 'Touch hyperplane peaks', 'Yes', 'No', 'Yes');            
                    peak_answer = 'Yes';
                    while(strcmp(peak_answer,'Yes'))
                        if strcmp(answer,'Yes')
                            while true
                                datacursormode(h1,'on')
                                hp_peak_ans = inputdlg({'Left peak','Right peak'},'What are the peak points?',1,{'',''},optionsFin);
                                if ~isempty(hp_peak_ans) && ~isempty(hp_peak_ans{1}) && ~isempty(hp_peak_ans{2})
                                    break
                                end
                            end
                            hp_peaks{iservo, idist} = [str2double(hp_peak_ans{1}) str2double(hp_peak_ans{2})];                    

                            % Final confirmation
                            % project the peak hyperplanes and all coordinates onto psi1 psi2 view
                            %%
                            h2 = figure('units','normalized','outerposition',[0 0 1 1]); 
                            if psi2Flip
                                A = viewmtx(psi1(iservo,idist),-90+psi2(iservo,idist));
                            else
                                A = viewmtx(psi1(iservo,idist),90-psi2(iservo,idist));
                            end
                            intersect_4d = [intersect_3d_total, ones(size(intersect_3d_total,1),1)]';
                            intersect_2d = A*intersect_4d;
                            intersect_2d = unique(round(intersect_2d(1:2,:)',2),'rows');
                            th_4d1 = [xyz_psi2(1,:) + hp_peaks{iservo, idist}(1);xyz_psi2(2:3,:);ones(1,size(xyz_psi2,2))];
                            th_2d1 = A*th_4d1;
                            th_2d1 = unique(th_2d1(1:2,:)','rows');
                            th_4d2 = [xyz_psi2(1,:) + hp_peaks{iservo, idist}(2);xyz_psi2(2:3,:);ones(1,size(xyz_psi2,2))];
                            th_2d2 = A*th_4d2;
                            th_2d2 = unique(th_2d2(1:2,:)','rows');
                            scatter(intersect_2d(:,1),intersect_2d(:,2),'k.'), hold on, scatter(th_2d1(:,1), th_2d1(:,2),'r.'), scatter(th_2d2(:,1), th_2d2(:,2),'r.')
                            %%
                            answer2 = MFquestdlg([0.5, 0.3], 'Is the result REALLY correct?', 'Final Confirmation', 'Yes', 'No', 'Yes');
                            switch answer2
                                case 'Yes'
                                    close(h2);
                                    peak_answer = 'No'; % get out of peak while
                                    hp_decision = 'Yes';
                                    done_flag = 0; % the whole precedure is finally done. get out of while(done_flag) loop.
                                case 'No'
                                    answer3 = MFquestdlg([0.5, 0.3], 'Do you want to change the steps?', 'Touch hyperplane sweep steps', 'Yes', 'No', 'Yes');
                                    switch answer3
                                        case 'Yes'                                    
                                            step_boundary_cell = inputdlg({'First step','Last step'},'What are the sweep boundaries?',1,{'',''},options);
                                            steps_hp{iservo, idist} = str2double(step_boundary_cell{1}):str2double(step_boundary_cell{2});
                                            peak_answer = 'No'; % get out of peak while
                                            close all
                                        case 'No'
                                            peak_answer = MFquestdlg([0.5, 0.3], 'Do you want to change the peak points?', 'Touch hyperplane peaks', 'Yes', 'No', 'Yes');
                                            if strcmp(peak_answer,'Yes')                                        
                                                close(h2);
                                            else
                                                psi1_return_answer = MFquestdlg([0.5, 0.3], 'Do you want to draw region for psi1?', 'Return to psi1 region selection', 'Yes', 'No', 'Yes');
                                                if strcmp(psi1_return_answer, 'Yes')
                                                    close all
                                                    psi1_polygon_answer = 'No';
                                                    peak_answer = 'No';
                                                    hp_decision = 'Yes';
                                                else
                                                    psi1_adjust = inputdlg('How much change in psi1? (Negative for CW)', 'psi1 manual adjustment', 1, {''}, options);
                                                    if psi1_adjust{1}
                                                        psi1(iservo,idist) = psi1(iservo,idist) + str2double(psi1_adjust{1});
                                                    else
                                                        error(['No optimal hyperplanes peaks found in ' sessionName ' of ' mouseName])
                                                    end
                                                end
                                            end
                                    end
                            end                    
                        else % answer = 'No' to question 'Does the result look correct?'
                            answer3 = MFquestdlg([0.5, 0.3], 'Do you want to change the steps?', 'Touch hyperplane sweep steps', 'Yes', 'No', 'Yes');
                            switch answer3
                                case 'Yes'                            
                                    step_boundary_cell = inputdlg({'First step','Last step'},'What are the sweep boundaries?',1,{'',''},options);
                                    steps_hp{iservo, idist} = str2double(step_boundary_cell{1}):str2double(step_boundary_cell{2});
                                    peak_answer = 'No';
                                    answer7 = 'No';
                                    close all
                                case 'No'
                                    answer7 = MFquestdlg([0.5, 0.3], 'Do you want to flip psi2?', 'Wierd psi2 error', 'Yes', 'No', 'Yes');
                                    switch answer7
                                        case 'Yes'
                                            psi2Flip = 1-psi2Flip;
                                            psi2(iservo, idist) = -psi2(iservo, idist);
                                            peak_answer = 'No';
                                        case 'No'
        %                                     questTitle = 'Return to psi1 polygon'; 
        %                                     start(timer('StartDelay',1,'TimerFcn',@(o,e)set(findall(0,'Tag',questTitle),'WindowStyle','normal')));
                                            psi1_return_answer = MFquestdlg([0.5, 0.3], 'Do you want to draw a region for psi1?', 'Return to psi1 region selection', 'Yes', 'No', 'Yes');
                                            if strcmp(psi1_return_answer, 'Yes')
                                                close all
                                                psi1_polygon_answer = 'No';
                                                peak_answer = 'No';
                                                hp_decision = 'Yes';
                                            else
                                                psi1_adjust = inputdlg('How much change in psi1?', 'psi1 manual adjustment', 1, {''}, options);
                                                if psi1_adjust{1}
                                                    psi1(iservo,idist) = psi1(iservo,idist) + str2double(psi1_adjust{1});
                                                else
                                                    error(['No optimal hyperplanes peaks found in ' sessionName ' of ' mouseName])
                                                end
                                            end       
                                    end
                            end
                        end
                    end
                end
            end
        end
        %% save parameters
        close all
        num_points_in_hp{iservo, idist} = num_points;
        touch_hp{iservo, idist} = xyz_psi2; % Don't round them! (at least at this saving process)                
        th_2d = [th_2d1; th_2d2];
        cvh = convhull(th_2d);
        thPolygon{iservo, idist} = th_2d(cvh,:);

        fprintf('%s %s trial type #%d/%d processed\n',mouseName, sessionName, (iservo-1)*length(distance_values) + idist, length(servo_values) * length(distance_values))        
    end
end
%%
save([whisker_d mouseName sessionName '_touch_hp.mat'],'touch_hp','num_points_in_hp','steps_hp','hp_peaks', 'psi1', 'psi2', 'servo_distance_pair', 'thPolygon')
fprintf('%s %s hp_peaks saved\n', mouseName, sessionName)


