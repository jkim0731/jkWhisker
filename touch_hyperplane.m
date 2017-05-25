%% Setup whisker array builder 
behavior_base_dir = 'Z:\Data\2p\soloData\';
whisker_base_dir = 'Z:\Data\Video\JK\';

mice = {'AH0648','AH0650','AH0651','AH0652','AH0653'};

sessionNum = [4, 6:15];
for sessionInd = 1 : length(sessionNum)
% for sessionInd = 1
    mouseName = 'AH0648';
    sessionName = sprintf('S%02d',sessionNum(sessionInd));

    behavior_d = [behavior_base_dir mouseName '\'];
    whisker_d = [whisker_base_dir mouseName sessionName '\'];

    % if exist([whisker_d, 'touch_hp.mat'],'file')
    %     error('touch_hp.mat exists.')    
    % end

    if exist('b','var')
        if strcmp(b{1}.mouseName, mouseName)
            disp('using the same behavior file')
        else
            disp('loading a new behavior file')
            load([behavior_d 'behavior.mat']) % loading b of the mouse (all the sessions)
        end
    else
        load([behavior_d 'behavior.mat']) % loading b of the mouse (all the sessions)
    end
    b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
    b_session = b{b_ind};

    load_fn = [mouseName sessionName '_post.mat'];
    load([whisker_d load_fn]); % loading errorlist

    if ~isempty(b_ind) % try only the ones with behavior session
        % %%
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
        dirTrialNums = setdiff(dirTrialNums,errorlist);
        trialNums = sort(dirTrialNums);
        trialNums = trialNums(~isnan(trialNums));
        trialNums = intersect(trialNums,b{b_ind}.trialNums); % try only the ones with behavior trials

        includef=cell(size(trialNums,1),1);
        for i = 1: length(trialNums)
            includef{i} = num2str(trialNums(i));
        end
    end

    % %% Make whisker-pole touch space for each type of trial, from 10 randomly selected trials (of each type)
    % Currently, only dealing with 4 types of trials: 'rc', 'rf', 'lc', 'lf'
    % Should make something different for straight pole touch in S00. 
    % 2017/04/11 JK
    trial_types = {'rc', 'rf', 'lc', 'lf'};
    steps_hp = cell(1,length(trial_types));
    num_points_in_hp = cell(1,length(trial_types));
    tt_ind = cell(1,length(trial_types));
    wl_array = cell(1,length(trial_types));
    % touch_points = cell(1,length(trial_types));
    touch_hp = cell(1,length(trial_types)); % touch hyperplanes
    hp_peaks = cell(1,length(trial_types)); % touch hyperplane peak points. 2 points for each hyperplane
    % %%
    % load('wl_array.mat')

    for trial_type_num = 1 : length(trial_types)    
    % trial_type_num = 1
        tt_ind{trial_type_num} = find(cellfun(@(x) strcmp(x.trialType,trial_types{trial_type_num}),b_session.trials));
        temp_files = cell(length(tt_ind{trial_type_num}),1);
        for j = 1 : length(tt_ind{trial_type_num})
            temp_files{j} = num2str(tt_ind{trial_type_num}(j));
        end
        wl = Whisker.WhiskerTrialLiteArray(whisker_d,'include_files',temp_files);
        wl_array{trial_type_num} = wl;
    end
    %%
    for trial_type_num = 1:4
    % trial_type_num = 3  
        intersect_3d = [];
        wl = wl_array{trial_type_num};        
        for tnum = 1 : length(wl.trials)
            try        
                top_ind = find(~isnan(wl.trials{tnum}.intersect_coord(:,1)));
                front_ind = find(~isnan(wl.trials{tnum}.intersect_coord(:,2)));
                intersect_ind = intersect(wl.trials{tnum}.pole_available_timepoints,intersect(top_ind,front_ind));
                intersect_3d = [intersect_3d; wl.trials{tnum}.intersect_coord(intersect_ind,1), wl.trials{tnum}.intersect_coord(intersect_ind,2), ones(length(intersect_ind),1)*wl.trials{tnum}.pole_pos];
            catch
                fprintf('Skipping trial #%d because of index problems \n',tnum);        
            end
        end
        psi1_answer = 'Yes'; % for overall psi1 detection
        psi1_polygon_answer = 'Yes'; % for re-drawing of polygon for psi1
        while(strcmp(psi1_answer,'Yes'))                    
            while(strcmp(psi1_polygon_answer,'No'))
                h2 = figure('units','normalized','outerposition',[0 0 1 1]); plot(intersect_3d(:,1), intersect_3d(:,2), 'k.', 'MarkerSize', 0.1), hold on
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
                psi1_polygon_answer = questdlg('Is the polygon right?', 'Polygon pre-selection', 'Yes', 'No', 'Yes');
                switch psi1_polygon_answer
                    case 'Yes'
                        close all
                    case 'No'
                        close(h2)
                        h2 = figure('units','normalized','outerposition',[0 0 1 1]); plot(intersect_3d(:,1), intersect_3d(:,2), 'k.', 'MarkerSize', 0.1), hold on
                end
                in = inpolygon(intersect_3d(:,1), intersect_3d(:,2), pre_poly(:,1), pre_poly(:,2));
                intersect_3d = intersect_3d(in,:);
            end            
            h1 = figure; plot3(intersect_3d(:,1), intersect_3d(:,2), intersect_3d(:,3), 'k.', 'MarkerSize', 0.1)
            title(wl.trials{1}.trial_type), xlabel('Top-view intersection coord'), ylabel('Front-view intersection coord'), zlabel('Pole position')

            %% when interested in certain points in the figure
            % ttype = 4;
            % zvalue = 90050;
            % tnum = find(cellfun(@(x) abs(x.pole_pos - zvalue) < 10, wl_array{ttype}.trials))
            % wl_array{ttype}.trials{tnum(1)}.trackerFileName
            % figure, plot3(wl_array{ttype}.trials{tnum(1)}.intersect_coord(:,1), wl_array{ttype}.trials{tnum(1)}.intersect_coord(:,2), 1:length(wl_array{ttype}.trials{tnum(1)}.intersect_coord(:,1)))
            % xlabel('Top-view intersection coord'), ylabel('Front-view intersection coord'), zlabel('Frame #')
            %% Calculate psi1 % takes ~ 15 sec
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Manual selection
            ind_opt = 1; % optimal peak index. Starting from 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            disp('Calculating psi1')

            max_zeros = 0;
            angle_steps_pre = 0:180;
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
                for j = 1 : length(x2d_pix)
                    x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) = x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) + 1;
                end
                stds_pre(i) = std(x2d_proj(find(x2d_proj(:))));
        %         stds_pre(i) = std(x2d_proj(:));
            end
            [P, I] = findpeaks(smooth(smooth(stds_pre)),'sortstr','descend');    

            answer = 'No';
            while(strcmp(answer,'No'))
                max_psi1_pre = angle_steps_pre(I(ind_opt));
                h2 = figure; subplot(1,2,1), plot(1:length(stds_pre), smooth(smooth(stds_pre,5))), hold on, plot(I(ind_opt),stds_pre(I(ind_opt)),'ro')
                subplot(1,2,2), plot(1:length(stds_pre), stds_pre)
                max_std = 0;
                psi1 = 0;
                angle_steps = max_psi1_pre-4.99:0.01:max_psi1_pre+4.99;
                stds = zeros(length(angle_steps),1);
                x2d_final = zeros(size(x2d_proj));
                for i = 1:length(angle_steps)
                    A = viewmtx(angle_steps(i),0);
                    x4d = [intersect_3d, ones(size(intersect_3d,1),1)]';
                    x2d = A*x4d;
                    x2d_pix = [floor(x2d(1,:));floor(x2d(2,:)/100)];
                    x2d_dim = [max(x2d_pix(2,:)) - min(x2d_pix(2,:)) + 1, max(x2d_pix(1,:)) - min(x2d_pix(1,:)) + 1];
                    x2d_proj = zeros(x2d_dim);

                    j_offset = min(x2d_pix(1,:)) - 1;
                    i_offset = min(x2d_pix(2,:)) - 1;
                    for j = 1 : length(x2d_pix)
                        x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) = x2d_proj(x2d_pix(2,j) - i_offset, x2d_pix(1,j) - j_offset) + 1;
                    end
                    temp_std = std(x2d_proj(find(x2d_proj(:))));
        %             temp_std = std(x2d_proj(:));

                    stds(i) = temp_std;
                    if temp_std > max_std
                        max_std = temp_std;
                        psi1 = angle_steps(i);
                        x2d_final = x2d_proj;
                    end
                end
                psi1 = psi1-90;
                subplot(1,2,2), plot(1:length(stds), stds)

                A = viewmtx(psi1+90,0);
                x4d = [intersect_3d, ones(size(intersect_3d,1),1)]';
                x2d = A*x4d;
                h3 = figure; plot(x2d(1,:), x2d(2,:),'k.', 'MarkerSize',3)
                questTitle='whisker-pole intersection coordinate scatter'; start(timer('StartDelay',1,'TimerFcn',@(o,e)set(findall(0,'Tag',questTitle),'WindowStyle','normal')));         
                answer = questdlg('Is psi1 correct?', questTitle, 'Yes', 'No', 'Yes');
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
                            psi1_answer = questdlg('Do you want to change your polygon drawing?', 'psi1 polygon', 'Yes', 'No', 'Yes');
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
        while(strcmp(answer1,'No'))
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
                answer2 = questdlg('Is the polygon right?', 'Polygon pre-selection', 'Yes', 'No', 'Yes');
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
            figure, imagesc(x2d_edge)
    %%
            theta = 0:0.01:180;
    %         x2d_flip = flip(x2d_final,1);
            R = radon(x2d_edge, theta);
    %         [~, max_ind] = max(max(R));
            [~, max_ind] = max(std(R));
            psi2 = (max_ind-1)*0.01;
            psi2 = atand(tand(psi2)/100); % psi2 adjusted because it was calculated with pole position divided by 100

            linex = zeros(3,2); 
            linex(1,1) = 1; linex(2,1) = floor(size(x2d_flip,2)/2); linex(3,1) = size(x2d_flip,2);
            for i = 1 : 3
                linex(i,2) = linex(i,1) + tand((max_ind-1)*0.01)*size(x2d_flip,1);
            end
            liney = [1 size(x2d_flip,1);1 size(x2d_flip,1);1 size(x2d_flip,1)];
            figure, 
            subplot(121), imagesc(R), xlabel('Angle = 0:0.01:180'), ylabel('Projected values'), axis square
            subplot(122), imagesc(x2d_flip, [0 1000]), axis square, hold on, 
            for i = 1 : 3
                line(linex(i,:),liney(i,:), 'LineWidth', 3, 'Color', [1 1 1])
            end

            questTitle='whisker-pole intersection coordinate side-view'; start(timer('StartDelay',1,'TimerFcn',@(o,e)set(findall(0,'Tag',questTitle),'WindowStyle','normal')));         
            answer1 = questdlg('Is psi2 correct?', questTitle, 'Yes', 'No', 'Yes');
            switch answer1
                case 'Yes'
                    close all
                case 'No'                 
                    answer3 = questdlg('Do you want to select the polygon again?', ' ', 'Yes', 'No', 'Yes');
                    close all
                    switch answer3
                        case 'Yes'
                            answer2 = 'No';
                        case 'No' 
                            error(['No optimal psi2 found in ' sessionName ' of ' mouseName])
                    end
            end
        end
        %% Calculate touch hyperplanes
        disp('Sweeping the hyperplane')
        xmin = -200; xmax = 200; zmin_data = min(intersect_3d(:,3)); zmax_data = max(intersect_3d(:,3)); zdiff = zmax_data - zmin_data; zmax = zmax_data + zdiff; zmin = zmin_data - zdiff;
        ymin_data = min(intersect_3d(:,2)); ymax_data = max(intersect_3d(:,2));
        z = zmin:zmax;
        xyz = zeros((length(z))*(xmax-xmin+1),3);
        for i = xmin:xmax
            xyz((i-xmin)*length(z)+1 : (i-xmin+1)*length(z),:) = [ones(length(z),1)*i, zeros(length(z),1), z'];
        end

        [xyz_psi1, ~, ~] = AxelRot(xyz',psi1,[0 0 1], 0); % rotate psi1 degrees counterclockwise around z axis

        zcenter = floor(mean([zmax_data, zmin_data]));
        x0 = [0 0 zcenter];
        u = [1 tand(psi1) 0];
        [xyz_psi2, ~, ~] = AxelRot(xyz_psi1, psi2, u, x0); 
        xyz_psi2(:,xyz_psi2(3,:) < zmin_data) = [];
        xyz_psi2(:,xyz_psi2(3,:) > zmax_data) = [];
        xyz_psi2(:,xyz_psi2(2,:) < ymin_data) = [];
        xyz_psi2(:,xyz_psi2(2,:) > ymax_data) = [];

        figure, plot3(intersect_3d(:,1),intersect_3d(:,2), intersect_3d(:,3),'k.', 'MarkerSize',3), xlabel('top'), ylabel('front'), zlabel('pos'), hold on
        plot3(xyz_psi2(1,:), xyz_psi2(2,:), xyz_psi2(3,:), 'r.', 'MarkerSize',3)

        %% ~ 0.5 min (depending on the length of "steps" and the size of xyz_psi2)
        %%%%%%%%%%%%%%%%%%%%%% manual selection
        switch trial_type_num
            case 1
                steps = 10:70;
            case 2
                steps = 20:80;
            case 3
                steps = 140:200;
            case 4
                steps = 140:200;
        end
        %%%%%%%%%%%%%%%%%%%%%% try as short as possible to reduce time next step
        %%
        hp_decision = 'No';
        while(strcmp(hp_decision,'No'))
            intersect_pix = round(intersect_3d);

            num_points = zeros(length(steps),1);
            parfor i = 1:length(steps) % this is time consuming...
                hp = round(xyz_psi2);
                hp(1,:) = hp(1,:)+steps(i);    
                num_points(i) = sum(ismember(intersect_pix, hp','rows'));
            end

            h1 = figure; plot(steps,num_points(:), 'k-', 'LineWidth', 3), xlabel('translocation (pix)'), ylabel('# intersection')

            options.WindowStyle = 'normal';
            questTitle='Touch hyperplane peaks'; start(timer('StartDelay',1,'TimerFcn',@(o,e)set(findall(0,'Tag',questTitle),'WindowStyle','normal')));         
            answer = questdlg('Does the result look correct?', questTitle, 'Yes', 'No', 'Yes');            
            peak_answer = 'Yes';
            while(strcmp(peak_answer,'Yes'))
                if strcmp(answer,'Yes')
                    datacursormode(h1,'on')
                    hp_peaks{trial_type_num} = str2num(cell2mat(inputdlg({'Left peak','Right peak'},'What are the peak points?',1,{'',''},options)))';                                    

                    % Final confirmation
                    h2 = figure('units','normalized','outerposition',[0 0 1 1]); 
                    subplot(121), plot3(intersect_3d(:,1),intersect_3d(:,2), intersect_3d(:,3),'k.', 'MarkerSize',1), xlabel('top'), ylabel('front'), zlabel('pos'), hold on
                    plot3(xyz_psi2(1,:) + hp_peaks{trial_type_num}(1)-2, xyz_psi2(2,:), xyz_psi2(3,:), 'r.', 'MarkerSize', 1) 
                    subplot(122), plot3(intersect_3d(:,1),intersect_3d(:,2), intersect_3d(:,3),'k.', 'MarkerSize',1), xlabel('top'), ylabel('front'), zlabel('pos'), hold on
                    plot3(xyz_psi2(1,:) + hp_peaks{trial_type_num}(2)+2, xyz_psi2(2,:), xyz_psi2(3,:), 'r.', 'MarkerSize', 1) 
                    options.WindowStyle = 'normal';
                    questTitle='Final Confirmation'; start(timer('StartDelay',1,'TimerFcn',@(o,e)set(findall(0,'Tag',questTitle),'WindowStyle','normal')));         
                    answer2 = questdlg('Is the result REALLY correct?', questTitle, 'Yes', 'No', 'Yes');
                    switch answer2
                        case 'Yes'
                            close all
                            peak_answer = 'No'; % get out of peak while
                            hp_decision = 'Yes';
                        case 'No'
                            answer3 = questdlg('Do you want to change the steps?', 'Touch hyperplane sweep steps', 'Yes', 'No', 'Yes');
                            switch answer3
                                case 'Yes'                                    
                                    step_boundary_cell = inputdlg({'First step','Last step'},'What are the sweep boundaries?',1,{'',''},options);
                                    steps = str2double(step_boundary_cell{1}):str2double(step_boundary_cell{2});
                                    peak_answer = 'No'; % get out of peak while
                                    close all
                                case 'No'
                                    peak_answer = questdlg('Do you want to change the peak points?', 'Touch hyperplane peaks', 'Yes', 'No', 'Yes');
                                    if strcmp(peak_answer,'Yes')                                        
                                        close(h2);
                                    else
                                        error(['No optimal hyperplanes peaks found in ' sessionName ' of ' mouseName])
                                    end
                            end
                    end                    
                else % answer = 'No' to question 'Does the result look correct?'
                    answer3 = questdlg('Do you want to change the steps?', 'Touch hyperplane sweep steps', 'Yes', 'No', 'Yes');
                    switch answer3
                        case 'Yes'                            
                            step_boundary_cell = inputdlg({'First step','Last step'},'What are the sweep boundaries?',1,{'',''},options);
                            steps = str2double(step_boundary_cell{1}):str2double(step_boundary_cell{2});
                            peak_answer = 'No';
                            close all
                        case 'No'
                            error(['No optimal hyperplanes peaks found in ' sessionName ' of ' mouseName])
                    end
                end
            end
        end

        %% save parameters
        steps_hp{trial_type_num} = steps;
        num_points_in_hp{trial_type_num} = num_points;
        touch_hp{trial_type_num} = xyz_psi2; % Don't round them! (at least at this saving process)
        sprintf('trial type #%d processed',trial_type_num)
    end
    %%
    save([whisker_d 'touch_hp.mat'],'touch_hp','num_points_in_hp','steps_hp','hp_peaks')
    disp('hp_peaks saved')
end