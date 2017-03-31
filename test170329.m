tic
session_names = {'AH0650S02', 'AH0651S02', 'AH0652S04', 'AH0653S03'};
% session_names = {'AH0648S02'};
r_in_mm = [0.5:0.5:5];

for session_ind = 1 : length(session_names)
%     d = ['Z:\Data\Video\JK\', session_names{session_ind}, '\'];
    d = ['/mnt/Data/Video/JK/', session_names{session_ind}, '/'];
    cd(d);
    load([session_names{session_ind},'mask.mat'])
    filelist=dir([d '*.measurements']);
    dirTrialNums=zeros(1,size(filelist,1));
    for i=1:length(filelist);
        dirTrialNums(i)=str2double(filelist(i).name(1:end-13)); % extract out the trial number from each measurements file present in directory
    end
    trialNums = sort(dirTrialNums);
    trialNums = trialNums(~isnan(trialNums));
    includef=cell(size(trialNums,1),1);
    for i = 1: length(trialNums)
        includef{i} = num2str(trialNums(i));
    end
    v = VideoReader([includef{1} '.mp4']);
    vv = read(v,1);
    vheight = size(vv,1);
    vwidth = size(vv,2);
    mouseName = session_names{session_ind}(1:6);
    sessionName = session_names{session_ind}(7:9);
    
    thab_0 = cell(length(trialNums),length(r_in_mm)); % thetaAtBase for tid 0
    thab_1 = cell(length(trialNums),length(r_in_mm)); % thetaAtBase for tid 1
    dkp_0 = cell(length(trialNums),length(r_in_mm)); % deltaKappa for tid 0
    dkp_1 = cell(length(trialNums),length(r_in_mm)); % deltaKappa for tid 1

    xy = cell(length(r_in_mm),2);
%     if session_ind ~= 1
%         Whisker.makeAllDirectory_WhiskerTrial(d,[0 1],'mask', {[maskx(1,:);masky(1,:)],[maskx(2,:);masky(2,:)]},...
%             'trial_nums',trialNums,'include_files',includef,...
%             'barRadius',15.3,'faceSideInImage', 'bottom', 'framePeriodInSec',.0032,...
%             'imagePixelDimsXY',[vwidth vheight],'pxPerMm',26.23,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward')
% 
%         Whisker.makeAllDirectory_WhiskerSignalTrial(d,'include_files',includef,'polyRoiInPix',[10 150],'follicleExtrapDistInPix',40);
%     end

    for i = 1 : length(r_in_mm)    
    %     Whisker.makeAllDirectory_WhiskerTrialLiteI(d,'include_files',includef(temp_tn),'r_in_mm',1.2,'calc_forces',false);
        Whisker.makeAllDirectory_WhiskerTrialLiteI(d,'include_files',includef,'r_in_mm',r_in_mm(i),'calc_forces',false);
        wl = Whisker.WhiskerTrialLiteArray(d);
        for j = 1 : length(trialNums)
            thab_0{j,i} = wl.trials{j}.thetaAtBase{1};
            thab_1{j,i} = wl.trials{j}.thetaAtBase{2};
            dkp_0{j,i} = wl.trials{j}.deltaKappa{1};
            dkp_1{j,i} = wl.trials{j}.deltaKappa{2};
        end
    %     load([num2str(temp_tn+1),'_WST.mat'])
    %     tp = [0 4];
    %     if i == 1
    %         ws.plot_fitted_whisker_time_projection(0,'k',tp), grid on, hold on;
    %         ws.plot_fitted_whisker_time_projection(1,'k',tp)
    %         ws.plot_fitted_whisker_ROI_time_projection(0,'r',tp)
    %         ws.plot_fitted_whisker_ROI_time_projection(1,'r',tp)
    %         ws.plot_mask(0,'g',tp)
    %         ws.plot_mask(1,'g',tp)
    %     end
    %     [~,~,y,x,~] = get_theta_kappa_at_roi_point(ws,0,r_in_mm(i));
    %     plot(x,y,'b.')
    %     [~,~,y,x,~] = get_theta_kappa_at_roi_point(ws,1,r_in_mm(i));
    %     plot(x,y,'b.')
    %     tid = [0 1]; % Set trajectory ID to view
    % Whisker.viewdouble_WhiskerTrialLiteArray(wl,tid)
    end
    %%
    % f = [1:length(thab{1,1})]; % because I know that there is no non-traced frame
    % 
    % figure,
    % subplot(211)
    % for i = 1 : length(r_in_mm)
    %     p = plot(f,thab{i,1});
    %     color_value = 0.5 + (i'AH0648S02',-1)/(length(r_in_mm)-1)*0.5;
    %     set(p, 'Color', [color_value, 0, 0]);
    %     hold on
    %     p = plot(f,thab{i,2});
    %     set(p, 'Color', [0, 0, color_value])
    % end
    % ylabel('deg'), title(' Theta at base, top view (red) and front view (blue)'), hold off
    % subplot(212)
    % for i = 1 : length(r_in_mm)
    %     p = plot(f,dkp{i,1});
    %     color_value = 0.5 + (i-1)/(length(r_in_mm)-1)*0.5;
    %     set(p, 'Color', [color_value, 0, 0]);
    %     hold on
    %     p = plot(f,dkp{i,2});
    %     set(p, 'Color', [0, 0, color_value])
    % end
    % ylabel('1/mm'), title('Change in kappa, top view (red) and front view (blue)'),
    % xlabel('Time (frames)')

    %%
    % load([num2str(temp_tn+1),'_WST.mat'])
    % tp = [0 4];
    % figure;ws.plot_fitted_whisker_time_projection(0,'k',tp), grid on, hold on;
    % ws.plot_fitted_whisker_time_projection(1,'k',tp)
    % ws.plot_fitted_whisker_ROI_time_projection(0,'r',tp)
    % ws.plot_fitted_whisker_ROI_time_projection(1,'r',tp)
    % ws.plot_mask(0,'g',tp)
    % ws.plot_mask(1,'g',tp)
    % [~,~,y,x,~] = get_theta_kappa_at_roi_point(ws,0,1);
    % plot(x,y,'b.')
    % [~,~,y,x,~] = get_theta_kappa_at_roi_point(ws,1,1);
    % plot(x,y,'b.')

    %%
    % ave_dkp = zeros(length(r_in_mm),2);
    % for i = 1 : length(r_in_mm)
    %     ave_dkp(i,1) = mean(dkp{i,1});
    %     ave_dkp(i,2) = mean(dkp{i,2});
    % end
    % figure, plot(r_in_mm,ave_dkp(:,1),'r-'), hold on, plot(r_in_mm,ave_dkp(:,2),'b-'), xlabel('arc length from the mask (mm)'), ylabel('average delta_kappa (1/mm)') 
    % [~, ind] = min(ave_dkp);
    % max_dkp = [r_in_mm(ind(1)), r_in_mm(ind(2))];
    % disp([num2str(r_in_mm(ind(1))), ' mm (top-view), ', num2str(r_in_mm(ind(2))), ' mm (front-view)'])

    %%
    var_dkp_0 = zeros(lenght(trialNums),length(r_in_mm));
    var_dkp_1 = zeros(lenght(trialNums),length(r_in_mm));
    for i = 1 : length(trialNums)
        for j = 1 : length(r_in_mm)
            var_dkp_0(i,j) = var(dkp_0{i,j});
            var_dkp_1(i,j) = var(dkp_1{i,j});
        end
    end
    % figure, plot(r_in_mm,var_dkp(:,1),'r-'), hold on, plot(r_in_mm,var_dkp(:,2),'b-'), xlabel('arc length from the mask (mm)'), ylabel('maximum variance in delta_kappa (1/mm^2)') 
%     [~, ind] = max(var_dkp);
%     maxvar_dkp = [r_in_mm(ind(1)), r_in_mm(ind(2))];
    % disp([num2str(r_in_mm(ind(1))), ' mm (top-view), ', num2str(r_in_mm(ind(2))), ' mm (front-view)'])
    save('test170329.mat','thab*','dkp*','var_dkp*')
end
toc
