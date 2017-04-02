%% Setup whisker array builder 
mouseName = 'AH0648'
sessionName ='S02'
videoloc = 'JK'
% optional = 'Spont'

if exist('optional','var')
    d = (['Z:\Data\Video\' videoloc filesep mouseName sessionName filesep optional filesep])
else
    d = (['Z:\Data\Video\' videoloc filesep mouseName sessionName filesep])
end
% load(['Z:\Users\Jon\DATA\BehaviorArrays\solo_' mouseName '_' sessionName '.mat'])

cd(d)

%%
excludef = jkmeasurements_dir(); % just the numbers 03/31/2017 JK
%%
filelist=dir([d '*.measurements']);

dirTrialNums=zeros(1,size(filelist,1));
% trialNums=[];  % enter which trial nums to process 

%%
% Assign the trial numbers to existing .measurements files in the directory
% NOTE : This assumes that the .measurements files have leading numbers
% corresponding to trial number in string positions 1:end-13 of the file
% name. These index numbers may need to be changed to match up to the
% numerical code of the trial number.  (2016/09/05 JK)

for i=1:length(filelist);
    dirTrialNums(i)=str2double(filelist(i).name(1:end-13)); % extract out the trial number from each measurements file present in directory
end
dirTrialNums = setdiff(dirTrialNums,excludef);
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
vavg = zeros(vheight,vwidth);

video_ind = linspace(includef{1}, includef{end},10);
for j = 1 : 10    
    v = VideoReader([includef{floor(video_ind(j))} '.mp4']);
    for i = 1 : v.NumberOfFrames
        vtemp = read(v,i);    
        vtemp = double(vtemp(:,:,1));
        vavg = vavg + vtemp/v.NumberOfFrames; % average
%         vavg = max(vavg,vtemp); % maximum
    end
end
%%

% Assign the trial numbers to existing .measurements files in the directory
% NOTE : This DOES NOT assume that the .measurements files have leading numbers
% corresponding to trial number in string positions 1:end-13 of the file
% name. It can be any form. For Spont analysis (2016/11/08 JK)

% for i=1:length(filelist);
%     dirTrialNums(i)=str2double(filelist(i).name(1:end-13)); % extract out the trial number from each measurements file present in directory
% end
% trialNums = sort(dirTrialNums);
% trialNums = trialNums(~isnan(trialNums));
% 
% includef=cell(size(trialNums,1),1);
% for i = 1: length(trialNums)
%     includef{i} = num2str(trialNums(i));
% end
% 
% v = VideoReader([includef{1} '.mp4']);
% vv = read(v,1);
% vheight = size(vv,1);
% vwidth = size(vv,2);
% vavg = zeros(vheight,vwidth);
% for i = 1 : v.NumberOfFrames
%     vtemp = read(v,i);    
%     vtemp = double(vtemp(:,:,1));
%     vavg = vavg + vtemp/v.NumberOfFrames;
% end

%% Optional section for cross correlating behavior and video trials, if you didn't pay attention to trial numbers
% vv = nan(max(dirTrialNums),1);
% 
% for i = 1:length(dirTrialNums) 
%     
%     
%     if ~isempty(find(dirTrialNums == i,1)); %create matrix showing trial location
%         fidx = find(dirTrialNums == i); %show index in matrix where trial is 
% 
%         tmp = load([filelist(fidx).name(1:end-13) '.bar']);
%          vv(i) = tmp(1,2);
%     else
%     end
%     
% end
% 
% gngThreshold = nanmean(vv) % for continuous pole postion with equal width go/nogo ranges, mean vv = the transition between go/nogo.
% vv2 = vv >= gngThreshold; % threshold for pole position on 1 side or the other
% vv3 = vv < gngThreshold;
% vvDiff = vv2-vv3;
% figure
% plot([vvDiff],'.')
% hold on
% plot(vv3,'ro')

% %%  Build behavior number vector
% bv = zeros(max(b.trialNums),1);
% 
% bv(b.trialNums) = b.trialTypes*2-1; %1=GO and -1=NOGO and 0=NAN
% [c, lags] = xcorr(bv,vvDiff); %correlation between bv (behavioral GO and NOGO) and vv2 (video GO and NOGO trials)
% [mc, mx] = max(c); %find max c value and max x value
% lag_shift = lags(mx) %find lag in between both bv and vv2
% 
% for i = 1: length(filelist)
%     includef{i} = filelist(i).name(1:end-13);
% end
% 
% % self inputted values since tossed first 132 trials and know lag shift
% % lag_shift = -1
% 
% trialNums = dirTrialNums+lag_shift;  % correct the trial numbers
% 
% hold on
% plot(bv,'go')
% 
% % restrict processing to trials...
% startTrial = 1;
% endTrial = 99999999;
% includef = includef(trialNums >= startTrial & trialNums <=endTrial);
% trialNums =  trialNums(trialNums >= startTrial & trialNums <=endTrial);

%% MASK
figure, imshow(mat2gray(vavg)), axis off, axis image, hold all;
% plot(x1,y1,'r.',x2,y2,'r.')
button = 1;
masknum = str2num(cell2mat(inputdlg({'How many trajectories?','How many points?'},'Trajectories',1,{'2','3'})));
maskx = zeros(masknum(1),masknum(2));
masky = zeros(masknum(1),masknum(2));
for i = 1 : masknum(1)
    for j = 1 : masknum(2)
        [x, y, button] = ginput(1);
        x = round(x); y = round(y);        
        scatter(x,y,'mo');
        maskx(i,j) = x;
        masky(i,j) = y;
        if j > 1
            plot(maskx(i,j-1:j), masky(i,j-1:j))
        end
    end
    
end
hold off;

%% Step 2 - Run with mask! for the first 8 trials and see how it looks

temp_tn = [50];

Whisker.makeAllDirectory_WhiskerTrial(d,[0 1],'mask', {[maskx(1,:);masky(1,:)],[maskx(2,:);masky(2,:)]},...
    'trial_nums',trialNums(temp_tn),'include_files',includef(temp_tn),...
    'barRadius',15.3,'faceSideInImage', 'bottom', 'framePeriodInSec',.0032,...
    'imagePixelDimsXY',[vwidth vheight],'pxPerMm',26.23,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward')

Whisker.makeAllDirectory_WhiskerSignalTrial(d,'include_files',includef(temp_tn),'polyRoiInPix',[40-20 40+20],'follicleExtrapDistInPix',33);
Whisker.makeAllDirectory_WhiskerTrialLiteI(d,'include_files',includef(temp_tn),'r_in_mm',1,'calc_forces',false);
wl = Whisker.WhiskerTrialLiteArray(d);

tid = [0 1]; % Set trajectory ID to view
Whisker.viewdouble_WhiskerTrialLiteArray(wl,tid)


%% check mask - load .WST file first
tp = [0 4];
figure;ws.plot_fitted_whisker_time_projection(0,'k',tp), grid on
hold on; ws.plot_fitted_whisker_time_projection(1,'k',tp)
hold on; ws.plot_fitted_whisker_ROI_time_projection(0,'r',tp)
hold on; ws.plot_fitted_whisker_ROI_time_projection(1,'r',tp)
hold on; ws.plot_mask(0,'g',tp)
hold on; ws.plot_mask(1,'g',tp)

%% Step 3 - run everything
% select matching files
%tmp = cellfun(@(x)str2num(x(15:end)),includef);
%incf_idx = find(tmp>= 12 & tmp <=85);

Whisker.makeAllDirectory_WhiskerTrial(d,[0 1],'mask', {[maskx(1,:);masky(1,:)],[maskx(2,:);masky(2,:)]},...
    'trial_nums',trialNums,'include_files',includef,...
    'barRadius',15.3,'faceSideInImage', 'top', 'framePeriodInSec',.0032,...
    'imagePixelDimsXY',[460 270],'pxPerMm',26.23,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','leftward')

Whisker.makeAllDirectory_WhiskerSignalTrial(d,'include_files',includef,'polyRoiInPix',[99-33 99+33],'follicleExtrapDistInPix',33);
Whisker.makeAllDirectory_WhiskerTrialLiteI(d,'include_files',includef,'r_in_mm',3,'calc_forces',true,'whisker_radius_at_base', 36.5,'whisker_length', 18,'baseline_time_or_kappa_value',0);
wl = Whisker.WhiskerTrialLiteArray(d);
save([d mouseName sessionName '-WTLIA.mat'],'wl');

tid = 0; % Set trajectory ID to view
Whisker.view_WhiskerTrialLiteArray(wl,tid)