%% Setup whisker array builder 
base_angle = 21;

behavior_base_dir = 'Z:\Data\2p\soloData\';
whisker_base_dir = 'Z:\Data\Video\JK\';

sessions = {'AH0648S02','AH0650S02','AH0651S02','AH0652S04','AH0653S03'};
[mouseName, sessionName] = strtok(sessions{1},'S');
% mouse = 'AH0648';
% session = 'S02';

behavior_d = [behavior_base_dir mouseName '\'];
whisker_d = [whisker_base_dir mouseName sessionName '\'];
cd(behavior_d)
load('behavior.mat') % loading b of the mouse (all the sessions)

bb = cellfun(@(x) x.sessionName,b,'UniformOutput',false);
b_ind = find(cellfun(@(x) strcmp(x,sessionName),bb));
b_session = b{b_ind};

cd(whisker_d)

load_fn = [mouseName sessionName '_post.mat'];
load(load_fn); % loading errorlist
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
dirTrialNums = setdiff(dirTrialNums,errorlist);
trialNums = sort(dirTrialNums);
trialNums = trialNums(~isnan(trialNums));

includef=cell(size(trialNums,1),1);
for i = 1: length(trialNums)
    includef{i} = num2str(trialNums(i));
end

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

%%  Build behavior number vector
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

%% Step 2 - Run with mask! for the first 8 trials and see how it looks
% 
% temp_tn = [50];
% 
% Whisker.makeAllDirectory_WhiskerTrial(d,[0 1],'mask', {[maskx(1,:);masky(1,:)],[maskx(2,:);masky(2,:)]},...
%     'trial_nums',trialNums(temp_tn),'include_files',includef(temp_tn),...
%     'barRadius',15.3,'faceSideInImage', 'bottom', 'framePeriodInSec',.0032,...
%     'imagePixelDimsXY',[vwidth vheight],'pxPerMm',26.23,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward')
% 
% % Whisker.makeAllDirectory_WhiskerSignalTrial(d,'include_files',includef(temp_tn),'polyRoiInPix',[40-20 40+20],'follicleExtrapDistInPix',26);
% Whisker.makeAllDirectory_WhiskerSignalTrial(d,'include_files',includef(temp_tn),'polyRoiInPix',[40-20 40+20]);
% Whisker.makeAllDirectory_WhiskerTrialLiteI(d,'include_files',includef(temp_tn),'r_in_mm',1,'calc_forces',false);
% wl = Whisker.WhiskerTrialLiteArray(d);
% 
% tid = [0 1]; % Set trajectory ID to view
% Whisker.viewdouble_WhiskerTrialLiteArray(wl,tid)

%% check mask - load .WST file first
% tp = [0 4];
% figure;ws.plot_fitted_whisker_time_projection(0,'k',tp), grid on
% hold on; ws.plot_fitted_whisker_time_projection(1,'k',tp)
% hold on; ws.plot_fitted_whisker_ROI_time_projection(0,'r',tp)
% hold on; ws.plot_fitted_whisker_ROI_time_projection(1,'r',tp)
% hold on; ws.plot_mask(0,'g',tp)
% hold on; ws.plot_mask(1,'g',tp)

%% Step 3 - run everything
% select matching files
%tmp = cellfun(@(x)str2num(x(15:end)),includef);
%incf_idx = find(tmp>= 12 & tmp <=85);

%% Make whisker-pole touch space for each type of trial, from 10 randomly selected trials (of each type)
% Currently, only dealing with 4 types of trials: 'rc', 'rf', 'lc', 'lf'
% Should make something different for straight pole touch in S00. 
% 2017/04/11 JK
trial_types = {'rc', 'rf', 'lc', 'lf'};
tt_ind = cell(1,length(trial_types));
touch_points = cell(1,length(trial_types));
%%
% for i = 1 : length(trial_types)    
i = 1
    tt_ind{i} = find(cellfun(@(x) strcmp(x.trialType,trial_types{i}),b_session.trials));
    if length(tt_ind{i}) > 10
        tt_ind{i} = sort(randsample(tt_ind{i},10));
    end
    temp_files = cell(length(tt_ind{i}),1);
    for j = 1 : length(tt_ind{i})
        temp_files{j} = num2str(tt_ind{i}(j));
    end
    Whisker.makeAllDirectory_WhiskerSignalTrial(whisker_d,'include_files',temp_files,'polyRoiInPix',[20 80]);
    Whisker.makeAllDirectory_WhiskerTrialLiteI(whisker_d,'include_files',temp_files,'r_in_mm',2,'calc_forces',false,'behavior',b_session);
    wl = Whisker.WhiskerTrialLiteArray(whisker_d,'include_files',temp_files);

    tid = [0 1]; % Set trajectory ID to view
    Whisker.viewdouble_WhiskerTrialLiteArray(wl,tid)
% end
%%
tt_ind{1}
%%
touch_points{1} = cell(1,length(wl.trialNums));
% [138,216,332,344,414,418,498,705,764,787];
touch_points{1}{1} = [592, 615, 649, 675, 689, 707, 722, 1121, 1146, 1163];
touch_points{1}{2} = [1130, 1143, 1153];
touch_points{1}{3} = [682, 772, 780, 856, 891, 905, 919, 944, 1014, 1020, 1032, 1146, 1275, 1344,1345];
touch_points{1}{4} = [596,980,1005,1223,2020,2085,2092,2419, 2434,2850,3188,3221,3340,3496,3508,3526];
touch_points{1}{5} = [1110,1130];
touch_points{1}{6} = [534:538];
touch_points{1}{7} = [1439,1559,1609,1625,2853,2858,2872,2899,2954];
touch_points{1}{8} = [612,740, 746,843,873,911,1095,1134,1190,1205];
touch_points{1}{9} = [589,605,664,679,716,725,732,764,885,917,1292,1536,1574,1592,1677];
% touch_points{1}{10} = [];

%%
for j = 1 : length(wl.trialNums)
    load([num2str(tt_ind{i}(j)),'_WL.mat']);
    
end


%%
Whisker.makeAllDirectory_WhiskerTrial(whisker_d,[0 1],'mask', {[maskx(1,:);masky(1,:)],[maskx(2,:);masky(2,:)]},...
    'trial_nums',trialNums,'include_files',includef,...
    'barRadius',15.3,'faceSideInImage', 'top', 'framePeriodInSec',.0032,...
    'imagePixelDimsXY',[width height],'pxPerMm',26.23,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','leftward')

Whisker.makeAllDirectory_WhiskerSignalTrial(whisker_d,'include_files',includef,'polyRoiInPix',[20 80]);
Whisker.makeAllDirectory_WhiskerTrialLiteI(whisker_d,'include_files',includef,'r_in_mm',3,'calc_forces',true,'whisker_radius_at_base', 36.5,'whisker_length', 18,'baseline_time_or_kappa_value',0);
wl = Whisker.WhiskerTrialLiteArray(whisker_d);
save([d mouseName sessionName '-WTLIA.mat'],'wl');

tid = 0; % Set trajectory ID to view
Whisker.view_WhiskerTrialLiteArray(wl,tid)
