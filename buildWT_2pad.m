function buildWT_2pad(mouseName, sessionName, d, bSession, videoFreq, ppm, barRadius)
curr_d = pwd;
whisker_d = [d, mouseName, sessionName, filesep];
cd(whisker_d)

tsflist = dir('*_timestamp.mat');
if ~isempty(tsflist) && length(tsflist) > 50 % at least 50 trials with time stamps
    xmlflist = dir('*.xml');
    if ~isempty(xmlflist)
        xmltn = zeros(length(xmlflist),1);
        for i = 1 : length(xmlflist)
            xmltn(i) = str2double(strtok(xmlflist(i).name, '.'));
        end
    else
        xmltn = [];
    end
    meanITI = zeros(length(tsflist),1);
    for i = 1 : length(tsflist)
        tn = str2double(strtok(tsflist(i).name, '_'));
        if ~ismember(tn, xmltn)
            load(tsflist(i).name) % loading tsSec, tsMilli, and tsMicro
            meanITI(i) = mean(diff(tsSec) + diff(tsMilli)/1000 + diff(tsMicro)/1000000);
        end
    end
    meanITI = meanITI(find(meanITI));
    videoFreq = 1/mean(meanITI);
end

load([whisker_d, mouseName, sessionName, '_post.mat']) % for maskx and masky. Ignoring includef from this one. It will be re-defined.
if size(maskx{1},1) > size(maskx{1},2)
    maskx{1} = maskx{1}'; masky{1} = masky{1}';
    maskx{2} = maskx{2}'; masky{2} = masky{2}';
end
filelist=dir([whisker_d '*.measurements']);
dirTrialNums=zeros(1,size(filelist,1));

% %%
% Assign the trial numbers to existing .measurements files in the directory
% NOTE : This assumes that the .measurements files have leading numbers
% corresponding to trial number in string positions 1:end-13 of the file
% name. These index numbers may need to be changed to match up to the
% numerical code of the trial number.  (2016/09/05 JK)

if contains(sessionName, 'spont')
    includef = cell(size(filelist,1),1);
    for i = 1 : length(includef)
        includef{i} = filelist(i).name(1:end-13);
    end
    trialNums = [];
else
    for i=1:length(filelist)
        dirTrialNums(i)=str2double(filelist(i).name(1:end-13)); % extract out the trial number from each measurements file present in directory
    end
    trialNums = sort(dirTrialNums);
    trialNums = trialNums(~isnan(trialNums));
    if ~isempty(bSession) % try only the ones with behavior session
        trialNums = intersect(trialNums,bSession.trialNums); % try only the ones with behavior trials
        %             
        %
%             tinds = (cellfun(@(x) (x.servoAngle == 90) * x.trialNum, bSession.trials));
%             trialNums = tinds(find(tinds));
        %             
        %             
    end
    includef=cell(size(trialNums,1),1);
    for i = 1: length(trialNums)
        includef{i} = num2str(trialNums(i));
    end
end
%%
if ~isempty(bSession)
    Whisker.makeAllDirectory_WhiskerTrial_2pad(whisker_d,[0 1],'mask', {[maskx{1};masky{1}],[maskx{2};masky{2}]},...
        'trial_nums',trialNums,'include_files',includef,...
        'barRadius',barRadius,'faceSideInImage', 'bottom', 'framePeriodInSec',videoFreq,...
        'imagePixelDimsXY',[width height],'pxPerMm',ppm,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward', 'behavior', bSession);
else
    Whisker.makeAllDirectory_WhiskerTrial_2pad(whisker_d,[0 1],'mask', {[maskx{1};masky{1}],[maskx{2};masky{2}]},...
        'trial_nums',trialNums,'include_files',includef,...
        'barRadius',barRadius,'faceSideInImage', 'bottom', 'framePeriodInSec',videoFreq,...
        'imagePixelDimsXY',[width height],'pxPerMm',ppm,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward');
end
cd(curr_d)

end