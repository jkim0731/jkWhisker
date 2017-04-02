function postmeasurements(mouseName,sessionName,videoloc,varargin)

%% Setup whisker array builder 
% mouseName = 'AH0653'
% sessionName ='S03'
% videoloc = 'JK'
% optional = 'Spont'
if nargin > 3
    optional = varargin{4};
elseif nargin > 4
    error('Too many input arguments')
end

if exist('optional','var')
    d = (['Z:\Data\Video\' videoloc filesep mouseName sessionName filesep optional filesep])
else
    d = (['Z:\Data\Video\' videoloc filesep mouseName sessionName filesep])
end
% load(['Z:\Users\Jon\DATA\BehaviorArrays\solo_' mouseName '_' sessionName '.mat'])

cd(d)

%%
excludef = jkmeasurements_dir();
%%
filelist=dir([d '*.measurements']);

dirTrialNums=zeros(1,size(filelist,1));
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

%% save mask
maskfn = [mouseName sessionName '_post.mat'];
save(maskfn,'maskx','masky', 'includef')