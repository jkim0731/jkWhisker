function postmeasurements(mouseName,sessionName,videoloc,ppm,varargin)

%% Setup whisker array builder 
% mouseName = 'AH0653'
% sessionName ='S03'
% videoloc = 'JK'
% flag_skip = 'skip'
% optional = 'Spont'
if nargin > 4
    flag_skip = varargin{1};
end
if nargin > 5
    optional = varargin{2};
elseif nargin > 6
    error('Too many input arguments')
end

if exist('optional','var')
    d = ([videoloc filesep mouseName sessionName filesep optional filesep])
else
    d = ([videoloc filesep mouseName sessionName filesep])
end
% load(['Z:\Users\Jon\DATA\BehaviorArrays\solo_' mouseName '_' sessionName '.mat'])

currD = pwd;
try
    cd(d)
catch
    disp(['Directory ', mouseName, sessionName, ' does not exist. Skipping'])
    return
end


savefn = [mouseName sessionName '_post.mat'];
if exist(savefn,'file')
    if exist('flag_skip', 'var') && strcmp(flag_skip, 'skip')
        disp(['Already post-measured. Skip ', mouseName, ' ', sessionName])
        return
    else
        disp(['Over-writing already present post-measurement in ', mouseName, ' ', sessionName])
    end    
end
%%
excludef = jkmeasurements_dir(ppm);
%%
filelist=dir([d '*.measurements']);

dirTrialNums=zeros(1,size(filelist,1));
%%
% Assign the trial numbers to existing .measurements files in the directory
% NOTE : This assumes that the .measurements files have leading numbers
% corresponding to trial number in string positions 1:end-13 of the file
% name. These index numbers may need to be changed to match up to the
% numerical code of the trial number.  (2016/09/05 JK)

for i=1:length(filelist)
    dirTrialNums(i)=str2double(filelist(i).name(1:end-13)); % extract out the trial number from each measurements file present in directory
end
dirTrialNums = setdiff(dirTrialNums,excludef);
trialNums = sort(dirTrialNums);
trialNums = trialNums(~isnan(trialNums));

includef=cell(size(trialNums,1),1);
for i = 1: length(trialNums)
    includef{i} = num2str(trialNums(i));
end

%% save information from postmeasurements
savefn = [mouseName sessionName '_post.mat'];
loadfn1 = dir('*follicle_n_mask.mat');
loadfn2 = dir('remeasure*.mat');
load(loadfn1.name)
load(loadfn2.name)
save(savefn,'maskx','masky', 'includef', 'errorlist', 'width', 'height')
cd(currD)