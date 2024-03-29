%% Setup whisker array builder 
mouseName = 'AH0653'
sessionName ='S03'
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

v = VideoReader([includef{1} '.mp4']);
vv = read(v,1);
vheight = size(vv,1);
vwidth = size(vv,2);
%%
vavg = zeros(vheight,vwidth);

video_ind = linspace(1,length(includef),10);
for j = 1 : 10    
    v = VideoReader([includef{floor(video_ind(j))} '.mp4']);
    for i = 1 : v.NumberOfFrames
        vtemp = read(v,i);    
        vtemp = double(vtemp(:,:,1));
        vavg = vavg + vtemp/v.NumberOfFrames; % average
%         vavg = min(vavg,vtemp); % minimum
    end
end

%% Mask
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


%% save mask
maskfn = [mouseName sessionName 'mask.mat'];
save(maskfn,'maskx','masky', 'includef')