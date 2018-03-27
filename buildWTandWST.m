function buildWTandWST(mouseName, sessionName, d, b_session, ppm)
    curr_d = pwd;
    whisker_d = [d, mouseName, sessionName, filesep];   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    cd(whisker_d)
    delete *_WT.mat
    delete *_WST.mat
    delete *_errorWST.mat
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    load([whisker_d, mouseName, sessionName, '_post.mat']) % for maskx and masky. Ignoring includef from this one. It will be re-defined.
    
    if ~isempty(b_session) % try only the ones with behavior session

        filelist=dir([whisker_d '*.measurements']);
        dirTrialNums=zeros(1,size(filelist,1));

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
        trialNums = intersect(trialNums,b_session.trialNums); % try only the ones with behavior trials

        includef=cell(size(trialNums,1),1);
        for i = 1: length(trialNums)
            includef{i} = num2str(trialNums(i));
        end
    end
%%
if size(maskx{1},1) > size(maskx{1},2)
    Whisker.makeAllDirectory_WhiskerTrial(whisker_d,[0 1],'mask', {[maskx{1}';masky{1}'],[maskx{2}';masky{2}']},...
        'trial_nums',trialNums,'include_files',includef,...
        'barRadius',6,'faceSideInImage', 'bottom', 'framePeriodInSec',0.003225806451613,...
        'imagePixelDimsXY',[width height],'pxPerMm',ppm,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward');
else
    Whisker.makeAllDirectory_WhiskerTrial(whisker_d,[0 1],'mask', {[maskx{1};masky{1}],[maskx{2};masky{2}]},...
        'trial_nums',trialNums,'include_files',includef,...
        'barRadius',6,'faceSideInImage', 'bottom', 'framePeriodInSec',0.003225806451613,...
        'imagePixelDimsXY',[width height],'pxPerMm',ppm,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','rightward');
end
Whisker.makeAllDirectory_WhiskerSignalTrial_2pad(whisker_d,'include_files',includef,'polyRoiInPix',[ppm 6*ppm]);
    
%                 Whisker.makeAllDirectory_WhiskerTrialLite_2pad(whisker_d,'include_files',includef,'r_in_mm',3,'calc_forces',false,'behavior',b_session);

cd(curr_d)

end