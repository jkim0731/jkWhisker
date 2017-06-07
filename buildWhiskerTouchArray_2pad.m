% For 2pad (2-port angle distance), WT, WST, and WL is generated during prebuild_JK.m
% build WL again using information from touch_hyperplane.
tic
touch_boundary_thickness = 2; % include # of pixels outer space from the hp_peak boundary. Default = 2.

behavior_base_dir = 'Z:\Data\2p\soloData\';
whisker_base_dir = 'Z:\Data\Video\JK\';
% mice = {'AH0648','AH0650','AH0651','AH0652','AH0653'};
mice = {'AH0650'};
sessionNum = [1,2,4:6,8:10,19:21];
for mi = 1 : size(mice,2)
    mouseName = mice{mi};
    behavior_d = [behavior_base_dir mouseName '\'];
    try
        load([behavior_d 'behavior.mat']) % loading b of the mouse (all the sessions)
    catch
        error(['No behavior.mat file in ' mouseName ' ' sessionName])
    end
    
    for si = 1 : length(sessionNum)
        sessionName = sprintf('S%02d',sessionNum(si));        
        whisker_d = [whisker_base_dir mouseName sessionName '\'];
        
        try
            load([whisker_d mouseName sessionName '_post.mat'])
        catch
            error(['No _post.mat file in ' mouseName ' ' sessionName])
        end
                
        try 
            load([whisker_d mouseName sessionName '_touch_hp.mat']);
        catch
            error(['No touch_hp.mat file in ' mouseName ' ' sessionName]);
        end

        b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
        b_session = b{b_ind};        
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

            Whisker.makeAllDirectory_WhiskerTrialLite_2pad(whisker_d,'include_files',includef,'r_in_mm',2,'calc_forces',false,'behavior',b_session,'touch_hp',touch_hp,'hp_peaks',hp_peaks,'touch_boundary_thickness',touch_boundary_thickness,'trial_types',trial_types);
        end
    end
end
toc