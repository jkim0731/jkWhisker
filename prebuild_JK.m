

%% basic information
% mice = {'AH0648', 'AH0650', 'AH0651', 'AH0652', 'AH0653'};
mice = {'AH0652', 'AH0653'};
videoloc = 'JK';
d = (['Z:\Data\Video\' videoloc filesep]);

% comment out when doing for all of the sessions in the mouse directory
% sessions_exc = {{'S02'}, {'S02'}, {'S02'}, {'S04'}, {'S03'}}; 
sessions_exc = {{'S04'}, {'S03'}}; 

%% Define follicle points and masks

% comment out when doing for selected sessions in each mouse directory
for i = 1 : size(mice,2)
    cd(d)
    sn = dir([mice{i},'*']);
    for j = 1 : length(sn)
%         if sn(j).isdir
        if sn(j).isdir && ~strcmp(sn(j).name,strcat(mice{i},sessions_exc{i}))
            [mouseName, sessionName] = strtok(sn(j).name,'S');            
            follicle_n_mask(mouseName,sessionName,videoloc)
        end
        close all
    end
end


% comment out when doing for all of the sessions in the mouse directory
% for i = 1 : size(mice,2)
%     cd(d)
%     sn = sessions{i};
%     for j = 1 : length(sn)
%         mouseName = mice{i};
%         sessionName = sn{j};
%         follicle_n_mask(mouseName,sessionName,videoloc)
%         close all
%     end
% end


%% Re-measure .measurements file before building whisker arrays 
% comment out when doing for selected sessions in each mouse directory
for i = 1 : size(mice,2)
    cd(d)
    sn = dir([mice{i},'*']);
    for j = 1 : length(sn)
%         if sn(j).isdir
        if sn(j).isdir && ~strcmp(sn(j).name,strcat(mice{i},sessions_exc{i}))
            [mouseName, sessionName] = strtok(sn(j).name,'S');
            postmeasurements(mouseName,sessionName,videoloc)
        end
    end
end

% comment out when doing for all of the sessions in the mouse directory
% for i = 1 : size(mice,2)
%     sn = sessions{i};
%     for j = 1 : length(sn)
%         mouseName = mice{i};
%         sessionName = sn{j};
%         postmeasurements(mouseName,sessionName,videoloc)
%     end
% end

%%



%%
trialNum = size(dkp_0,1);
baseline_frames = 1:450;
figure, subplot(121),
for i = 1 : trialNum
    if length(wl.trials{i}.time{1}) ~= length(thab_0{i,1}) || length(thab_0{i,1}) ~= length(dkp_0{i,1})
        disp(['trial Num ' num2str(i) ' does not match in length'])
    else
        plot(thab_0{i,1}(baseline_frames),dkp_0{i,1}(baseline_frames),'k.'), xlabel('Theta at base'), ylabel('kappa'), title('Top-view theta vs kappa during free whisking'), hold on
    end
end
