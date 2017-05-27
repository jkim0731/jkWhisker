

%% basic information
% mice = {'AH0648', 'AH0650', 'AH0651', 'AH0652', 'AH0653'};
mice = {'AH0650'};
videoloc = 'JK';
d = (['Z:\Data\Video\' videoloc filesep]);

% comment out when doing for all of the sessions in the mouse directory
% sessions_exc = {{'S02'}, {'S02'}, {'S02'}, {'S04'}, {'S03'}}; 
sessions = {{'S19','S20','S21'}}; 

%% Define follicle points and masks

% comment out when doing for selected sessions in each mouse directory
% for i = 1 : size(mice,2)
%     cd(d)
%     sn = dir([mice{i},'*']);
%     for j = 1 : length(sn)
% %         if sn(j).isdir
%         if sn(j).isdir && ~strcmp(sn(j).name,strcat(mice{i},sessions_exc{i}))
%             [mouseName, sessionName] = strtok(sn(j).name,'S');            
%             follicle_n_mask(mouseName,sessionName,videoloc)
%         end
%         close all
%     end
% end


% comment out when doing for all of the sessions in the mouse directory
for i = 1 : size(mice,2)
    cd(d)
    sn = sessions{i};
    for j = 1 : length(sn)
        mouseName = mice{i};
        sessionName = sn{j};
        follicle_n_mask(mouseName,sessionName,videoloc)
        close all
    end
end


%% Re-measure .measurements file before building whisker arrays 
% comment out when doing for selected sessions in each mouse directory
% for i = 1 : size(mice,2)
%     cd(d)
%     sn = dir([mice{i},'*']);
%     for j = 1 : length(sn)
% %         if sn(j).isdir
%         if sn(j).isdir && ~strcmp(sn(j).name,strcat(mice{i},sessions_exc{i}))
%             [mouseName, sessionName] = strtok(sn(j).name,'S');
%             postmeasurements(mouseName,sessionName,videoloc)
%         end
%     end
% end

% comment out when doing for all of the sessions in the mouse directory
for i = 1 : size(mice,2)
    sn = sessions{i};
    for j = 1 : length(sn)
        mouseName = mice{i};
        sessionName = sn{j};
        postmeasurements(mouseName,sessionName,videoloc)
    end
end


