curr_dir = pwd;
video_dir = 'Z:\Data\Video\JK';
cd(video_dir)
dir_list = dir;
cd(curr_dir)
% dir_list = dir_list(
% mice = {'AH0648', 'AH0650', 'AH0651', 'AH0652', 'AH0653'};
% for i = 1 : size(mice,1)
%     
% end
%%
st_ind = []; % second_trial index
for i = 1 : size(dir_list,1)
    try
        if str2double(dir_list(i).name(4:6)) > 647 && strcmp(dir_list(i).name(end-1:end),'_2')
            st_ind = [st_ind;i];
        end
    catch
    end
end

%%
maxtn = [318,30];

%%
merge_dir_ind = st_ind - 1;
% for i = 1 : size(st_ind,1)
for i = 2
    cd([video_dir,filesep,dir_list(st_ind(i)).name])
    source_mp4 = ls('*.mp4');     
    for j = 1 : size(source_mp4,1)        
        stn = str2double(strtok(source_mp4(j,:),'.')); % source trial number
        FileRename([video_dir,filesep,dir_list(st_ind(i)).name,filesep,   source_mp4(j,:)], ...
            [video_dir,filesep,dir_list(merge_dir_ind(i)).name,filesep,   num2str(maxtn(i)+stn-1), '.mp4'])        
    end
    source_whiskers = ls('*.whiskers');     
    for j = 1 : size(source_whiskers,1)
        stn = str2double(strtok(source_whiskers(j,:),'.')); % source trial number
        FileRename([video_dir,filesep,dir_list(st_ind(i)).name,filesep,   source_whiskers(j,:)], ...
            [video_dir,filesep,dir_list(merge_dir_ind(i)).name,filesep,   num2str(maxtn(i)+stn-1), '.whiskers'])        
    end
    source_measurements = ls('*.measurements');     
    for j = 1 : size(source_measurements,1)
        stn = str2double(strtok(source_measurements(j,:),'.')); % source trial number
        FileRename([video_dir,filesep,dir_list(st_ind(i)).name,filesep,   source_measurements(j,:)], ...
            [video_dir,filesep,dir_list(merge_dir_ind(i)).name,filesep,   num2str(maxtn(i)+stn-1), '.measurements'])        
    end
    cd('..')
end

%%
% k = 9;
% stn = str2double(strtok(source_mp4(k,:),'.mp4'))
% maxtn(i)+stn-1
% stnright = str2double(strtok(source_mp4(k,:),'.'))
% maxtn(i)+stnright-1