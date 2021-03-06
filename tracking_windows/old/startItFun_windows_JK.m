%% This command manages the overall pipeline of the whisker tracking process, from untracked .seq to tracked .mp4
%
% Created by J. Sy, 19 November 2014
% This variation made by J. Sy, 11 Dece,ber 2014, to allow for queueing of
% directories
%
% Dependencies: mp4_converterJS3.py, whiskerTrackerParforLinux.m, Janelia Farm Whisker Tracking package 
%
%% Section 1: Directories
% We've got to set a few directories straight here. For usability, input
% commands will be used
function [conv_time, track_time, copy_time, nfiles] = startItFun_windows_JK(startDir, endDir)
% options.Resize ='on';
% options.WindowStyle='normal';
% options.Interpreter='none';

% trackingInfo = inputdlg({'Enter directory where your .seq files are found',...
%     'Enter directory where you would like to place your tracked .mp4s (Should be on NAS mounted under /mnt)',...
%     'Pixel density:', 'Number of Whiskers'},'Please input information', [1 100; 1 100; 1 8; 1 8], {'','','0.016',''},options); %Displays input prompt for directory info

% trackingInfo{1} = Directory files are found 
% trackingInfo{2} = output directory
% trackingInfo{3} = mm per pixel. Defaults to 0.016
% trackingInfo{4} = Number of whiskers on mice, used to determine which
% version of the whisker tracker to use 

% seqDir = input('Please give the directory of the seq files to be converted: ','s');
% trackingDir = input('Please give the directory where you would like to move the .mp4 vidoes to be tracked \nNote: This should be on the NAS: ','s');
%disp('Thank you, no further inputs required') 
startTimeStat = datestr(now); 
tic;

% Section 2: Convert .seq files to .mp4 files
% Uses python script mp4_converterJS3.py, originally written by
% S. Peron, to convert files
cd(startDir);
if exist('.seq','file')
    system('del .seq')
end
mp4_list = dir('*.mp4');
seq_list = dir('*.seq');
nfiles = length(seq_list);
if length(mp4_list) < length(seq_list)
    disp('CONVERTING TO .MP4')    
    mp4_converter_parallel_windows_JK;
    system('dir')
    disp('Finished converting')        
else
    disp('Already converted')
end
mp4files = dir('*.mp4');
if length(mp4files) ~= nfiles
        error('mp4 conversion error: mistmatch in file number')
else        
    for i = 1 : length(mp4files)
        if strcmp(strtok(mp4files(i).name,'.'), strtok(seq_list(i).name,'.'))
        else
            error('mp4 conversion error: file name mismatch in %s', seq_list(i).name)
        end        
    end
end
conv_time = toc;
% Section 3: Upload data to NAS and track the .mp4 files
% Uses 'whiskerTrackerParforLinux.m' which in turn utilizes scripts
% written by Janelia Farm and included in the Whisker Tracking package

% system (['cp ' trackingInfo{1} '/*.mp4 /home/hireslab/WhiskerVideos/Transit']) %Move to temporary directory so we can change permissions below
% system (['chmod ugo+rwx /home/hireslab/WhiskerVideos/Transit/*.mp4']) %Change permissions to avoid write errors later in code
tic;
disp('STARTING WHISKER TRACKING')
disp('Note: This will take some time') 

if ~strcmp(startDir(end),filesep)
    startDir = [startDir, filesep];
end
copy_answer = system(['copy C:\Users\shires\Documents\GitHub\jkWhisker\default.parameters ', startDir]);
if copy_answer ~= 0
    error('No default.parameters')
end
whisker_tracker_true_parallel_JK

disp('Finished tracking')
track_time = toc;

tic;
nf_mp4 = length(dir('*.mp4'));
nf_whiskers = length(dir('*.whiskers'));
nf_measurements = length(dir('*.measurements'));
if nf_mp4 ~= nf_whiskers || nf_mp4 ~= nf_measurements
    disp('Number of files do not match')
    [diff1, diff2] = when_number_of_files_dont_match;
    if ~isempty(diff1) % not traced
        %% (1) TRACE: Uses Janelia Farm's whisker tracking software to track all whiskers in a directory         
        parfor n=1:length(diff1)            
            system(['trace ' num2str(diff1(n)) '.mp4 ' num2str(diff1(n))])
            disp([num2str(diff1(n)) '.mp4 has been traced'])
        end
        diff2 = union(diff1,diff2);
    end
    if ~isempty(diff2) % not measured
        parfor n=1:length(diff2)            
            system(['measure ' '--face ' 'bottom ' num2str(diff2(n)) '.whiskers ' num2str(diff2(n)) '.measurements']);
            disp([num2str(diff2(n)) '.whiskers has been measured'])
        end
    end    
end

%%
system(['copy ', startDir, '*.mp4 ', endDir]) 
system(['copy ', startDir, '*.whiskers ', endDir])
system(['copy ', startDir, '*.measurements ', endDir]) 
system(['copy ', startDir, 'default.parameters ', endDir]) 
system(['copy ', startDir, '*.detectorbank ', endDir]) 
system(['copy ', startDir, '*.xml ', endDir]) 
copy_time = toc;

system(['del ', startDir, '*.mp4']) 
system(['del ', startDir, '*.whiskers'])
system(['del ', startDir, '*.measurements']) 
system(['del ', startDir, 'default.parameters']) 
system(['del ', startDir, '*.detectorbank']) 

% Section 4: Finish program and display time statistics
% Indicates end of program and displays start time, end time, and elapsed
% time 
disp('FINISHED PROGRAM') 
endTimeStat = datestr(now);
disp(['Program started: ' startTimeStat])
disp(['Program ended: ' endTimeStat]) 
end 