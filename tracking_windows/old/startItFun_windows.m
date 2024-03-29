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
function [conv_time, track_time, copy_time, nfiles] = startItFun_windows(startDir, endDir, pixDen, whiskerNumber)
options.Resize ='on';
options.WindowStyle='normal';
options.Interpreter='none';

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
mp4_list = dir('*.mp4');
seq_list = dir('*.seq');
nfiles = length(seq_list);
if length(mp4_list) < length(seq_list)
    disp('CONVERTING TO .MP4')    
    mp4_converter_parallel_windows;
    system('dir')
    disp('Finished converting')
    
    mp4files = dir('*.mp4');
    if length(mp4files) ~= nfiles
        error('mp4 conversion error (mistmatch in file number)')
    elseif        
        for i = 1 : length(mp4files)
            if strcmp(mp4files(i).name
            
else
    disp('Already converted')
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

if exist('default.parameters','file')
    whiskerTrackerParfor_JK_windows(pixDen, whiskerNumber) % Uses 'classify' for multiple whisker tracking. 
else
    try
        system('copy C:\Users\shires\Documents\GitHub\jkWhisker\default.parameters startDir')
    catch
        error('No default.parameters')
    end
end
disp('Finished tracking')
track_time = toc;

tic;
system(['copy ', startDir, '\*.mp4 ', endDir]) 
system(['copy ', startDir, '\*.whiskers ', endDir]) 
system(['copy ', startDir, '\*.measurements ', endDir]) 
system(['copy ', startDir, '\default.parameters ', endDir]) 
system(['copy ', startDir, '\*.detectorbank ', endDir]) 
copy_time = toc;

% Section 4: Finish program and display time statistics
% Indicates end of program and displays start time, end time, and elapsed
% time 
disp('FINISHED PROGRAM') 
endTimeStat = datestr(now);
disp(['Program started: ' startTimeStat])
disp(['Program ended: ' endTimeStat]) 
end 