%% whiskerQueue- queues up different directories for tracking assignments 
%pause on 
% 
% whiskerNumber = '2';
% pixDen = '0.056'; % Telecentric lens of tpm two-view, no binning
delete(gcp('nocreate'))
cd('Y:\JK_temp\whisker')
dirlist = dir('JK*');
conv_time = zeros(length(dirlist)-1,1);
track_time = zeros(length(dirlist)-1,1);
copy_time = zeros(length(dirlist)-1,1);
nfiles = zeros(length(dirlist)-1,1);
for i = 3 : length(dirlist)        
    delete(gcp('nocreate'))
    tic
    if dirlist(i).isdir
        startDir = ['Y:\JK_temp\whisker\', dirlist(i).name];
        endDir = ['Y:\JK_temp\whisker\tracked\', dirlist(i).name];
        system(['mkdir ', endDir])        
        [conv_time(i), track_time(i), copy_time(i), nfiles(i)] = startItFun_windows_JK(startDir, endDir);
    end
    session_time = toc;
end

% %%
% startDir = 'Y:\JK_temp\whisker\JK018S20';
% endDir = 'Y:\JK_temp\whisker\tracked\JK018S20';
% system(['mkdir ', endDir])   
% cd(startDir)
% whiskerTrackerParfor_JK_windows(pixDen, whiskerNumber)
% 
% system(['copy ', startDir, '\*.mp4 ', endDir]) 
% system(['copy ', startDir, '\*.whiskers ', endDir]) 
% system(['copy ', startDir, '\*.measurements ', endDir]) 
% system(['copy ', startDir, '\default.parameters ', endDir]) 
% system(['copy ', startDir, '\*.detectorbank ', endDir]) 