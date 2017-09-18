%% whiskerQueue- queues up different directories for tracking assignments 
%pause on 
tic
whiskerNumber = '2';
pixDen = '0.056'; % Telecentric lens of tpm two-view, no binning

startDir = 'Y:\JK_temp\whisker\JK017S25';
endDir = 'Y:\JK_temp\whisker\tracked\JK017S25';
system(['mkdir ', endDir])
[conv_time, track_time, copy_time, nfiles] = startItFun_windows(startDir, endDir, pixDen, whiskerNumber);
toc
% startDir = '/media/hireslab/AHHD_006/Data/AH0167/150705/Camera1';
% endDir = '/mnt/Data/Video/AHHD_006/AH0167/150705';
% startItFun(startDir, endDir, pixDen, whiskerNumber)

% startDir = '/media/hireslab/AHHD_006/Data/AH0167/150706/Camera1';
% endDir = '/mnt/Data/Video/AHHD_006/AH0167/150706';
% startItFun(startDir, endDir, pixDen, whiskerNumber)