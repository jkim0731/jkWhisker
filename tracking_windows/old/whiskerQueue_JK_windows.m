%% whiskerQueue- queues up different directories for tracking assignments 
%pause on 
% 
% whiskerNumber = '2';
% pixDen = '0.056'; % Telecentric lens of tpm two-view, no binning
delete(gcp('nocreate'))

%%
cd('D:\')
dirlist = dir('JK*');
conv_time = zeros(length(dirlist)-1,1);
track_time = zeros(length(dirlist)-1,1);
copy_time = zeros(length(dirlist)-1,1);
nfiles = zeros(length(dirlist)-1,1);

%%
for i = 1 : length(dirlist)        
    delete(gcp('nocreate'))
    tic
    if dirlist(i).isdir
        if isempty(strfind(dirlist(i).name,'spont'))
            startDir = ['D:\', dirlist(i).name];
            endDir = ['D:\tracked\', dirlist(i).name];
            system(['mkdir ', endDir])        
            [conv_time(i), track_time(i), copy_time(i), nfiles(i)] = startItFun_windows_JK(startDir, endDir);
        end
    end
    session_time = toc;
end

%%
cd('E:\JK036spont')
dirlist = dir('JK*');
conv_time = zeros(length(dirlist)-1,1);
track_time = zeros(length(dirlist)-1,1);
copy_time = zeros(length(dirlist)-1,1);
nfiles = zeros(length(dirlist)-1,1);

%%
for i = 1 : length(dirlist)        
    delete(gcp('nocreate'))
    tic
    if dirlist(i).isdir
        if isempty(strfind(dirlist(i).name,'spont'))
            startDir = ['E:\JK036spont\', dirlist(i).name];
            endDir = ['E:\JK036spont\tracked\', dirlist(i).name];
            system(['mkdir ', endDir])        
            [conv_time(i), track_time(i), copy_time(i), nfiles(i)] = startItFun_windows_JK(startDir, endDir);
        end
    end
    session_time = toc;
end
