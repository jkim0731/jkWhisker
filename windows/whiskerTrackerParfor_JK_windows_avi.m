%% Used to do whisker tracking on raw video data
%Attempts to process in parallel
%In theory, sections 3 (CLASSIFY) and 3.5 (CLASSIFY-SINGLE-WHISKER) should
%not be in use at the same time, since they do the same thing
% Created by SAH and JS, 2014-08-04

%% (1) TRACE: Uses Janelia Farm's whisker tracking software to track all whiskers in a directory 
function whiskerTrackerParfor_JK_windows_avi(varargin)
if nargin > 1
    pixDen = varargin{1};
    if nargin > 2
        whiskerNum = varargin{2};
    end
end
% cd(['Z:\Data\Video\JON\AH0717\170901'])

% delete(gcp('nocreate')); %turns off all other parallel pool processes
% numCores = feature('numcores'); %identify number of cores available for MATLAB to use
% parpool('local',numCores); %parallel pool using max number of cores available

%%
traces = dir('*.avi'); %Searches only for .avi files, change if using other type (e.g. SEQ)

% for n=1
parfor n=1:length(traces)
    [~, outputFileName] = fileparts(traces(n).name);
%     system(['trace ' traces(n).name ' test'])
    sout = system(['trace ' traces(n).name ' ' outputFileName]);
    trial_ind = 0;
    while (sout~= 0 && trial_ind < 3) % sometimes there can be an error. Try again at least 3 times more
        sout = system(['trace ' traces(n).name ' ' outputFileName]);
        trial_ind = trial_ind + 1;
    end
    if sout == 0
        disp([traces(n).name ' has been traced'])   
    else
        disp(['Failed to trace ' traces(n).name])   
    end
end

%% (2) MEASURE: Generates measurements of traced shapes for later evaluation
measures = dir('*.whiskers');

parfor n=1:length(measures)
    [~, outputFileName] = fileparts(measures(n).name);
    sout = system(['measure ' '--face ' 'bottom ' measures(n).name ' ' outputFileName '.measurements']);
    trial_ind = 0;
    while (sout ~= 0 && trial_ind < 3)
        sout = system(['measure ' '--face ' 'bottom ' measures(n).name ' ' outputFileName '.measurements']);
        trial_ind = trial_ind + 1;
    end
    if sout == 0
        disp([measures(n).name ' has been measured'])
    else
        disp(['Failed to measure ' measures(n).name])
    end
end

%% (2)-1 Error check and re-do the analysis (2017/09/18 JK)
avi_flist = dir('*.avi');
whiskers_flist = dir('*.whiskers');
measure_flist = dir('*.measurements');
if length(measure_flist) < length(avi_flist)
    % 1) re-trace those failed to trace before
    if length(whiskers_flist) < length(avi_flist)
        avi_list = zeros(length(avi_flist),1);
        for i = 1 : length(avi_flist)
            avi_list(i) = str2double(strtok(avi_flist(i).name,'.')); % assume that all filenames are integers
        end
        whiskers_list = zeros(length(whiskers_flist),1);
        for i = 1 : length(whiskers_flist)
            whiskers_list(i) = str2double(strtok(whiskers_flist(i).name,'.'));
        end
        trace_errorlist = setdiff(avi_list,whiskers_list);    
        parfor i = 1 : length(trace_errorlist)
            temp_fname = num2str(trace_errorlist(i));
            sout = system(['trace ' temp_fname '.avi ' temp_fname]);
            trial_ind = 0;
            while (sout ~= 0 && trial_ind < 3)
                sout = system(['trace ' temp_fname '.avi ' temp_fname]);
                trial_ind = trial_ind + 1;
            end
            if sout == 0
                disp([temp_fname '.avi traced successfully'])
            else
                disp(['Failed to trace ' temp_fname])
            end
        end
    end
    % 2) re-measure    
    measure_list = zeros(length(measure_flist),1);
    for i = 1 : length(measure_flist)
        measure_list(i) = str2double(strtok(measure_flist(i).name,'.'));
    end
    measure_errorlist = setdiff(avi_list,measure_list);
    parfor i = 1 : length(measure_errorlist)
        temp_fname = num2str(measure_errorlist(i));
        sout = system(['measure ' '--face ' 'bottom ' temp_fname '.whiskers ' temp_fname '.measurements']);
        trial_ind = 0;
        while (sout ~= 0 && trial_ind < 3)
            sout = system(['measure ' '--face ' 'bottom ' temp_fname '.whiskers ' temp_fname '.measurements']);
            trial_ind = trial_ind + 1;
        end
        if sout == 0
            disp([temp_fname '.whiskers has been measured'])
        else
            disp(['Failed to measure ' temp_fname '.whiskers'])
        end
    end
end
%% 3) for those still remaining not measured, trace again and then  measure them
measure_flist = dir('*.measurements');
if length(measure_flist) < length(avi_flist)
    measure_list = zeros(length(measure_flist),1);
    for i = 1 : length(measure_flist)
        measure_list(i) = str2double(strtok(measure_flist(i).name,'.'));
    end
    measure_errorlist = setdiff(avi_list,measure_list);
    parfor i = 1 : length(measure_errorlist)
        temp_fname = num2str(measure_errorlist(i));
        sout = system(['trace ' temp_fname '.avi ' temp_fname]);
        trial_ind = 0;
        while (sout ~= 0 && trial_ind < 3)
            sout = system(['trace ' temp_fname '.avi ' temp_fname]);
            trial_ind = trial_ind + 1;
        end
        if sout == 0
            disp([temp_fname '.avi traced successfully'])
        else
            disp(['Failed to trace ' temp_fname])
        end
    end
    parfor i = 1 : length(measure_errorlist)
        temp_fname = num2str(measure_errorlist(i));
        sout = system(['measure --face bottom ' temp_fname '.whiskers ' temp_fname '.measurements']);
        trial_ind = 0;
        while (sout ~= 0 && trial_ind < 3)
            sout = system(['measure --face bottom ' temp_fname '.whiskers ' temp_fname '.measurements']);
            trial_ind = trial_ind + 1;
        end
        if sout == 0
            disp([temp_fname '.whiskers measured successfully'])
        else
            disp(['Failed to measure ' temp_fname])
        end
    end
end
%% (3) CLASSIFY: Helps refine tracing to more accurately determine which shapes are whiskers
% Use for multiple whiskers
classes = dir('*.measurements');

parfor n=1:length(classes)
    [~, outputFileName] = fileparts(classes(n).name);
    system(['classify ' classes(n).name ' ' outputFileName '.measurements ' 'bottom ' '--px2mm ' pixDen ' -n ' whiskerNum]);
    display([classes(n).name ' has been classified'])
end


    
%% (3.5) CLASSIFY-SINGLE-WHISKER: A variation of the above code designed for one whisker
% %Comment out if not in use, use for single whiskers 
% classes = dir('*.measurements');
% 
% parfor n=1:length(classes)
%     [~, outputFileName] = fileparts(classes(n).name);
%     system(['classify-single-whisker ' classes(n).name ' ' outputFileName '.measurements']);
%     display([classes(n).name ' has been classified'])
% end

%% (4) RECLASSIFY: Refines previous step
% classes = dir('*.measurements');
% 
% parfor n=1:length(classes)
%     [~, outputFileName] = fileparts(classes(n).name);
%     system(['reclassify ' classes(n).name ' ' outputFileName '.measurements' ' ' '-n ' whiskerNum]);
%     display([classes(n).name ' has been reclassified'])
%     display([classes(n).name ' completed'])
% end
%%
%Please visit http://whiskertracking.janelia.org/wiki/display/MyersLab/Whisker+Tracking+Tutorial
%for more information
%   Clack NG, O'Connor DH, Huber D, Petreanu L, Hires A., Peron, S., Svoboda, K., and Myers, E.W. (2012) 
%   Automated Tracking of Whiskers in Videos of Head Fixed Rodents.
%   PLoS Comput Biol 8(7):e1002591. doi:10.1371/journal.pcbi.1002591