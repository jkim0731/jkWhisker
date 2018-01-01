% WHISKER_TRACKER_TRUE_PARALLEL(PIXDEN, WHISKERNUM, FACE, VIDEOTYPE, MPARALLEL)
%performs parallel processing on video files using the Janelia Farm whisker
%tracker.

%% Change Log
% Change                             Editor            Date
%-------------------------------------------------------------------------------
% Initial creation                   J. Sy             2014-08-04
% Initial edits                      S.A. Hires        2014-08-04
% Added error-check                  J. Kim            2017-09-18
% Changed directory recognition      J. Cheung         2017-09-27
% Parallelization w/ MParallel       J. Sy             2017-12-13
% Function renaming and options      J. Sy             2017-12-15

% Notes ------------------------------------------------------------------------
%In theory, sections 3 (CLASSIFY) and 3.5 (CLASSIFY-SINGLE-WHISKER) should
%not be in use at the same time, since they do the same thing
% Created by SAH and JS, 2014-08-04


function [traceTime, totalTime] = whisker_tracker_true_parallel_JK(varargin)
%% (0) SET PARAMETERS: set parameters based on inputs
totalTStart = tic;

%% DEFAULTS
ppm = 17.81002608/2; % for TPM telecentric lens  
pixDen = 1/ppm; %This is the default pixel density
whiskerNum = 2;
vidType = '.mp4';
faceSide = 'bottom';
numCores = feature('numcores'); %Number of cores to run on, by default, all of them
errorCheck = true;

%% Reset inputs
if nargin >= 1
    pixDen = varargin{1};
end
if nargin >= 2
    whiskerNum = varargin{2};
end
if nargin >= 3
    faceSide = varargin{3};
end
if nargin >= 4
    numCores = varargin{4};
end
if nargin >= 5
    errorCheck = varargin{5};
end

%% Find MParallel
if exist('MParallel.exe','file') == 2
    %MParallel is already on the path. Yay!
    mpPath = which('MParallel.exe');
else
    try
        %See if we can get MParallel from the NAS
        system('copy Z:\Software\WhiskerTracking\MParallel.exe startDir');
        mpPath = which('MParallel.exe');
    catch
        error('Cannot find MParallel.exe on the path')
    end
end
%% (1) TRACE: Uses Janelia Farm's whisker tracking software to track all whiskers in a directory
traceTStart = tic;
tracecmd = sprintf(['dir /b *%s | %s --count=%.00f --shell --stdin'...
' --pattern="trace {{0}} {{0:N}}.whiskers"'], vidType, mpPath, numCores);
system(tracecmd);
traceTime = toc(traceTStart);

%% (2) MEASURE: Generates measurements of traced shapes for later evaluation
measurecmd = sprintf(['dir /b *.whiskers | %s --count=%.00f --shell --stdin'...
' --pattern="measure --face %s {{0}} {{0:N}}.measurements"'], mpPath, numCores, faceSide);
system(measurecmd);

%% (2)-1 Error check and re-do the analysis (2017/09/18 JK)
if errorCheck == true
    mp4_flist = dir(['*' vidType]);
    mp4_list = zeros(length(mp4_flist),1);
    for i = 1 : length(mp4_flist)
        mp4_list(i) = str2double(strtok(mp4_flist(i).name,'.')); % assume that all filenames are integers
    end
    whiskers_flist = dir('*.whiskers');
    whiskers_list = zeros(length(whiskers_flist),1);
    for i = 1 : length(whiskers_flist)
        whiskers_list(i) = str2double(strtok(whiskers_flist(i).name,'.'));
    end
    measure_flist = dir('*.measurements');
    measure_list = zeros(length(measure_flist),1);
    for i = 1 : length(measure_flist)
        measure_list(i) = str2double(strtok(measure_flist(i).name,'.'));
    end
    if length(measure_flist) < length(mp4_flist)
        % 1) re-trace those failed to trace before
        if length(whiskers_flist) < length(mp4_flist)
            trace_errorlist = setdiff(mp4_list,whiskers_list);
            for i = 1 : length(trace_errorlist)
                temp_fname = num2str(trace_errorlist(i));
                sout = system(['trace ' temp_fname vidType ' ' temp_fname]);
                trial_ind = 0;
                while (sout ~= 0 && trial_ind < 3)
                    sout = system(['trace ' temp_fname vidType ' ' temp_fname]);
                    trial_ind = trial_ind + 1;
                end
                if sout == 0
                    disp([temp_fname vidType ' traced successfully'])
                else
                    disp(['Failed to trace ' temp_fname])
                end
            end
        end
        % 2) re-measure

        measure_errorlist = setdiff(mp4_list,measure_list);
        for i = 1 : length(measure_errorlist)
            temp_fname = num2str(measure_errorlist(i));
            sout = system(['measure ' '--face ' faceSide ' ' temp_fname '.whiskers ' temp_fname '.measurements']);
            trial_ind = 0;
            while (sout ~= 0 && trial_ind < 3)
                sout = system(['measure ' '--face ' faceSide ' ' temp_fname '.whiskers ' temp_fname '.measurements']);
                trial_ind = trial_ind + 1;
            end
            if sout == 0
                disp([temp_fname '.whiskers has been measured'])
            else
                disp(['Failed to measure ' temp_fname '.whiskers'])
            end
        end
    end
    % 3) for those still remaining not measured, trace again and then  measure them
    measure_flist = dir('*.measurements');
    if length(measure_flist) < length(mp4_flist)
        measure_list = zeros(length(measure_flist),1);
        for i = 1 : length(measure_flist)
            measure_list(i) = str2double(strtok(measure_flist(i).name,'.'));
        end
        measure_errorlist = setdiff(mp4_list,measure_list);
        for i = 1 : length(measure_errorlist)
            temp_fname = num2str(measure_errorlist(i));
            sout = system(['trace ' temp_fname '.mp4 ' temp_fname]);
            trial_ind = 0;
            while (sout ~= 0 && trial_ind < 3)
                sout = system(['trace ' temp_fname '.mp4 ' temp_fname]);
                trial_ind = trial_ind + 1;
            end
            if sout == 0
                disp([temp_fname '.mp4 traced successfully'])
            else
                disp(['Failed to trace ' temp_fname])
            end
        end
        for i = 1 : length(measure_errorlist)
            temp_fname = num2str(measure_errorlist(i));
            sout = system(['measure --face' faceSide ' ' temp_fname '.whiskers ' temp_fname '.measurements']);
            trial_ind = 0;
            while (sout ~= 0 && trial_ind < 3)
                sout = system(['measure --face' faceSide ' ' temp_fname '.whiskers ' temp_fname '.measurements']);
                trial_ind = trial_ind + 1;
            end
            if sout == 0
                disp([temp_fname '.whiskers measured successfully'])
            else
                disp(['Failed to measure ' temp_fname])
            end
        end
    end
end
%% (3) CLASSIFY: Helps refine tracing to more accurately determine which shapes are whiskers
%Use for multiple whiskers

classifycmd = sprintf(['dir /b *.measurements | %s --count=%.00f --shell --stdin'...
' --pattern="classify {{0}} {{0:N}}.measurements %s --px2mm %.02f -n %.00f"'], mpPath, numCores, faceSide, pixDen, whiskerNum);
system(classifycmd);

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
totalTime = toc(totalTStart);
end
%Please visit http://whiskertracking.janelia.org/wiki/display/MyersLab/Whisker+Tracking+Tutorial
%for more information
%   Clack NG, O'Connor DH, Huber D, Petreanu L, Hires A., Peron, S., Svoboda, K., and Myers, E.W. (2012)
%   Automated Tracking of Whiskers in Videos of Head Fixed Rodents.
%   PLoS Comput Biol 8(7):e1002591. doi:10.1371/journal.pcbi.1002591
