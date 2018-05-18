% TRACKING_QUE_WITH_STATS- queues up different directories for tracking assignments
% The spiritual sucessor to whiskerQueue

% =========================================================================
% ADJUSTABLE SETTINGS
% =========================================================================
% TRACKING PARAMETERS -----------------------------------------------------
NUMBER_OF_WHISKERS = 2;
ppm = 17.81002608; % for TPM telecentric lens w/o binning
PIXEL_DENSITY = 1/ppm;
FACE_LOCATION = 'bottom';
DISPLAY_STATISTICS = true;
CONVERT_ALL_VIDEOS = true;
NUMBER_OF_RESERVED_CORES = 1;
USE_ERROR_CHECK = true;
base_dir = 'D:\JK\';

% DIRECTORY LIST ----------------------------------------------------------
%Directory to find files to convert or track
cd(base_dir)
dirlist = dir('JK*');
convertVid = {true, false, true, false};
% =========================================================================
% Queue settings should only be changed above this line
% =========================================================================
%%
%TRACKING STARTS HERE!
allFileClock = tic;
sumFiles = 0;
sumTrackTime = 0;
sumConversionTime = 0;
sumCopyTime = 0;
for i = 1:length(dirlist)
    if dirlist(i).isdir
        startDir = [base_dir, dirlist(i).name];
        endDir = [base_dir, 'tracked\', dirlist(i).name];
        system(['mkdir ', endDir])
        if CONVERT_ALL_VIDEOS == true
            [copyTime, trackTime, convertTime, totalFiles] = start_whisker_tracking(...
                startDir, endDir, 1, NUMBER_OF_WHISKERS, PIXEL_DENSITY, ...
                FACE_LOCATION, NUMBER_OF_RESERVED_CORES, USE_ERROR_CHECK);
        else
            [copyTime, trackTime, convertTime, totalFiles] = start_whisker_tracking(...
                startDir, endDir, convertVid{i}, NUMBER_OF_WHISKERS, PIXEL_DENSITY, ...
                FACE_LOCATION, NUMBER_OF_RESERVED_CORES, USE_ERROR_CHECK);
        end
        sumFiles = sumFiles + totalFiles;
        sumTrackTime = sumTrackTime + trackTime;
        sumConversionTime = sumConversionTime + convertTime;
        sumCopyTime = sumCopyTime + copyTime;
    end
end

%STAT SECTION -------------------------------------------------------------
if DISPLAY_STATISTICS == true
  allFileTime = toc(allFileClock);
  totalHours = floor(allFileTime/3600);
  extraMinutes = floor(rem(allFileTime,3600)/60);
  extraSeconds = rem(rem(allFileTime,3600),60);
  %Some more math
  convPct = 100*(sumConversionTime/allFileTime);
  trackPct = 100*(sumTrackTime/allFileTime);
  timePerFile = allFileTime/sumFiles;
  %Display stats
  fprintf('Tracking statistics: \n')
  fprintf('Total time: %.00f hours %.00f minutes %.02f seconds \n', ...
  totalHours, extraMinutes, extraSeconds)
  fprintf('Of this time: \n')
  fprintf('Video conversion took %.02f seconds or %.02f percent of the total time \n', sumConversionTime, convPct)
  fprintf('Whisker tracking took %.02f seconds or %.02f percent of the total time \n', sumTrackTime, trackPct)
  fprintf('You tracked a total of %.00f files with an average time of %.02f seconds per file \n', sumFiles, timePerFile)
  fprintf('Aaaaaaaaaaaaand we''re done! \n')

end
