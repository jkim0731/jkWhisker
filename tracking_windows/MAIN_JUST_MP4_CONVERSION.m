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
NUMBER_OF_RESERVED_CORES = 0;
USE_ERROR_CHECK = true;
base_dir = 'E:\';

% DIRECTORY LIST ----------------------------------------------------------
%Directory to find files to convert or track
cd(base_dir)
dirlist = dir('JK052piezo*');
convertVid = {true, false, true, false};
% =========================================================================
% Queue settings should only be changed above this line
% =========================================================================
% %%
% %TRACKING STARTS HERE!
% allFileClock = tic;
% sumFiles1 = 0;
% sumTrackTime = 0;
% sumConversionTime1 = 0;
% sumCopyTime1 = 0;
% for i = 1:length(dirlist)
%     if dirlist(i).isdir
%         startDir = [base_dir, dirlist(i).name];
%         endDir = [base_dir, 'tracked\', dirlist(i).name];
%         system(['mkdir ', endDir])
%         if CONVERT_ALL_VIDEOS == true
%             [copyTime, convertTime, totalFiles] = start_mp4_conversion1(...
%                 startDir, endDir, 1, NUMBER_OF_WHISKERS, PIXEL_DENSITY, ...
%                 FACE_LOCATION, NUMBER_OF_RESERVED_CORES, USE_ERROR_CHECK);
%         else
%             [copyTime, convertTime, totalFiles] = start_mp4_conversion1(...
%                 startDir, endDir, convertVid{i}, NUMBER_OF_WHISKERS, PIXEL_DENSITY, ...
%                 FACE_LOCATION, NUMBER_OF_RESERVED_CORES, USE_ERROR_CHECK);
%         end
%         sumFiles1 = sumFiles1 + totalFiles;
%         sumConversionTime1 = sumConversionTime1 + convertTime;
%         sumCopyTime1 = sumCopyTime1 + copyTime;
%     end
% end
% 
% %STAT SECTION -------------------------------------------------------------
% if DISPLAY_STATISTICS == true
%   allFileTime = toc(allFileClock);
%   totalHours = floor(allFileTime/3600);
%   extraMinutes = floor(rem(allFileTime,3600)/60);
%   extraSeconds = rem(rem(allFileTime,3600),60);
%   %Some more math
%   convPct = 100*(sumConversionTime/allFileTime);  
%   timePerFile = allFileTime/sumFiles;
%   %Display stats
%   fprintf('Tracking statistics: \n')
%   fprintf('Total time: %.00f hours %.00f minutes %.02f seconds \n', ...
%   totalHours, extraMinutes, extraSeconds)
%   fprintf('Of this time: \n')
%   fprintf('Video conversion took %.02f seconds or %.02f percent of the total time \n', sumConversionTime, convPct)
%   fprintf('You tracked a total of %.00f files with an average time of %.02f seconds per file \n', sumFiles, timePerFile)
%   fprintf('Aaaaaaaaaaaaand we''re done! \n')
% 
% end
% 
allFileClock = tic;
sumFiles2 = 0;
sumTrackTime = 0;
sumConversionTime2 = 0;
sumCopyTime2 = 0;
for i = 1:length(dirlist)
    if dirlist(i).isdir
        startDir = [base_dir, dirlist(i).name];
        endDir = [base_dir, 'tracked2\', dirlist(i).name];
        system(['mkdir ', endDir])
        if CONVERT_ALL_VIDEOS == true
            [copyTime, convertTime, totalFiles] = start_mp4_conversion2(...
                startDir, endDir, 1, NUMBER_OF_WHISKERS, PIXEL_DENSITY, ...
                FACE_LOCATION, NUMBER_OF_RESERVED_CORES, USE_ERROR_CHECK);
        else
            [copyTime, convertTime, totalFiles] = start_mp4_conversion2(...
                startDir, endDir, convertVid{i}, NUMBER_OF_WHISKERS, PIXEL_DENSITY, ...
                FACE_LOCATION, NUMBER_OF_RESERVED_CORES, USE_ERROR_CHECK);
        end
        sumFiles2 = sumFiles2 + totalFiles;
        sumConversionTime2 = sumConversionTime2 + convertTime;
        sumCopyTime2 = sumCopyTime2 + copyTime;
    end
end
toc
% 
% %STAT SECTION -------------------------------------------------------------
if DISPLAY_STATISTICS == true
  allFileTime = toc(allFileClock);
  totalHours = floor(allFileTime/3600);
  extraMinutes = floor(rem(allFileTime,3600)/60);
  extraSeconds = rem(rem(allFileTime,3600),60);
  %Some more math
  convPct = 100*(sumConversionTime/allFileTime);  
  timePerFile = allFileTime/sumFiles;
  %Display stats
  fprintf('Tracking statistics: \n')
  fprintf('Total time: %.00f hours %.00f minutes %.02f seconds \n', ...
  totalHours, extraMinutes, extraSeconds)
  fprintf('Of this time: \n')
  fprintf('Video conversion took %.02f seconds or %.02f percent of the total time \n', sumConversionTime, convPct)
  fprintf('You tracked a total of %.00f files with an average time of %.02f seconds per file \n', sumFiles, timePerFile)
  fprintf('Aaaaaaaaaaaaand we''re done! \n')

end

% allFileClock = tic;
% sumFiles3 = 0;
% sumTrackTime = 0;
% sumConversionTime3 = 0;
% sumCopyTime3 = 0;
% for i = 1:length(dirlist)
%     if dirlist(i).isdir
%         startDir = [base_dir, dirlist(i).name];
%         endDir = [base_dir, 'tracked3\', dirlist(i).name];
%         system(['mkdir ', endDir])
%         if CONVERT_ALL_VIDEOS == true
%             [copyTime, convertTime, totalFiles] = start_mp4_conversion3(...
%                 startDir, endDir, 1, NUMBER_OF_WHISKERS, PIXEL_DENSITY, ...
%                 FACE_LOCATION, NUMBER_OF_RESERVED_CORES, USE_ERROR_CHECK);
%         else
%             [copyTime, convertTime, totalFiles] = start_mp4_conversion3(...
%                 startDir, endDir, convertVid{i}, NUMBER_OF_WHISKERS, PIXEL_DENSITY, ...
%                 FACE_LOCATION, NUMBER_OF_RESERVED_CORES, USE_ERROR_CHECK);
%         end
%         sumFiles3 = sumFiles3 + totalFiles;
%         sumConversionTime3 = sumConversionTime3 + convertTime;
%         sumCopyTime3 = sumCopyTime3 + copyTime;
%     end
% end
% 
% %STAT SECTION -------------------------------------------------------------
% if DISPLAY_STATISTICS == true
%   allFileTime = toc(allFileClock);
%   totalHours = floor(allFileTime/3600);
%   extraMinutes = floor(rem(allFileTime,3600)/60);
%   extraSeconds = rem(rem(allFileTime,3600),60);
%   %Some more math
%   convPct = 100*(sumConversionTime/allFileTime);  
%   timePerFile = allFileTime/sumFiles;
%   %Display stats
%   fprintf('Tracking statistics: \n')
%   fprintf('Total time: %.00f hours %.00f minutes %.02f seconds \n', ...
%   totalHours, extraMinutes, extraSeconds)
%   fprintf('Of this time: \n')
%   fprintf('Video conversion took %.02f seconds or %.02f percent of the total time \n', sumConversionTime, convPct)
%   fprintf('You tracked a total of %.00f files with an average time of %.02f seconds per file \n', sumFiles, timePerFile)
%   fprintf('Aaaaaaaaaaaaand we''re done! \n')
% 
% end
% 
% allFileClock = tic;
% sumFiles4 = 0;
% sumTrackTime = 0;
% sumConversionTime4 = 0;
% sumCopyTime4 = 0;
% for i = 1:length(dirlist)
%     if dirlist(i).isdir
%         startDir = [base_dir, dirlist(i).name];
%         endDir = [base_dir, 'tracked4\', dirlist(i).name];
%         system(['mkdir ', endDir])
%         if CONVERT_ALL_VIDEOS == true
%             [copyTime, convertTime, totalFiles] = start_mp4_conversion4(...
%                 startDir, endDir, 1, NUMBER_OF_WHISKERS, PIXEL_DENSITY, ...
%                 FACE_LOCATION, NUMBER_OF_RESERVED_CORES, USE_ERROR_CHECK);
%         else
%             [copyTime, convertTime, totalFiles] = start_mp4_conversion4(...
%                 startDir, endDir, convertVid{i}, NUMBER_OF_WHISKERS, PIXEL_DENSITY, ...
%                 FACE_LOCATION, NUMBER_OF_RESERVED_CORES, USE_ERROR_CHECK);
%         end
%         sumFiles4 = sumFiles4 + totalFiles;
%         sumConversionTime4 = sumConversionTime4 + convertTime;
%         sumCopyTime4 = sumCopyTime4 + copyTime;
%     end
% end
% 
% %STAT SECTION -------------------------------------------------------------
% if DISPLAY_STATISTICS == true
%   allFileTime = toc(allFileClock);
%   totalHours = floor(allFileTime/3600);
%   extraMinutes = floor(rem(allFileTime,3600)/60);
%   extraSeconds = rem(rem(allFileTime,3600),60);
%   %Some more math
%   convPct = 100*(sumConversionTime/allFileTime);  
%   timePerFile = allFileTime/sumFiles;
%   %Display stats
%   fprintf('Tracking statistics: \n')
%   fprintf('Total time: %.00f hours %.00f minutes %.02f seconds \n', ...
%   totalHours, extraMinutes, extraSeconds)
%   fprintf('Of this time: \n')
%   fprintf('Video conversion took %.02f seconds or %.02f percent of the total time \n', sumConversionTime, convPct)
%   fprintf('You tracked a total of %.00f files with an average time of %.02f seconds per file \n', sumFiles, timePerFile)
%   fprintf('Aaaaaaaaaaaaand we''re done! \n')
% 
% end
