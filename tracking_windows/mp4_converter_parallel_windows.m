% MP4_CONVERTER_PARALLEL converts .seqs to .mp4s within specified directory. If no directory is specified the converter will run in the current directory.
%{ 
Edit History
Editor         Date             Change
---------------------------------------------------------------------------
S. Peron       2014-03-01       Original python script recieved 
J. Sy          2017-08-31       MATLAB Function created from python script
J. Sy          2017-09-12       Initial functionality of script 
---------------------------------------------------------------------------
>>> All Hires Lab members please email Jonathan Sy for questions. <<<

 
Dependencies         Location
---------------------------------------------------------------------------
raw2tiff             http://gnuwin32.sourceforge.net/packages/tiff.htm
ffmpeg               https://www.ffmpeg.org/download.html
---------------------------------------------------------------------------

Called with format mp4_converter_parallel(directory, temp directory location, header, parallel [Y/N]) 
%} 
function [] = mp4_converter_parallel_windows(varargin)
parallelRunYes = 1; %Converter runs in parallel by default
headerSize = 8192; %Our current default header, changed with a  varargin call
tempPath = 'C:\Users\shires\Documents\scratch'; %Default scratch folder
raw2tiffPath = '\\dtsfs5\FSC_V2\BISC\usr\shires\Documents\Code\raw2tiff';
ffmpegPath = '\\dtsfs5\FSC_V2\BISC\usr\shires\Documents\Code\ffmpeg';

%First check where to look for .seqs
if nargin == 0
    fprintf('No directory specified, running with current directory \n')
    seqDir = pwd; %Set to whatever folder user is currently in
else
    if varargin{1} == 0
        fprintf('Running with current directory \n')
        seqDir = pwd;
    elseif ~ischar(varargin{1})
        error('Directory name must be a string! \n')
    else
        seqDir = varargin{1};
    end
end
%Then check where to make the tiff files (should have some space)
if nargin >= 2
    if varargin{2} == 0
        fprintf('Running with default scratch directory \n')
    elseif ~ischar(varargin{1})
        error('Directory name must be a string! \n')
    else
        tempPath = varargin{2};
    end
end 
%Check header size
if nargin >= 3
    if ischar(varargin{3})
        fprintf('Using default header size of %d \n', headerSize)
    else
        headerSize = varargin{3};
    end 
end 
%Check if running in parallel
if nargin >= 4
    parallelRunYes = varargin{4};
end 
%Tell people to stop being idiots if they call with more arguments
if nargin >= 5
    fprintf('The converter currently supports only four arguments, all additional arguments will be ignored')
end 


cd(seqDir)
seqList = dir('*.seq');
numSeq = numel(seqList);
if numSeq == 0
    error('No .seq files in directory')
end

if parallelRunYes == 1
    coreNum = feature('numcores');
    delete(gcp('nocreate'));
    parpool(coreNum);
    parfor i = 1:48 %CHANGE BACK CHANGE BACK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        process_to_mp4(seqList(i).name, seqDir, headerSize, tempPath, ffmpegPath, raw2tiffPath)
    end
else
    for i = 1:numSeq
        process_to_mp4(seqList(i).name, seqDir, headerSize, tempPath, ffmpegPath, raw2tiffPath)
    end
end

end

function process_to_mp4(myFile, myDir, header, temp, fPath, rPath)
    %This function does the actual processing of things
    myPath = [myDir filesep myFile];
    [~,fname] = fileparts(myFile); 
    tempDir = [temp filesep fname '_tmp'];  
    %Note: original code had hardcoded setting here on whether or not to
    %clear the seqs or tiffs. Add soft-code later if people care
    
    %system(['chmod -R 777' myDir filesep]) %Need Windows equivalent
    
    fID = fopen(myPath, 'r'); %Read info from relevant file 
    fseek(fID, 548, 0);
    %Extract video information, unpack as signed long integer,
    %little-endian
    width = fread(fID, 1, 'long', 'l');
    height = fread(fID, 1, 'long', 'l');
    depth = fread(fID, 1, 'long', 'l'); %Copied from python but we do jack all with depth
    fseek(fID, 12, 0);
    nFrames = fread(fID, 1, 'long', 'l');
    fseek(fID, 4, 0);
    trueSize = fread(fID, 1, 'long', 'l');
    fclose(fID); 
    
    fprintf('%f \n',width)
    fprintf('%f \n',height)
    
    system(['mkdir ' tempDir])
    fprintf('Created temporary directory: %s',tempDir)
    
    %Run raw2tiff to createa tiff stack out of the .seqs
    for i = 0:nFrames
        offset = header + (i*trueSize);
        raw2tiffCMD = sprintf('%s -H %d -c none -M -w %d -l %d %s %s\\%05d.tif', rPath, offset, width, height, myPath, tempDir, i+1);
        system(raw2tiffCMD);
    end 
    
    %Now run ffmpeg to create a .mp4 out of the tiff stack
    ffmpegOutPath = [myDir filesep fname '.mp4'];
    ffmpegCMD = sprintf('%s -y -i %s\\%%5d.tif -b:v 800k -codec:v mpeg4 %s', fPath, tempDir, ffmpegOutPath);
    system(ffmpegCMD);
    fprintf('Converted file %s successfully', fname) 
    
    %Clear the tiff stack this processor just created to save space
    system(['rmdir ' tempDir '/s /q'])    
end 

