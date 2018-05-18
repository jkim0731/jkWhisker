function norpix_seq_reader_jsy(seqIn)
% Reading seq file and converting to mp4 file directly, without using tiff
% stacks. Seq file info from norpix manual.
% When number of frames are larger than a threshold "nframes_threshold", it
% divides the file into multiple mp4 files. 
% Using ffmpeg to solve weird problem of number of frames being doubled in
% whisker tracker. 

%Open file
tic
seqID = fopen(seqIn);
%Read important header information
fseek(seqID, 548, 'bof');
iWidth = fread(seqID, [1], 'ulong');
iHeight = fread(seqID, [1], 'ulong');
iBitDepth = fread(seqID, [1], 'ulong');
iRBitDepth = fread(seqID, [1], 'ulong'); %Discard, useless
iSize = fread(seqID, [1], 'ulong');
fseek(seqID, 580, 'bof');
iTrueSize = fread(seqID, [1], 'ulong');

seqInfo = dir(seqIn);
seqSize = seqInfo.bytes;
totalFrameNumber = (seqSize-8192)/iTrueSize;
nframes_threshold = 1e5;
ffmpegPath = 'C:\Users\shires\Documents\GitHub\jkWhisker\tracking_windows\ffmpeg';

%Perform sanity checks
if (iWidth*iHeight) ~= iSize
    error('Nonsensical image size')
end

%Create mp4 name
[path,seqName,~] = fileparts(seqIn);
mp4Name = [path filesep seqName '_pre.mp4'];

% Create timestamp file name
tsName = [path, filesep, seqName, '_timestamp.mat'];
tsSec = [];
tsMilli = [];
tsMicro = [];

if totalFrameNumber > nframes_threshold
    frameNum = 0;
    for i = 1 : ceil(totalFrameNumber/nframes_threshold)        
        mp4Name = [path filesep seqName sprintf('_%02d_pre.mp4',i)];
        newVid = VideoWriter(mp4Name, 'MPEG-4');
        open(newVid)
        
        tsName = [path, filesep, seqName, sprintf('_%02dtimestamp.mat',i)];
        tsSec = [];
        tsMilli = [];
        tsMicro = [];

        temp_nframes = 0;
        while temp_nframes < nframes_threshold
            %Read frame from seq
            validPos = fseek(seqID, 8192+iTrueSize*frameNum, 'bof');
            if validPos == -1
                break
            end
            currentFrame = uint8(fread(seqID, [iWidth,iHeight], 'uint8'));
            currentSec = fread(seqID, 1, 'uint32');
            tsSec = [tsSec; currentSec];
            currentMilli = fread(seqID, 1, 'uint16');
            tsMilli = [tsMilli; currentMilli];
            currentMicro = fread(seqID, 1, 'uint16');
            tsMicro = [tsMicro; currentMicro];
            
            [checkW, checkH] = size(currentFrame);
            if currentFrame == -1
                break
            elseif checkW ~= iWidth
                break
            elseif checkH ~= iHeight
                break
            end
            currentFrame = rot90(currentFrame, -1);
            currentFrame = fliplr(currentFrame);
            writeVideo(newVid, currentFrame)
            frameNum = frameNum + 1;
            temp_nframes = temp_nframes + 1;
        end
        close(newVid)
        mp4NameFin = [path filesep seqName sprintf('_%02d.mp4',i)];
        mp4CorrectionCMD = sprintf('%s -i %s -b:v 800k -codec:v mpeg4 -c:a copy %s',ffmpegPath, mp4Name, mp4NameFin);
        system(mp4CorrectionCMD);
        
        save(tsName, 'tsSec', 'tsMicro', 'tsMilli')
    end
else
    %Create and open mp4 file
    newVid = VideoWriter(mp4Name, 'MPEG-4');
    open(newVid)

    %Jump to first frame
    fseek(seqID, 8192, 'bof');
    frameNum = 0;

    %Begin frame extraction and writing loop
    finishedWriting = 0;
    while finishedWriting == 0
        %Read frame from seq
        validPos = fseek(seqID, 8192+iTrueSize*frameNum, 'bof');
        if validPos == -1
            break
        end
        currentFrame = fread(seqID, [iWidth,iHeight], 'uint8');
        currentSec = fread(seqID, 1, 'int32');
        tsSec = [tsSec; currentSec];
        currentMilli = fread(seqID, 1, 'uint16');
        tsMilli = [tsMilli; currentMilli];
        currentMicro = fread(seqID, 1, 'uint16');
        tsMicro = [tsMicro; currentMicro];

        [checkW, checkH] = size(currentFrame);
        if currentFrame == -1
            break
        elseif checkW ~= iWidth
            break
        elseif checkH ~= iHeight
            break
        end
        currentFrame = rot90(currentFrame, -1);
        currentFrame = fliplr(currentFrame);
        %Scale range to be between 0 and 1
        scaleVal = max(max(currentFrame));
        currentFrame = currentFrame/scaleVal;
        writeVideo(newVid, currentFrame)
        frameNum = frameNum + 1;
    end
    close(newVid)
    mp4NameFin = [path filesep seqName '.mp4'];
    mp4CorrectionCMD = sprintf('%s -i %s -b:v 800k -codec:v mpeg4 -c:a copy %s', ffmpegPath, mp4Name, mp4NameFin);    
    system(mp4CorrectionCMD);
    delete(mp4Name);
    save(tsName, 'tsSec', 'tsMicro', 'tsMilli')
end
fclose(seqID);
toc

end

