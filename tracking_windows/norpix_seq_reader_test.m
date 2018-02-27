function [] = norpix_seq_reader_test(seqIn,nFramesThreshold)
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

%Perform sanity checks
if (iWidth*iHeight) ~= iSize
    error('Nonsensical image size')
end

% calculate number of frames. If it exceeds a certain threshold, chop it.
finfo = dir(seqIn);
nTotFrames = (finfo.bytes-8192)/iTrueSize;
if nTotFrames > nFramesThreshold
    nFiles = ceil(nTotFrames/nFramesThreshold);
end

if nFiles == 1
    %Create mp4 name
    [path,seqName,~] = fileparts(seqIn);
    mp4Name = [path filesep seqName '.mp4'];

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
    fclose(seqID);
    toc
else
    %Jump to first frame
    fseek(seqID, 8192, 'bof');
    frameNum = 0;    
    %Create mp4 name
    [path,seqName,~] = fileparts(seqIn);
    for ifile = 1 : nFiles
        %Create mp4 name        
        mp4Name = sprintf(['%s//%0',num2str(numel(num2str(fix(abs(nFiles))))),'d.mp4'],path,ifile);

        %Create and open mp4 file
        newVid = VideoWriter(mp4Name, 'MPEG-4');
        open(newVid)
        
        %Begin frame extraction and writing loop
        currNumFrames = 0;
        while currNumFrames < nFramesThreshold
            %Read frame from seq
            validPos = fseek(seqID, 8192+iTrueSize*frameNum, 'bof');
            if validPos == -1
                break
            end
            currentFrame = fread(seqID, [iWidth,iHeight], 'uint8');
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
            currNumFrames = currNumFrames + 1;
        end
        close(newVid)
    end
    fclose(seqID);
    toc
end
end

