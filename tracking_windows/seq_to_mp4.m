% SEQ_TO_MP4(FLOCATION,ITYPE) converts either a single seq or a directory to mp4
% This is a heavily modified version of Read_Seq_File.m by Sk Sahariyaz Zaman 

function seq_to_mp4(fLocation, iType, cores)
  switch iType
  case 'dir'
    %Running on directory
    convertList = dir([fLocation filesep '*.seq']);
    numSEQ = length(convertList);
    %Run in single or parallel depending on cores
    if cores == 1
      for i =1:numSEQ
        read_the_seq([fLocation filesep convertList(i).name])
      end
    else
      delete(gcp('nocreate')); %turns off all other parallel pool processes
      pool = parpool(cores);
%       parfor i =1:numSEQ
      for i =1:numSEQ
        read_the_seq([fLocation filesep convertList(i).name])
      end
    end

  case 'single'
    %Single file, running on file
    read_the_seq(fLocation)
  end
end

function read_the_seq(seqFullPath)
  %Extract out the seq's name and path 
  [seqPath,seqFile,~] = fileparts(seqFullPath);
  outputMP4 = [seqPath filesep seqFile '.mp4'];

  %Call MATLAB's videoWriter function, tell it to write to mp4
  mp4Writer = VideoWriter(outputMP4, 'MPEG-4');
  mp4Writer.Quality = 100; %Quality could perhaps be higher? More experimentation needed

  open(mp4Writer);
  frameWindow = 1; %Amount of frames to move by, should be 1 unless you want 
  %to intentionally skip frames

  % Changed from the original Zaman code to just read from 0:Inf, since
  % that's all our lab actually needs
  firstFrame=0;
  lastFrame=inf;

  %Calls seqIo to read the SEQ file. seqIo, and its dependencies,
  %seqReaderPlugin and seqWriterPlugin, are sourced from Piotr's toolbox,
  %found at https://github.com/pdollar/toolbox/tree/master/videos
  seqID = seqIo(seqFullPath, 'reader' );
  info=seqID.getinfo();
  lastFrame=min(lastFrame,info.numFrames-1);
  allFrames=firstFrame:frameWindow:lastFrame;


  for currentFrame=allFrames
%         seqID.seek(currentFrame);
%         [frame2write, ~] =seqID.getframe();
%         frame2write = circshift(frame2write, [-16, 33]);        
%         %
%         moveTheseBits = circshift(frame2write(end-15:end,:), 321, 2);
%         frame2write(end-15:end,:) = moveTheseBits;
%         moveTheseBits = circshift(frame2write(:,1:33), 1, 1);
%         frame2write(:,1:33) = moveTheseBits;
%         writeVideo(mp4Writer, frame2write);
%         clear frame2write moveTheseBits;
%         
        
        seqID.seek(currentFrame);
        [initialFrame, ~] = seqID.getframe();
        adjustedFrame = initialFrame;
        %
        moveTheseBits = initialFrame(1:16,:);
        moveTheseBits = circshift(moveTheseBits, 354, 2);
        adjustedFrame(1:16,:) = [];
        adjustedFrame = circshift(adjustedFrame, 33,2);
        adjustedFrame = [adjustedFrame; moveTheseBits];
        %
        
        writeVideo(mp4Writer, adjustedFrame);
        clear initialFrame adjustedFrame;
        
  end

  seqID.close();
  close(mp4Writer);

end
