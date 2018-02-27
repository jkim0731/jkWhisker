% SEQ_TO_MP4(FLOCATION,ITYPE) converts either a single seq or a directory to mp4
% This is a heavily modified version of Read_Seq_File.m by Sk Sahariyaz Zaman 

function seq_to_mp4_test(fLocation, iType, cores, nFramesThreshold)
  switch iType
  case 'dir'
    %Running on directory
    convertList = dir([fLocation filesep '*.seq']);
    numSEQ = length(convertList);
    %Run in single or parallel depending on cores
    if cores == 1
      for i =1:numSEQ
        norpix_seq_reader_test([fLocation filesep convertList(i).name], nFramesThreshold)
      end
    else
      delete(gcp('nocreate')); %turns off all other parallel pool processes
      pool = parpool(cores);
      parfor i =1:numSEQ
        norpix_seq_reader_test([fLocation filesep convertList(i).name], nFramesThreshold)
      end
    end

  case 'single'
    %Single file, running on file
    norpix_seq_reader_test(fLocation, nFramesThreshold)
  end
end
