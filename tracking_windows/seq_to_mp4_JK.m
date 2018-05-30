% SEQ_TO_MP4(FLOCATION,ITYPE) converts either a single seq or a directory to mp4
% This is a heavily modified version of Read_Seq_File.m by Sk Sahariyaz Zaman 

function seq_to_mp4_JK(fLocation, iType)
    switch iType
        case 'dir'
        %Running on directory
        convertList = dir([fLocation filesep '*.seq']);
        numSEQ = length(convertList);
        
        % Running in single core is better than parfor
        % Haven't tested mparallel.exe yet. 
        for i =1:numSEQ
            norpix_seq_reader_JK([fLocation filesep convertList(i).name])
        end

        case 'single'
            %Single file, running on file
            norpix_seq_reader_JK(fLocation)
    end
end
