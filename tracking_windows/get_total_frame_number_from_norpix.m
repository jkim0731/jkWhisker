function totalFrameNumber = get_total_frame_number_from_norpix(seqIn)

seqID = fopen(seqIn);
seqInfo = dir(seqIn);
seqSize = seqInfo.bytes;
%Read important header information
fseek(seqID, 548, 'bof');
iWidth = fread(seqID, [1], 'ulong');
iHeight = fread(seqID, [1], 'ulong');
iBitDepth = fread(seqID, [1], 'ulong');
iRBitDepth = fread(seqID, [1], 'ulong'); %Discard, useless
iSize = fread(seqID, [1], 'ulong');
fseek(seqID, 580, 'bof');
iTrueSize = fread(seqID, [1], 'ulong');


totalFrameNumber = (seqSize-8192)/iTrueSize;
fclose(seqID);