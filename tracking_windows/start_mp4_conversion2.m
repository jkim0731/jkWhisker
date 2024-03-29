function [copyTime, convertTime, totalFiles] = start_mp4_conversion2(varargin)
    % START_WHISKER_TRACKING An improved version of the mp4 converter and Janelia farm whisker tracking script
    % START_WHISKER_TRACKING(STARTDIR, ENDDIR) specifies location of current files and where to send them, defaults to working directory
    % START_WHISKER_TRACKING(STARTDIR, ENDDIR, CONVERTVIDEO) CONVERTVIDEO = 1 will convert seq -> mp4, = 0 will skip conversion
    % START_WHISKER_TRACKING(STARTDIR, ENDDIR, CONVERTVIDEO, WHISKERNUM, PIXDEN) Adjust whisker number or pixel density as appropriate
    % START_WHISKER_TRACKING(STARTDIR, ENDDIR, CONVERTVIDEO, WHISKERNUM, PIXDEN, TRANSFERVIDEO) TRANSFERVIDEO = 1 copies files to end directory, =0 will not copy

    %% CHANGE LOG -----------------------------------------------------------

    %  Change                            Date            Editor
    % -----------------------------------------------------------------
    %  Initial Script                   2014-12-13      J. Sy
    %  Variation for Jon                2017-10-06      J. Cheung
    %  Added mp4 converter selector     2017-10-13      J. Sy
    %  Function name change             2017-12-13      J. Sy
    %  Added new parameters             2017-12-15      J. Sy

    %% SECTION 1: INPUT HANDLING -------------------------------------------
    %Default Settings:
    startDir = pwd;
    endDir = pwd;
    convertVideo = 1; %mp4 converter will run
    whiskerNum = 2; %Default number of whiskers is 2
    ppm = 17.81002608/2; % for TPM telecentric lens  
    pixDen = 1/ppm; %This is the default pixel density
    faceSide = 'bottom';
    transferVideo = true; %Default to transfering to end directory
    reservedCores = 0; %Standard is to not reserve any cores for parallel
    errorCheck = true; %Default to performin JK error check

    if nargin == 0
        fprintf('Function called empty, setting all parameters to default')
        fprintf('Everything will be processed in the current directory')
    else
        %Handle variable argument processing in subfunction so all the code people don't care about can be at the bottom
        [startDir, endDir, convertVideo, whiskerNum, pixDen, faceSide, reservedCores, errorCheck] = varInputHandler(varargin);
    end
    %If startDir and endDir are the same, don't transfer Video
    if strcmp(startDir, endDir)
        transferVideo = false;
    end
    numCores = feature('numcores'); %identify number of cores available for MATLAB to use
    numCores = numCores - reservedCores; 
  
    %% SECTION 2: CONVERT .SEQ TO .MP4 --------------------------------------
    mp4List = dir([startDir filesep '*mp4']);
    seqList = dir([startDir filesep '*seq']);
    if convertVideo == 1
        if isempty(seqList) 
            return;
        end
        
        if ~isempty(mp4List)
            delete([startDir,'\*.mp4']);
        end
        
        tic
        fprintf('STARTING MP4 CONVERSION OF %s \n', startDir)
        seq_to_mp4_JK1(startDir,'dir', 1)
        fprintf('FINISHED MP4 CONVERSION \n')
        convertTime = toc;
    else
        fprintf('Skipping video conversion in %s \n', startDir)
    end
    totalFiles = length(dir('*.mp4'));
    %% SECTION 4: COPY FILES ------------------------------------------------

    if transferVideo
        tic;
        system(['copy ', startDir, '\*.mp4 ', endDir]);
        copyTime = toc;
        delete([startDir, '\*.mp4']);
    end

end


    %% ------------------------------------------------------------------------------
    %And here is a subfunction to handle the varargin stuff because most people don't want to see this
function [startDir, endDir, convertVideo, whiskerNum, pixDen, faceSide, rCores, eCheck] = varInputHandler(vInputs)
    %Defaults
    startDir = pwd;
    endDir = pwd;
    convertVideo = 1; %mp4 converter will run
    whiskerNum = 2; %Default number of whiskers is 1
    ppm = 17.81002608/2; % for TPM telecentric lens  
    pixDen = 1/ppm; %This is the default pixel density
    rCores = 0; %Standard is to not reserve any cores for parallel
    eCheck = true; %Default to perform JK error check

    %PARAMETER 1: Where to find files
    if length(vInputs) >= 1 %Change start directory
    %Empty call
        if ~isempty(vInputs{1})
            %Check if actually a directory
            if ischar(vInputs{1})
                startDir = vInputs{1};
            else
                error('Directory name must be string')
            end
            %Check if exist
            if exist(vInputs{1},'dir') ~= 7
                error('Cannot find input directory %c', vInputs{1})
            end
        end
    end
    %PARAMETER 2: Where to transfer files
    if length(vInputs) >= 2 %Change end directory
    %Empty call
        if ~isempty(vInputs{2})
            %Check if actually a directory
            if ischar(vInputs{2})
                endDir = vInputs{2};
            else
                error('Directory name must be string')
            end
                %Check if exist
            if exist(vInputs{2},'dir') ~= 7
                error('Cannot find output directory %c', vInputs{2})
            end
        end
    end
    %PARAMETER 3: Do we need to convert SEQs?
    if length(vInputs) >= 3
        if vInputs{3} == 1 || strcmpi(vInputs{3},'true')
            convertVideo = 1;
        elseif vInputs{3} == 0 || strcmpi(vInputs{3},'false')
            convertVideo = 0;
        else
            fprintf('Improper entry for video parameter, defaulting to convertVideo = TRUE')
        end
    end
    %PARAMETER 4: Number of whiskers to track
    if length(vInputs) >= 4
        if ~isempty(vInputs{4})
            whiskerNum = vInputs{4};
        end
    end
    %PARAMETER 5: Pixel density
    if length(vInputs) >= 5
        if ~isempty(vInputs{5}) && isnumeric(vInputs{5})
            pixDen = vInputs{5};
        end
    end
    %PARAMETER 6: Where is the face in the image?
    if length(vInputs) >= 6
        if ischar(vInputs{6})
            faceSide = vInputs{6};
        else
            warning('Invalid face side input, defaulting to "top"')
        end
    end
    %PARAMETER 7: Should we save any CPU cores for non-tracking?
    if length(vInputs) >= 7
        if isnumeric(vInputs{7})
            rCores = vInputs{7};
        end
    end
    %PARAMETER 8: Whether or not to perform Jinho's error check (the check will not be parallel)
    if length(vInputs) >= 8
        switch vInputs{8}
            case 1
                eCheck = true;
            case 0
                eCheck = false;
            otherwise
                warning('Invalid error check input')
                return
        end
    end
    %EXTRA PARAMETERS
    if length(vInputs) > 8
        fprintf('You called this function with too many variables, ignoring extras')
    end
end
    