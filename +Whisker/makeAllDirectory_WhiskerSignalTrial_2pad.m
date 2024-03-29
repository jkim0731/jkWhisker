function makeAllDirectory_WhiskerSignalTrial_2pad(d,varargin)
%
%
%
%   USAGE:
%
%
%   INPUTS:
%
%   d: Directory path name as string.
%
%
%   Optional parameter/value pair arguments:
%
%           'include_files': Optional cell array of strings giving file name prefixes
%                           of files in directory 'd' to process. Files will be processed
%                           in the order they are given within this cell array. *NOTE*: If
%                           this argument is not given, *all* '_WT.mat' files in directory
%                           'd' will be processed.
%
%           'ignore_files': Optional cell array of strings giving file name prefixes
%                           (i.e., file names without the '_WT.mat' suffix/extension) to ignore.
%                           Trumps 'include_files' argument.
%
%           'polyRoiInPix':
%               Sets arc-length limits (in pixels) on which to perform secondary curve fitting.
%               This argument can be given in two forms:
%               (1) an 1x2 vector that gives the ROI for *all* whiskers; or
%               (2) a cell array where first element is a vector of trajectory IDs
%                   (of length N) and subsequent elements comprise N 1x2 vectors
%                   giving ROIs for the trajectory IDs specified in the first
%                   element (respectively).
%               Limits are inclusive.
%
%           'barPosOffset': [xOffset yOffset].
%
%           'follicleExtrapDistInPix': Distance to extrapolate past the end of the tracked whisker
%                              in order to estimate follicle coordinates. If this argument is
%                              given, follicle position will be estimated. Default is not to estimate.
%                              Presently can only give one value, which will be used for all TIDs. Need to
%                              improve this.
%
%           'polyFitsMask': [xpoints; ypoints]. See Whisker.WhiskerSignalTrial polyFitsMask property.
%
%
%   DESCRIPTION:
%
%   Requires WhiskerTrial objects to be saved, as .mat files, in the
%   directory specified by argument 'd'.  These files are read in one at a time and
%   converted to WhiskerSignalTrial objects, which are then saved to disk in the same directory
%   as .mat files with a '_WST.mat' suffix/extension.
%
%   Processes all trajectory IDs within each WhiskerTrial.
%
%
% 3/10, DHO.
%

% Sometimes there is an error with VideoReader during parallel processing (~0.2%). 
% In that or any other case when WhiskerTrial does not work correctly, throw an error
% (fname_errorWST.mat file with trial_num in it)
% 2017/05/15 JK

p = inputParser;

p.addRequired('d', @ischar);
p.addParameter('include_files', {}, @(x) all(cellfun(@ischar,x)));
p.addParameter('ignore_files', {}, @(x) all(cellfun(@ischar,x)));
p.addParameter('polyRoiInPix', [0 200], @(x) numel(x)==2 && x(1)>=0);
p.addParameter('barPosOffset', [0 0], @(x) numel(x)==2);
p.addParameter('follicleExtrapDistInPix', NaN, @isnumeric);
p.addParameter('polyFitsMask', [], @isnumeric);

p.parse(d,varargin{:});

disp 'List of all arguments:'
disp(p.Results)


if ~strcmp(d(end), filesep)
    d = [d filesep];
end

currentDir = pwd;
cd(d)

fnall = arrayfun(@(x) x.name(1:(end-7)), dir([d '*_WT.mat']),'UniformOutput',false);

if ~isempty(p.Results.include_files) % Make sure files are found. If not, ignored.
    ind = ismember(p.Results.include_files,fnall);
    fnall = p.Results.include_files(ind);
    if sum(ind) ~= numel(ind)
        disp('The following files in ''include_files'' were not found in directory ''d'' and will be skipped:')
        disp(p.Results.include_files(~ind))
    end
end

if ~isempty(p.Results.ignore_files)
    ind = ~ismember(fnall,p.Results.ignore_files);
    fnall = fnall(ind);
end

inBoth = intersect(p.Results.include_files,p.Results.ignore_files);
if ~isempty(inBoth)
    disp('The following files were given in BOTH ''include_files'' and ''ignore files'' and will be ignored:')
    disp(inBoth)
end


nfiles = length(fnall);


% Calculating the mean slope and origin for pole axis top
% From the same radial distance
% Except for 90 degrees pole angle (calculate its slope from grand average
% of all other angles)

wta = {};
for fi = 1 : nfiles
    load([fnall{fi}, '_WT.mat']);
    if ~strcmp(w.trialType, 'oo') && ~isempty(w.trialType) % trial type is empty when it's piezo stimulation.
        wta{end+1} = w;
    end
end
if size(wta,1) < size(wta,2)
    wta = wta';
end

if ~isempty(wta)
    angles = unique(cellfun(@(x) x.angle, wta));
    rds = unique(cellfun(@(x) x.radialDistance, wta)); % radial distances
    rds = rds(find(rds));
    meanAxisTop = cell(length(angles), length(rds));
    for ai = 1 : length(angles)
        for ri = 1 : length(rds)        
            axesUpX = cellfun(@(x) (x.angle == angles(ai) && x.radialDistance == rds(ri)) * x.poleAxesUp{1}(1,:), wta, 'uniformoutput', false);
            axesUpY = cellfun(@(x) (x.angle == angles(ai) && x.radialDistance == rds(ri)) * x.poleAxesUp{1}(2,:), wta, 'uniformoutput', false);
            axesUpX = cell2mat(axesUpX);
            axesUpY = cell2mat(axesUpY);
            inds = find(mean(axesUpX,2));
            axesUpX = axesUpX(inds,:);
            axesUpY = axesUpY(inds,:);
            meanAxisTop{ai, ri} = [mean(axesUpX);mean(axesUpY)];
        end    
    end
    if ismember(90, angles) && length(unique(angles)) > 1 % when there is 90 degrees pole and also there are other angles in the same session (for example, pre sessions have only 90 degrees)    
        ind = find(angles == 90);
        rds90 = unique(cellfun(@(x) (x.angle == 90) * x.radialDistance, wta));
        rds90 = rds90(find(rds90));
        for ri = 1 : length(rds90)
            rdi = find(rds == rds90(ri));
            axesUpX = [];
            axesUpY = [];
            for ai = 1 : length(angles)
                if angles(ai) ~= 90
                    axesUpX = [axesUpX; meanAxisTop{ai,rdi}(1,:)];
                    axesUpY = [axesUpY; meanAxisTop{ai,rdi}(2,:)];
                end
            end
            if ~isempty(axesUpX) && ~isempty(axesUpY)
                meanAxisTop{ind,ri} = [mean(axesUpX); mean(axesUpY)];
            end
        end
    end
else
    angles = []; rds = []; meanAxisTop = {};
end







% fnall = {'105','296','297','302','303','305','306','308','309','310','312','313','314','316','317','318','322','323','328','329','331','334','335','336','339','341','345','348','350','351','352','356','358','359','360','363','367','382','385','389','390','391','507','508','510','527','641'};







nfiles = length(fnall);

if ~isempty(fnall)
    if exist('parfor','builtin') % Parallel Computing Toolbox is installed.
        parfor k=1:nfiles
%         for k=1:nfiles
            fn = fnall{k};
            disp(['Processing ''_WT.mat'' file '  fn ', ' int2str(k) ' of ' int2str(nfiles)])
            
            w = pctload([fn '_WT.mat']);

%                 puf_ind = find(puf_fn_list == str2double(fn)); % ?? What
%                 is this?? 2018/03/17 JK

            ai = find(angles == w.angle);
            ri = find(rds == w.radialDistance);
            if ~isempty(ai) && ~isempty(ri)
                ws = Whisker.WhiskerSignalTrial_2pad(w,meanAxisTop{ai,ri},'polyRoiInPix',p.Results.polyRoiInPix);
            else
                ws = Whisker.WhiskerSignalTrial_2pad(w,[],'polyRoiInPix',p.Results.polyRoiInPix);
            end
            if ~isempty(p.Results.polyFitsMask)
                x = p.Results.polyFitsMask(1,:);
                y = p.Results.polyFitsMask(2,:);
                tidList = ws.trajectoryIDs;
                for tid=tidList
                    ws.set_mask_from_points(tid,x,y);
                end
            end

            if ~isnan(p.Results.follicleExtrapDistInPix)
                ws.recompute_cached_follicle_coords(p.Results.follicleExtrapDistInPix,ws.trajectoryIDs); % Right now fits even "contact detection" tids, need to change format***
            end

            outfn = [fn '_WST.mat'];

            pctsave(outfn,ws)
        end
    else
        for k=1:nfiles
            fn = fnall{k};
            disp(['Processing ''_WT.mat'' file '  fn ', ' int2str(k) ' of ' int2str(nfiles)])
            
            load([fn '_WT.mat'],'w');

            ai = find(angles == w.angle);
            ri = find(rds == w.radialDistance);
            if ~isempty(ai) && ~isempty(ri)
                ws = Whisker.WhiskerSignalTrial_2pad(w,meanAxisTop{ai,ri},'polyRoiInPix',p.Results.polyRoiInPix);
            else
                ws = Whisker.WhiskerSignalTrial_2pad(w,[],'polyRoiInPix',p.Results.polyRoiInPix);
            end

            if ~isempty(p.Results.polyFitsMask)
                x = p.Results.polyFitsMask(1,:);
                y = p.Results.polyFitsMask(2,:);
                tidList = ws.trajectoryIDs;
                for tid=tidList
                    ws.set_mask_from_points(tid,x,y);
                end
            end

            if ~isnan(p.Results.follicleExtrapDistInPix)
                ws.recompute_cached_follicle_coords(p.Results.follicleExtrapDistInPix,ws.trajectoryIDs); % Right now fits even "contact detection" tids, need to change format***
            end

            outfn = [fn '_WST.mat'];

            save(outfn,'ws');
        end
    end
end

cd(currentDir)
end

function pctsave(outfn,ws)
save(outfn,'ws');
end

function w = pctload(fn)
load(fn,'w');
end






