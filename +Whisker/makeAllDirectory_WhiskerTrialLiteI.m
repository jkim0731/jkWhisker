function makeAllDirectory_WhiskerTrialLite(d,varargin)
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
%                           this argument is not given, *all* '_WST.mat' files in directory
%                           'd' will be processed.
%
%           'ignore_files': Optional cell array of strings giving file name prefixes
%                           (i.e., file names without the '_WST.mat' suffix/extension) to ignore.
%                           Trumps 'include_files' argument.
%
%           'r_in_mm': The arc-length along whisker at which to measure kappa. Units of mm. Defaults to 1 mm.
%
%           'calc_forces': Either true or false. Requires the pole position
%                   to be tracked (i.e., barPos property of WhiskerSignalTrial must
%                   not be empty). Default is false. If true, will calculate the following timeseries:
%                       -M0:  Moment at the follicle. In Newton-meters.
%                       -Faxial: Axial force into follice. In Newtons.
%                       -deltaKappa: Change from baseline curvature, at point specified by r_point. In 1/mm.
%                       -Fnorm: The force on the whisker normal to the contacted object. In Newtons.
%                       -thetaAtBase: The whisker angle nearest the follicle. In degrees.
%                       -thetaAtContact: The whisker angle nearest the point of contact. I.e., nearest the center of the pole. In degrees.
%                       -distanceToPoleCenter: The closest distance between the whisker and the center of the pole. In mm.
%                       -meanKappa: The mean of kappa over the entire secondary polynomial fitted ROI. In 1/mm.
%
%   The following optional parameter/value pair arguments are ignored if 'calc_forces'
%   is not true:
%
%           'whisker_radius_at_base': Given in microns. Defaults is 33.5 microns.
%
%           'whisker_length': Given in mm. Default is 16 mm.
%
%           'youngs_modulus': In Pa. Default is 5e9 Pa.
%
%           'baseline_time_or_kappa_value': Either (1) a 1x2 vector giving starting and stopping times (inclusive) for measuring baseline whisker curvature, in sec;
%                                            or (2) a scaler giving a baseline kappa value (measured by the user separately) to directly subtract from kappa
%                                             timeseries, in 1/mm. Default is [0 0.1].
% NOTES:
%   Still need make these arguments settable on a whisker-by-whisker (trajectory ID by trajectory ID) basis.
%
%
%
%   DESCRIPTION:
%
%   Requires WhiskerSignalTrial objects to be saved, as .mat files, in the
%   directory specified by argument 'd'.  These files are read in one at a time and
%   converted to WhiskerTrialLite objects, which are then saved to disk in the same directory
%   as .mat files with a '_WL.mat' suffix/extension.
%
%   Processes all trajectory IDs within each WhiskerSignalTrial.
%
%
% 3/10, DHO.
%

p = inputParser;

p.addRequired('d', @ischar);
p.addParamValue('include_files', {}, @(x) all(cellfun(@ischar,x)));
p.addParamValue('ignore_files', {}, @(x) all(cellfun(@ischar,x)));
p.addParamValue('r_in_mm', 1, @(x) isnumeric(x) && numel(x)==1);
p.addParamValue('calc_forces', false, @islogical);

p.addParamValue('whisker_radius_at_base', 33.5, @isnumeric);
p.addParamValue('whisker_length', 16, @isnumeric);
p.addParamValue('youngs_modulus', 5e9, @isnumeric);
p.addParamValue('baseline_time_or_kappa_value', [0 0.1], @isnumeric);
p.addParamValue('proximity_threshold', -1, @isnumeric);

p.parse(d,varargin{:});

disp 'List of all arguments:'
disp(p.Results)


if ~strcmp(d(end), filesep)
    d = [d filesep];
end

currentDir = pwd;
cd(d)

fnall = arrayfun(@(x) x.name(1:(end-8)), dir([d '*_WST.mat']),'UniformOutput',false);

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

if ~isempty(fnall)
    if exist('parfor','builtin') % Parallel Computing Toolbox is installed
        parfor k=1:nfiles
            fn = fnall{k};
            disp(['Processing ''_WST.mat'' file '  fn ', ' int2str(k) ' of ' int2str(nfiles)])

            ws = pctload([fn '_WST.mat']);
            wl = Whisker.WhiskerTrialLiteI(ws,'r_in_mm',p.Results.r_in_mm,'calc_forces',p.Results.calc_forces,...
                'whisker_radius_at_base',p.Results.whisker_radius_at_base,...
                'whisker_length',p.Results.whisker_length,'youngs_modulus',p.Results.youngs_modulus,...
                'baseline_time_or_kappa_value',p.Results.baseline_time_or_kappa_value,...
                'proximity_threshold',p.Results.proximity_threshold);

            outfn = [fn '_WL.mat'];

            pctsave(outfn,wl);
        end
    else
        for k=1:nfiles
            fn = fnall{k};
            disp(['Processing ''_WST.mat'' file '  fn ', ' int2str(k) ' of ' int2str(nfiles)])

            load([fn '_WST.mat'],'ws');
            wl = Whisker.WhiskerTrialLiteI(ws,'r_in_mm',p.Results.r_in_mm,'calc_forces',p.Results.calc_forces,...
                'whisker_radius_at_base',p.Results.whisker_radius_at_base,...
                'whisker_length',p.Results.whisker_length,'youngs_modulus',p.Results.youngs_modulus,...
                'baseline_time_or_kappa_value',p.Results.baseline_time_or_kappa_value,...
                'proximity_threshold',p.Results.proximity_threshold);

            outfn = [fn '_WL.mat'];

            save(outfn,'wl');
        end
    end
end

cd(currentDir)
end

function pctsave(outfn,wl)
save(outfn,'wl');
end

function ws = pctload(loadfn)
load(loadfn,'ws');
end








