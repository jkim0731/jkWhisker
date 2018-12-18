function buildWL_2pad(whisker_d, rInMm, varargin)

p = inputParser;

p.addRequired('whisker_d', @ischar);
p.addRequired('rInMm', @isnumeric);
p.addParameter('b_session', [], @(x) isa(x,'Solo.BehavTrial2padArray'));
p.addParameter('touchKappaSTDthreshold', 2, @isnumeric);
p.addParameter('whiskingAmpThreshold', 2.5, @isnumeric);
p.addParameter('stdHistogramThreshold', 3, @isnumeric);
p.addParameter('distanceHistogramBin', 0.2, @isnumeric);
p.addParameter('touchBoundaryThickness', 0.25, @isnumeric);
p.addParameter('touchBoundaryBuffer', 0.1, @isnumeric);
p.addParameter('maxPointsNearHyperplane', 10, @isnumeric);
            
p.parse(whisker_d, rInMm,varargin{:});

cd(whisker_d)
wstFnList = dir('*_WST.mat');
wstTrialFn = cell(length(wstFnList),1);
for i = 1 : length(wstFnList)
    wstTrialFn{i} = wstFnList(i).name(1:end-8);
end


% wlFnList = dir('*_WL_2pad.mat');
% wlTrialFn = cell(length(wlFnList),1);
% for i = 1 : length(wlFnList)
%     wlTrialFn{i} = wlFnList(i).name(1:end-12);
% end
% fnind = find(~ismember(wstTrialFn, wlTrialFn));
% wstTrialFn = wstTrialFn(fnind);

load(wstFnList(1).name)

loadfn = [ws.mouseName, ws.sessionName, '_touch_hp.mat'];
load(loadfn)

if ~isempty(p.Results.b_session)
    Whisker.makeAllDirectory_WhiskerTrialLite_2pad(whisker_d, 'include_files', wstTrialFn, 'calc_forces', false, 'rInMm', rInMm, ...
        'behavior', p.Results.b_session, 'hp_peaks', hp_peaks, 'touch_hp', touch_hp, 'psi1', psi1, 'psi2', psi2, ...
        'servo_distance_pair', servo_distance_pair, 'touchKappaSTDthreshold', p.Results.touchKappaSTDthreshold, 'whiskingAmpThreshold', p.Results.whiskingAmpThreshold, 'stdHistogramThreshold', p.Results.stdHistogramThreshold, ...
        'distanceHistogramBin', p.Results.distanceHistogramBin, 'touchBoundaryThickness', p.Results.touchBoundaryThickness, 'touchBoundaryBuffer', p.Results.touchBoundaryBuffer, 'maxPointsNearHyperplane', p.Results.maxPointsNearHyperplane);
else
    Whisker.makeAllDirectory_WhiskerTrialLite_2pad(whisker_d, 'include_files', wstTrialFn, 'calc_forces', false, 'rInMm', rInMm);
end