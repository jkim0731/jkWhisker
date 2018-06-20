function buildWL(whisker_d, b_session, rInMm)

cd(whisker_d)
wstFnList = dir('*_WST.mat');
wstTrialFn = cell(length(wstFnList),1);
for i = 1 : length(wstFnList)
    wstTrialFn{i} = wstFnList(i).name(1:end-8);
end
loadfn = dir('*_touch_hp.mat');
load(loadfn(1).name)
Whisker.makeAllDirectory_WhiskerTrialLite_2pad(whisker_d, 'include_files', wstTrialFn, 'calc_forces', false, 'rInMm', rInMm, ...
    'behavior', b_session, 'hp_peaks', hp_peaks, 'touch_hp', touch_hp, 'psi1', psi1, 'psi2', psi2, 'touch_boundary_thickness', 2, ...
    'servo_distance_pair', servo_distance_pair)