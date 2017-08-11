% close all
% clear
base_angle = 21;
behavior_base_dir = 'Y:\JK_temp\SoloData\';
whisker_base_dir = 'Y:\JK_temp\whisker\tracked\';
trial_types = {'rc', 'rf', 'lc', 'lf'};
mice = {'JK017','JK018','JK020'};
% sessionNum = {[1:4,6:15],[1,2,4:6,8:10],[2,4:6,8:19],[6,8:13],[2,4:8,13:19]};

sessions{1} = 1:10; % 017
sessions{2} = 1:10; % 018
sessions{3} = 1:10; % 020

mousensession = '13'; % defining which mouse, which session. (1) mouse, (2) before/after (3) session number within before/after
try
    mouseName = mice{str2double(mousensession(1))};
    sessionName = sprintf('S%02d', sessions{str2double(mousensession(1))}(str2double(mousensession(2))));
catch
    error('No matched mouse or session')
end

whisker_d = [whisker_base_dir mouseName sessionName '\'];
if ~exist('b','var') || ~iscell(b) || ~isfield(b{1},'mouseName') || ~strcmp(b{1}.mouseName,mouseName)
    load([behavior_base_dir mouseName filesep 'behavior.mat']) % load b
end
%%
b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
b_session = b{b_ind};
tt_ind = cell(1,length(trial_types));
wl = cell(1,length(trial_types));

for ti = 1 : length(trial_types)    
    tt_ind{ti} = find(cellfun(@(x) strcmp(x.trialType,trial_types{ti}),b_session.trials));
    temp_files = cell(length(tt_ind{ti}),1);
    for j = 1 : length(tt_ind{ti})
        temp_files{j} = num2str(tt_ind{ti}(j));
    end
    wl{ti} = Whisker.WhiskerTrialLiteArray_2pad(whisker_d,'include_files',temp_files);     
end