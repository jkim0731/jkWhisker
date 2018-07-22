function errorlist = jkmeasurements_dir(ppm, varargin)

if nargin > 2
    skip_flag = varargin{1};
else
    skip_flag = 'no';
end

save_filename = ['remeasure_', date, '.mat'];

if exist(save_filename,'file')
    if strcmp(skip_flag,'skip')    
        disp([save_filename ' already exists. Errorlist is returned without jkmeasurements.'])
        load(save_filename)
        return
    elseif strcmp(skip_flag,'no')
        disp('Overwriting remeasure file');
    else
        error('Not identified skip_flag argument in');
    end    
end

load_fn = ls('*follicle_n_mask.mat');
load(load_fn)
w = width;
h = height;
f = follicle_first;


flist = dir('*.measurements');
% v = VideoReader([flist(1).name(1:end-13),'.mp4']);
errorlist = zeros(size(flist,1),1); % List of files having error(s) in whisker tracking. (Consider re-tracking if the list is too long, about >10% of total trials in a session)
% 03/31/2017 JK

%% Listing error files
follicle_threshold = ppm*2; % 2 X 4 mm threshold for follicle
if contains(pwd,'piezo')       
    length_threshold = ppm*0.5; % 0.5 mm threshold for tracker data length
else
    length_threshold = ppm*3; % 3 mm threshold for tracker data length
end
    
parfor i = 1 : size(flist,1)
    error_flag = jkmeasurements(flist(i).name(1:end-13), w, h, f, follicle_threshold, length_threshold);
    if error_flag == 1
%         errorlist = [errorlist; str2double(flist(i).name(1:end-13))];
        errorlist(i) = str2double(flist(i).name(1:end-13));
    end
end

%%
prev_result = ls('remeasure_*.mat');
if ~isempty(prev_result)
    for i = 1: size(prev_result,1)
        delete(prev_result(i,:))
    end
end
save(save_filename,'errorlist') 