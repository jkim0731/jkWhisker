function errorlist = jkmeasurements_dir()

save_filename = ['remeasure_', date, '.mat'];
if exist(save_filename,'file')
    disp([save_filename ' already exists. Errorlist is returned without jkmeasurements.'])
    load(save_filename)
    return
end

load_fn = ls('*follicle_n_mask.mat');
load(load_fn)
w = width;
h = height;
f = follicle_first;

% Adjust this accordingly
follicle_threshold = 20;
length_threshold = 30;



flist = dir('*.measurements');
% v = VideoReader([flist(1).name(1:end-13),'.mp4']);
errorlist = zeros(size(flist,1),1); % List of files having error(s) in whisker tracking. (Consider re-tracking if the list is too long, about >10% of total trials in a session)
% 03/31/2017 JK

%% Listing error files

parfor i = 1 : size(flist,1)
    error = jkmeasurements(flist(i).name(1:end-13), w, h, f, follicle_threshold, length_threshold);
    if error == 1
%         errorlist = [errorlist; str2double(flist(i).name(1:end-13))];
        errorlist(i) = str2double(flist(i).name(1:end-13));
    end
end

%%

save(save_filename,'errorlist') 