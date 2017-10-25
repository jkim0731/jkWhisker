function jkmeasurements_single_whisker_dir(varargin)
if nargin < 1
    flist = dir('*.measurements');
else
    if isstring(varargin{1})               
        if varargin{1}(end) == filesep
            flist = dir([varargin{1},'*.measurements']);
        else
            flist = dir([varargin{1},filesep,'*.measurements']);
        end
    else
        error('Wrong input - it should be a directory name');
    end
end

for i = 1 : length(flist)
    jkmeasurements_single_whisker(flist(i).name);
end