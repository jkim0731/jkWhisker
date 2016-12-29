function avi2tiff_convertDirectory(d)
%
%  d: directory path name as string.
%   
%   Converts each AVI file in directory d to a TIFF
%       file with same name (but .tiff extension).
%
%   Will not overwrite existing TIFF files.  These files are skipped.
%
%
% DHO, 12/08.
%

% d='C:\Documents and Settings\dan\My Documents\work\full\';

if ~strcmp(d(end), '\')
    d = [d '\'];
end

currentDir = pwd;
cd(d) 

fnall = dir([d '*avi']);

tiffnames = dir([d '*tif']);
tiffnames = arrayfun(@(x) x.name, tiffnames, 'UniformOutput',false);

nfiles = length(fnall);

if ~isempty(fnall)
    for k=1:nfiles
        disp(['Processing file ' int2str(k) ' of ' int2str(nfiles)]) 
        fn = fnall(k).name;
        outfn = [fn(1:(end-3)) 'tif'];
        if ~isempty(tiffnames)
            if sum(strcmp(outfn, tiffnames)) < 1
                Whisker.avi2tiff(fn);
            %------------------------------
            % --- DELETE ORIGINAL AVI FILE:
            delete(fn)
            %------------------------------
            end
        else
            Whisker.avi2tiff(fn);
            %------------------------------
            % --- DELETE ORIGINAL AVI FILE:
            delete(fn)
            %------------------------------
        end
    end
end

cd(currentDir)






