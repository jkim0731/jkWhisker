function jkmeasurements_single_whisker(fn)
% To find out the maximum length whisker from .measurements,
% especially when 'classify' was not commanded

if isempty(strfind(fn,'.measurements'))
    fn = [fn, '.measurements'];
end
b = LoadMeasurements(fn);

%% Finding the longest whisker
fids = zeros(length(b),1);
for i = 1 : length(b)
    fids(i) = b(i).fid;
end

% making exactly one top-view and front-view whisker tracked.
for i = 1 : max(fids)+1
    tind = find(fids == i-1);            
    t0max = 0;
    maxind = 0;
    for j = 1 : length(tind)
        b(tind(j)).label = -1;
        if b(tind(j)).length > t0max
            t0max = b(tind(j)).length;
            maxind = j;
        end
    end
    if maxind == 0
        error('no top-view whisker found at file %s frame % d \n', fn(1:end-11), i)
    else        
        b(tind(maxind)).label = 0;
    end
end

SaveMeasurements(fn,b)

