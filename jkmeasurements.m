function error = jkmeasurements(fn,width, height, follicle_first,follicle_threshold, length_threshold)
error = 0; % indicates if there was an error in whisker tracker or not. Used to tag files to exclude from further analysis (if the list is too long, probably need to re-track whiskers)
% 03/31/2017 JK
%%
fn = strsplit(fn,' ');
fn = fn{1};
if length(fn) < 13
    fn = [fn, '.measurements'];
end
if ~strcmp(fn(end-12:end), '.measurements')
    fn = [fn, '.measurements'];
end

b = LoadMeasurements(fn);


%% Setting some constraints from the first frame
for i = 1 : length(b)
    if b(i).tip_y > b(i).follicle_y
        temp = b(i).tip_y;
        b(i).tip_y = b(i).follicle_y;
        b(i).follicle_y = temp;
    end
    
    if abs(b(i).follicle_y - follicle_first(1,1)) < follicle_threshold && b(i).follicle_x > width/2 && b(i).length > length_threshold % top-view whisker
        b(i).label = 0;
    elseif abs(b(i).follicle_y - follicle_first(2,1)) < follicle_threshold && b(i).follicle_x < width/2 && b(i).length > length_threshold % front-view whisker
        b(i).label = 1;
    else
        b(i).label = -1;
    end
end

fids = zeros(length(b),1);
labels = zeros(length(b),1);
for i = 1 : length(b)
    fids(i) = b(i).fid;
    labels(i) = b(i).label;
end

% making exactly one top-view and front-view whisker tracked.
for i = 1 : max(fids)+1
    tind = find(fids == i-1);
    n0ind = find(labels(tind) == 0); % top-view whisker, to have tid 0
    n1ind = find(labels(tind) == 1); % front-view, to have tid 1
        
    if length(n0ind) > 1 % when more than one whisker is found to be top-view whisker, take the longest one.
        t0max = 0;        
        maxind = 0;
        for j = 1 : length(n0ind)
            if b(tind(n0ind(j))).length > t0max
                t0max = b(tind(n0ind(j))).length;
                maxind = j;
            end
        end
        for j = 1 : length(n0ind)
            if j ~= maxind
                b(tind(n0ind(j))).label = -1;
            end
        end
    elseif isempty(n0ind) % when a top-view whisker is not found, try searching again, with 2x follicle y movement
        t0max = 0;
        maxind = 0;
        for j = 1 : length(tind)
            if abs(b(tind(j)).follicle_y - follicle_first(1,1)) < follicle_threshold * 2 && b(tind(j)).follicle_x > width/2 && b(tind(j)).length > t0max
                t0max = b(tind(j)).length;
                maxind = j;
            end
        end
        if maxind == 0
            fprintf('no top-view whisker found at file %s frame % d \n', fn(1:end-11), i)
            error = 1;
        else
            b(tind(maxind)).label = 0;
        end
    end

    if length(n1ind) > 1 % same as in n0ind
        t1max = 0;
        maxind = 0;
        for j = 1 : length(n1ind)
            if b(tind(n1ind(j))).length > t1max
                t1max = b(tind(n1ind(j))).length;
                maxind = j;
            end
        end
        for j = 1 : length(n1ind)
            if j ~= maxind
                b(tind(n1ind(j))).label = -1;
            end
        end
    elseif isempty(n1ind)
        t0max = 0;
        maxind = 0;
        for j = 1 : length(tind)
            if abs(b(tind(j)).follicle_y - follicle_first(2,1)) < follicle_threshold * 2 && b(tind(j)).follicle_x < width/2 && b(tind(j)).length > t0max
                t0max = b(tind(j)).length;
                maxind = j;
            end
        end
        if maxind == 0
            fprintf('no front-view whisker found at file %s frame % d \n', fn(1:end-11), i)
            error = 1;
        else
            b(tind(maxind)).label = 1;
        end
    end

end


i = 1 ;
while (1)
    old_fn = [fn(1:end-13),'_old',str2double(i),'.measurements'];
    if exist(old_fn,'file')
        i = i + 1;
    else
        FileRename(fn,old_fn);
        break
    end
end

SaveMeasurements(fn,b)

