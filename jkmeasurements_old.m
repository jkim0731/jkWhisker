function jkmeasurements(fn)
%%
fn = '2';
v = VideoReader([fn,'.mp4']);

follicle_roi = [10 5]; % 10 & 5 pixels movement of follicle from the first frame for x and y, respectively.
follicle_first = zeros(2,2);
vv = read(v,1);
figure, imshow(vv), axis off, axis image, hold all
i = 1;
while (i < 3)
    [y, x] = ginput(1);
    scatter(y, x, 'mo');
    if i == 1
        button = questdlg('is this correct?', 'Top-view follicle', 'Yes', 'No', 'Cancel', 'Yes');
        switch button
            case 'Yes'
                follicle_first(i,2)= y; follicle_first(i,1) = x;
                i = i + 1;
            case 'No' 
                continue
            case 'Cancel'
                disp('measurements adjustment aborted')
                return
        end
    elseif i == 2
        button = questdlg('is this correct?', 'Front-view follicle', 'Yes', 'No', 'Cancel', 'Yes');
        switch button
            case 'Yes'
                follicle_first(i,2)= y; follicle_first(i,1) = x;
                i = i + 1;
            case 'No' 
                continue
            case 'Cancel'
                disp('measurements adjustment aborted')
                return
        end
    end
end

%%
follicle_y_threshold = v.height/3*2;
follicle_x_threshold = v.width/2;
length_threshold = v.height/5;


fn = strsplit(fn,' ');
fn = fn{1};
if length(fn) < 13
    fn = [fn, '.measurements'];
end
if ~strcmp(fn(end-12:end), '.measurements')
    fn = [fn, '.measurements'];
end

b = LoadMeasurements(fn);

fids = zeros(length(b),1);
labels = zeros(length(b),1);
for i = 1 : length(b)
    fids(i) = b(i).fid;
    labels(i) = b(i).label;
end

%% Assume that first frame (fids = 0) always gives correct answer
% tind = find(fids == 0);
% for i = 1 : length(tind)
%     if b(i).tip_y > b(i).follicle_y
%         temp = b(i).tip_y;
%         b(i).tip_y = b(i).follicle_y;
%         b(i).follicle_y = temp;
%     end
%     if b(i).length < 80
%         b(i).label = -1;
%     else
%         b(i).label = 0;
%     end
%     if b(i).follicle_y < 300
%         b(i).label = -1;
%     end
%     if b(i).label == 0 && b(i).follicle_x < 350
%         b(i).label = 1;
%     end
% end
% 
% for i = length(tind)
%     n0ind = find(labels(tind) == 0);
%     n1ind = find(labels(tind) == 1);
%     if length(n0ind) > 1
%         t0max = b(tind(n0ind(1))).length;
%         maxind = 1;
%         for j = 2 : length(n0ind)
%             if b(tind(n0ind(j))).length > t0max
%                 t0max = b(tind(n0ind(j))).length;
%                 maxind = j;
%             end
%         end
%         for j = 1 : length(n0ind)
%             if j ~= maxind
%                 b(tind(n0ind(j))).label = -1;
%             end
%         end
%     elseif isempty(n0ind)
%         t0max = 0;
%         maxind = 0;
%         for j = 1 : length(tind)
%             if b(tind(j)).follicle_y > 300 && b(tind(j)).follicle_x > 350 && b(tind(j)).length > t0max
%                 t0max = b(tind(j)).length;
%                 maxind = j;
%             end
%         end
%         if maxind == 0
%             fprintf('no top-view whisker found at file %s frame % d', fn(1:end-11), i+1)
%         end
%         b(tind(maxind)).label = 0;
%     elseif length(n0ind) == 1
%         maxind = n0ind;
%     end
%     
%     top_ref_ind = tind(maxind);
%     
%     
%     if length(n1ind) > 1
%         t1max = b(tind(n1ind(1))).length;
%         maxind = 1;
%         for j = 2 : length(n1ind)
%             if b(tind(n1ind(j))).length > t1max
%                 t1max = b(tind(n1ind(j))).length;
%                 maxind = j;
%             end
%         end
%         for j = 1 : length(n1ind)
%             if j ~= maxind
%                 b(tind(n1ind(j))).label = -1;
%             end
%         end
%     elseif isempty(n1ind)
%         t0max = 0;
%         maxind = 0;
%         for j = 1 : length(tind)
%             if b(tind(j)).follicle_y > 300 && b(tind(j)).follicle_x < 350 && b(tind(j)).length > t0max
%                 t0max = b(tind(j)).length;
%                 maxind = j;
%             end
%         end
%         if maxind == 0
%             fprintf('no front-view whisker found at file %s frame % d', fn(1:end-11), i+1)
%         end
% 
%         b(tind(maxind)).label = 1;
%     elseif length(n0ind) == 1
%         maxind = n0ind;        
%     end
%     
%     
% end
% 
% 
%% Setting some constraints from the first frame
for i = 1 : length(b)
    if b(i).tip_y > b(i).follicle_y
        temp = b(i).tip_y;
        b(i).tip_y = b(i).follicle_y;
        b(i).follicle_y = temp;
    end
    if b(i).length < length_threshold
        b(i).label = -1;
    else
        b(i).label = 0;
    end
    if b(i).follicle_y < follicle_y_threshold
        b(i).label = -1;
    end
    if b(i).label == 0 && b(i).follicle_x < follicle_x_threshold
        b(i).label = 1;
    end
end

% making exactly one top-view and front-view whisker tracked.
for i = 1 : max(fids)+1
    tind = find(fids == i-1);
    n0ind = find(labels(tind) == 0); % top-view whisker, to have tid 0
    n1ind = find(labels(tind) == 1); % front-view, to have tid 1
        
    if length(n0ind) > 1
        t0max = b(tind(n0ind(1))).length;
        maxind = 1;
        for j = 2 : length(n0ind)
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
    elseif isempty(n0ind)
        t0max = 0;
        maxind = 0;
        for j = 1 : length(tind)
            if b(tind(j)).follicle_y > follicle_y_threshold && b(tind(j)).follicle_x > follicle_x_threshold && b(tind(j)).length > t0max
                t0max = b(tind(j)).length;
                maxind = j;
            end
        end
        if maxind == 0
            fprintf('no top-view whisker found at file %s frame % d \n', fn(1:end-11), i)
        else
            b(tind(maxind)).label = 0;
        end
    end
        
    if length(n1ind) > 1
        t1max = b(tind(n1ind(1))).length;
        maxind = 1;
        for j = 2 : length(n1ind)
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
            if b(tind(j)).follicle_y > follicle_y_threshold && b(tind(j)).follicle_x < follicle_x_threshold && b(tind(j)).length > t0max
                t0max = b(tind(j)).length;
                maxind = j;
            end
        end
        if maxind == 0
            fprintf('no front-view whisker found at file %s frame % d \n', fn(1:end-11), i)
        else
            b(tind(maxind)).label = 1;
        end
    end
end

SaveMeasurements(fn,b)
save('',)
