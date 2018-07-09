function follicle_n_mask(mouseName,sessionName,videoloc, ppm, transmm, facePosition, varargin)

%% Setup whisker array builder 
% mouseName = 'AH0653'
% sessionName ='S03'
% videoloc = 'JK'
% facePosition = 'bottom'
% optional = 'noskip'
optional = 'skip';
if nargin > 6
    optional = varargin{1};
elseif nargin > 7
    error('Too many input arguments')
end

transpix = round(ppm * transmm);

% if exist('optional','var')
%     d = ([videoloc filesep mouseName sessionName filesep optional filesep])
% else
    d = ([videoloc filesep mouseName sessionName filesep])
    currD = pwd;
% end
% load(['Z:\Users\Jon\DATA\BehaviorArrays\solo_' mouseName '_' sessionName '.mat'])

cd(d)

maskfn = [mouseName sessionName 'follicle_n_mask.mat'];
if exist(maskfn,'file') && strcmp(optional, 'skip')
    disp([maskfn ' already exists. Skipping.'])    
    return
elseif exist(maskfn,'file') && strcmp(optional, 'noskip')
    disp([maskfn ' already exists. Overriding.'])    
end

if contains(sessionName, 'spont')
    number_of_random_trials = 1; % for averaging for mask detection
elseif contains(sessionName, 'piezo')
    number_of_random_trials = 3; % for averaging for mask detection
else
    number_of_random_trials = 10; % for averaging for mask detection
end

%% Follicle
flist = dir('*.mp4');
center_ind = round(length(flist)/2);
v = VideoReader(flist(center_ind).name);
% length_threshold = 40;
% follicle_threshold = 40; % 40 pixels movement of follicle in x and y is tolerable
follicle_first = zeros(2,2);
width = v.Width;
height = v.Height;
vavg = zeros(v.Height,v.Width);
nof = fix(v.FrameRate*v.Duration); % number of frames
while hasFrame(v)
    vtemp = readFrame(v);    
    vtemp = double(vtemp(:,:,1));
    vavg = vavg + vtemp/nof;
end
vavg = mat2gray(vavg);
figure('units','normalized','outerposition',[0 0 1 1]), imshow(vavg,'InitialMagnification','fit'), axis off, axis image, title({[mouseName, ' ', sessionName];'Select follicle points, first top-view, and then front-view'}), hold all
i = 1;
while (i < 3)
    [y, x] = ginput(1);
    obj_h = scatter(y, x, 'mo');
    if i == 1
        button = questdlg('is this correct?', 'Top-view follicle', 'Yes', 'No', 'Cancel', 'Yes');
    else
        button = questdlg('is this correct?', 'Front-view follicle', 'Yes', 'No', 'Cancel', 'Yes');
    end
    switch button
        case 'Yes'
            follicle_first(i,2)= y; follicle_first(i,1) = x;
            i = i + 1;
            scatter(y, x, 'go');
        case 'No' 
            delete(obj_h);
            continue
        case 'Cancel'
            disp('measurements adjustment aborted')
            return
    end
end

drawnow;

%% Mask
temp_vavg = zeros(v.Height,v.Width, number_of_random_trials);
randlist = randi(length(flist), 1, number_of_random_trials);
parfor i = 1 : number_of_random_trials
    v = VideoReader(flist(randlist(i)).name);    
    nof = fix(v.FrameRate*v.Duration); % number of frames
    while hasFrame(v)
        vtemp = readFrame(v);    
        vtemp = double(vtemp(:,:,1));
        temp_vavg(:,:,i) = temp_vavg(:,:,i) + vtemp/nof;
    end
    temp_vavg(:,:,i) = imgaussfilt(temp_vavg(:,:,i),3);    
%     width = min(size(vavg,2),size(temp_vavg,2)); % for error in JK027 S17.
%     vavg = vavg(:,1:width) + temp_vavg(:,1:width)/number_of_random_trials; % for error in JK027 S17.
end
vavg = mean(temp_vavg,3);
vavg_filt = imgaussfilt(vavg,3);
maskx = {[],[]};
masky = {[],[]};

BW = edge(vavg_filt);
switch facePosition
    case 'bottom'
        BW = circshift(BW, -transpix, 1);
    case 'top'
        BW = circshift(BW, transpix, 1);
    case 'right'
        BW = circshift(BW, -transpix, 2);
    case 'left'
        BW = circshift(BW, transpix, 2);
end

[edge_i,edge_j] = ind2sub(size(BW),find(BW));
top_ind = find(edge_j >= size(BW,2)/2);
edge_i = edge_i - size(BW,2) + size(vavg,2);
edge_j(top_ind) = edge_j(top_ind) - (size(BW,1) - size(vavg,1))*(sind(21)/sind(45)); % 21 degrees tilted
scatter(edge_j,edge_i,3,'mo');

i = 1;
while (i < 3)            
    temp_ind = zeros(2,1);
    for j = 1 : 2
        [x, y] = ginput(1);
        temp_dist = (edge_j - x).^2 + (edge_i - y).^2;
        [~,temp_ind(j)] = min(temp_dist);        
    end
    temp_j = floor(edge_j(min(temp_ind):max(temp_ind)));
    temp_j(temp_j < 1) = 1;
    temp_j(temp_j > size(vavg,2)) = size(vavg,2);    
    temp_i = floor(edge_i(min(temp_ind):max(temp_ind)));
    temp_i(temp_i < 1) = 1;
    temp_i(temp_i > size(vavg,1)) = size(vavg,1);
    temp_bw = zeros(size(vavg));
    temp_ind = sub2ind(size(temp_bw),temp_i,temp_j);
    temp_bw(temp_ind) = 1;
    bl = bwlabel(temp_bw);        
    [mask_i,mask_j] = ind2sub(size(vavg),find(bl == bl(temp_i(1),temp_j(1))));

    obj_h = scatter(mask_j,mask_i,3,'bo');
    
    qnum = length(mask_j);
    polyDegree = 2;
    mask_j = mask_j'; mask_i = mask_i';
    q = (0:(qnum-1))./(qnum-1);
    px = Whisker.polyfit(q,mask_j,polyDegree);
    py = Whisker.polyfit(q,mask_i,polyDegree);
    q = linspace(-0.3,1.3);
    maskx{i} = polyval(px,q);
    masky{i} = polyval(py,q);
    plot_h = plot(maskx{i},masky{i},'g-','LineWidth',2);
    
    drawnow;
    if i == 1
        button = questdlg('is this correct?', 'Top-view mask', 'Yes', 'No', 'Cancel', 'Yes');
    else
        button = questdlg('is this correct?', 'Front-view mask', 'Yes', 'No', 'Cancel', 'Yes');
    end
    switch button
        case 'Yes'
            i = i + 1;
        case 'No'         
            answer = questdlg('Do you want manual point selection?', 'Manual selection', 'Yes', 'No', 'Cancel', 'Yes');
            switch answer
                case 'Yes'
                    delete(obj_h);
                    delete(plot_h);
                    maskPoints = []; % points of the polygon                    
                    temp_point = ginput(1);
                    j = 1;
                    while(~isempty(temp_point)) % finish drawing polygon by pressing "enter"
                        maskPoints = [maskPoints; temp_point];
                        plot(maskPoints(j,1), maskPoints(j,2), 'bo', 'MarkerSize', 3)                        
                        temp_point = ginput(1);
                        j = j + 1;
                    end
                    mask_i = maskPoints(:,1); mask_j = maskPoints(:,2);
                    qnum = length(mask_j);
                    polyDegree = 2;
                    mask_j = mask_j'; mask_i = mask_i';
                    q = (0:(qnum-1))./(qnum-1);
                    px = Whisker.polyfit(q,mask_i,polyDegree);
                    py = Whisker.polyfit(q,mask_j,polyDegree);
                    q = linspace(-0.3,1.3);
                    maskx{i} = polyval(px,q);
                    masky{i} = polyval(py,q);
                    plot_h = plot(maskx{i},masky{i},'g-','LineWidth',2);

                    drawnow;
                    if i == 1
                        answer2 = questdlg('is this correct?', 'Top-view mask', 'Yes', 'No', 'Cancel', 'Yes');
                    else
                        answer2 = questdlg('is this correct?', 'Front-view mask', 'Yes', 'No', 'Cancel', 'Yes');
                    end
                    switch answer2
                        case 'Yes'
                            i = i + 1;
                        case 'No' 
                            delete(obj_h);
                            delete(plot_h);
                            continue
                        case 'Cancel'
                            disp('measurements adjustment aborted')
                        return
                    end
                case 'No'
                    delete(obj_h);
                    delete(plot_h);
                    continue
                case 'Cancel'
                    disp('measurements adjustment aborted')
                    return
            end
        case 'Cancel'
            disp('measurements adjustment aborted')
                    return
    end
end

%% save mask

save(maskfn,'maskx','masky','width', 'height', 'follicle_first')
cd(currD)
