function follicle_n_mask(mouseName,sessionName,videoloc,varargin)

%% Setup whisker array builder 
% mouseName = 'AH0653'
% sessionName ='S03'
% videoloc = 'JK'
% optional = 'Spont'
if nargin > 3
    optional = varargin{4};
elseif nargin > 4
    error('Too many input arguments')
end

if exist('optional','var')
    d = (['Z:\Data\Video\' videoloc filesep mouseName sessionName filesep optional filesep])
else
    d = (['Z:\Data\Video\' videoloc filesep mouseName sessionName filesep])
end
% load(['Z:\Users\Jon\DATA\BehaviorArrays\solo_' mouseName '_' sessionName '.mat'])

cd(d)

maskfn = [mouseName sessionName 'follicle_n_mask.mat'];
if exist(maskfn,'file')
    disp([maskfn ' already exists. Skipping.'])
    return
end

%% Follicle
flist = dir('*.mp4');
v = VideoReader(flist(1).name);
% length_threshold = 40;
% follicle_threshold = 40; % 40 pixels movement of follicle in x and y is tolerable
follicle_first = zeros(2,2);
vv = read(v,1);
width = size(vv,2);
height = size(vv,1);
vavg = zeros(height,width);
for i = 1 : v.NumberOfFrames
    vtemp = read(v,i);    
    vtemp = double(vtemp(:,:,1));
    vavg = vavg + vtemp/v.NumberOfFrames;
end
vavg = mat2gray(vavg);
figure('units','normalized','outerposition',[0 0 1 1]), imshow(vavg), axis off, axis image, title('Select follicle points, first top-view, and then front-view'), hold all
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
                scatter(y, x, 'go');
            case 'No' 
                delete(findobj(gca, 'type', 'patch'));
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
                delete(findobj(gca, 'type', 'patch'));
                continue
            case 'Cancel'
                disp('measurements adjustment aborted')
                return
        end
    end
end

drawnow;

%% Mask
button = 1;
masknum = str2num(cell2mat(inputdlg({'How many trajectories?','How many points?'},'Trajectories',1,{'2','4'})));
maskx = zeros(masknum(1),masknum(2));
masky = zeros(masknum(1),masknum(2));
i = 1; 
while (i <= masknum(1))
    j = 1;
    while (j <= masknum(2))
        [x, y, button] = ginput(1);
        x = round(x); y = round(y);        
        scatter(x,y,'mo');
        maskx(i,j) = x;
        masky(i,j) = y;
        j = j + 1;
    end
    if j > 1            
        qnum = length(maskx(i,:));
        if qnum < 2
            error('Must define at least 2 points.')
        elseif qnum < 6
            polyDegree = qnum-1;
        else
            polyDegree = 5;
        end
        q = (0:(qnum-1))./(qnum-1);
        px = Whisker.polyfit(q,maskx(i,:),polyDegree);
        py = Whisker.polyfit(q,masky(i,:),polyDegree);
        q = linspace(0,1);
        x = polyval(px,q);
        y = polyval(py,q);
        plot(x,y,'g-','LineWidth',2)
    end

    button = questdlg('is this correct?', 'Mask points', 'Yes', 'No', 'Cancel', 'Yes');
    switch button
        case 'Yes'            
            i = i + 1;            
        case 'No' 
            delete(findobj(gca, 'type', 'patch'));
            delete(findobj(gca, 'type', 'line'));
            continue
        case 'Cancel'
            disp('masking aborted')
            return
    end
end
hold off;
drawnow;

%% save mask

save(maskfn,'maskx','masky','width', 'height', 'follicle_first')