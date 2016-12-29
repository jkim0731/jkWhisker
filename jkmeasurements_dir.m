function jkmeasurements_dir()

flist = dir('*.measurements');
% v = VideoReader([flist(1).name(1:end-13),'.mp4']);
v = VideoReader('2.mp4');
length_threshold = 40;

follicle_threshold = 40; % 40 pixels movement of follicle in x and y is tolerable
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
figure, imshow(vavg), axis off, axis image, hold all
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

%%
parfor i = 1 : size(flist,1)
% for i = 1
    jkmeasurements(flist(i).name(1:end-13), width, height, follicle_first, follicle_threshold, length_threshold);
end
save_filename = ['re_adj_jk_', date, '.mat'];
save(save_filename,'follicle_first')