function [pole_edge, varargout] = pole_edge_detection(video_fn)

% Automatic detection of pole_edge from video. Both front and top view.
% Currently only for 2pad 2017 JK

% Added automatic pole_available_frames 2018/03/06 JK

%% Manually set factors. Should be changed individually for different settings
width_factor = 1;
height_factor = 0.7;
%%
if isnumeric(video_fn)
    v = VideoReader([num2str(video_fn),'.mp4']);
else
    v = VideoReader([video_fn,'.mp4']);
end
nof = fix(v.FrameRate*v.Duration);
h = floor(v.height*height_factor);
w =  floor(v.width*width_factor);
vavg = zeros(v.height,v.width);
pole_img_ts = zeros(h,w,nof);
for i = 1 : nof
    temp = readFrame(v);
    vavg = vavg + double(temp)/nof;    
    pole_img_ts(:,:,i) = adapthisteq(temp(1:h,1:w,1),'NumTiles',[5,20]);    
end
vavg = mat2gray(vavg);
%%    
pole_edge_pic = edge(vavg,'Roberts');

front_pole_edge = pole_edge_pic; front_pole_edge(floor(v.height*0.7):end,:) = 0; front_pole_edge(:, floor(v.width*0.5):end) = 0; front_pole_edge(:,1:10) = 0;    
[edge_i,edge_j] = find(front_pole_edge == 1);  
result_i = zeros(max(edge_j)-min(edge_j)+1,1);
result_j = zeros(size(result_i));
for i = min(edge_j):max(edge_j) % just take the edge detected at the lower part, i.e., closer to the face
    temp_samej_ind = find(edge_j == i);
    if ~isempty(temp_samej_ind)
        inner_edge_ind = find(edge_i(temp_samej_ind) == max(edge_i(temp_samej_ind)));        
        result_i(i-min(edge_j)+1) = edge_i(inner_edge_ind + min(temp_samej_ind)-1);
        result_j(i-min(edge_j)+1) = edge_j(inner_edge_ind + min(temp_samej_ind)-1);
    end
end
noedge_ind = find(result_i == 0);
if ~isempty(noedge_ind) % connect them if there is a nick
    edge_ind = setdiff(1:length(result_i),noedge_ind);
    itp_result_i = interp1(edge_ind, result_i(edge_ind),1:length(result_i));
    itp_result_j = interp1(edge_ind, result_j(edge_ind),1:length(result_j));
    result_i = round(itp_result_i);
    result_j = round(itp_result_j);
end

front_pole_edge = zeros(size(front_pole_edge));
front_pole_edge(sub2ind(size(front_pole_edge),result_i,result_j)) = 1;

CC = bwconncomp(front_pole_edge); % take the longest (or largest) connected component
ccc = cellfun(@(x) length(x), CC.PixelIdxList);
ci = find(ccc == max(ccc));
C = CC.PixelIdxList{ci};
[result_i, result_j] = ind2sub(size(front_pole_edge),C);

[result_j, i_sort] = sort(result_j);
result_i = result_i(i_sort);

if length(result_j) > v.width*0.2 % linear fit of the front pole
    p = polyfit(result_j(1:floor(v.width*0.2)),result_i(1:floor(v.width*0.2)),1);
else
    p = polyfit(result_j,result_i,1);
end

q = linspace(1,v.width*0.5);
pole_edge{2} = [result_i, result_j];
varargout{1}{2} = [polyval(p,q); q];

%%
top_pole_edge = pole_edge_pic; top_pole_edge(floor(v.height*0.7):end,:) = 0; top_pole_edge(:, 1:floor(v.width*0.5)) = 0;
[edge_i,edge_j] = find(top_pole_edge == 1);   
result_i = zeros(max(edge_j)-min(edge_j)+1,1);
result_j = zeros(size(result_i));
for i = min(edge_j):max(edge_j)
    temp_samej_ind = find(edge_j == i);
    if ~isempty(temp_samej_ind)
        inner_edge_ind = find(edge_i(temp_samej_ind) == max(edge_i(temp_samej_ind)));        
        result_i(i-min(edge_j)+1) = edge_i(inner_edge_ind + min(temp_samej_ind)-1);
        result_j(i-min(edge_j)+1) = edge_j(inner_edge_ind + min(temp_samej_ind)-1);
    end
end
noedge_ind = find(result_i == 0);
if ~isempty(noedge_ind)
    edge_ind = setdiff(1:length(result_i),noedge_ind);
    itp_result_i = interp1(edge_ind, result_i(edge_ind),1:length(result_i));
    itp_result_j = interp1(edge_ind, result_j(edge_ind),1:length(result_j));
    result_i = round(itp_result_i);
    result_j = round(itp_result_j);
end

top_pole_edge = zeros(size(top_pole_edge));
top_pole_edge(sub2ind(size(top_pole_edge),result_i,result_j)) = 1;    

CC = bwconncomp(top_pole_edge);
ccc = cellfun(@(x) length(x), CC.PixelIdxList);
ci = find(ccc == max(ccc));
C = CC.PixelIdxList{ci};
[result_i, result_j] = ind2sub(size(front_pole_edge),C);

[result_j, i_sort] = sort(result_j, 'descend');
result_i = result_i(i_sort);

if length(result_j) > v.width*0.2 % linear fit of the front pole
    p = polyfit(result_j(1:floor(v.width*0.2)),result_i(1:floor(v.width*0.2)),1);
else
    p = polyfit(result_j,result_i,1);
end

q = linspace(v.width,v.width*0.5);
pole_edge{1} = [result_i, result_j];
varargout{1}{1} = [polyval(p,q); q];
varargout{2} = vavg; 

%% Calculating pole available time (based on image correlation)
ref_img = mean(pole_img_ts(:,:,floor(nof/2 - nof/10):floor(nof/2+nof/10)),3);
[pole_edge_ref, edge_thresh] = edge(imgaussfilt(ref_img,3),'Prewitt','nothinning');

im_corr = zeros(nof,1);
for i = 1 : nof
    pole_edge_ts = edge(imgaussfilt(pole_img_ts(:,:,i),3),'Prewitt', 'nothinning', edge_thresh);
    im_corr(i) = sum(sum(xcorr2(single(pole_edge_ts), single(pole_edge_ref))));
end                
corr_upper = im_corr(im_corr > (max(im_corr)+min(im_corr))/2);

threshold = prctile(corr_upper(floor(length(corr_upper)/5) : length(corr_upper) - floor(length(corr_upper)/5)),20);

pole_available_frames = find(im_corr >= threshold,1,'first') : find(im_corr >= threshold,1,'last');
varargout{3} = pole_available_frames;
end