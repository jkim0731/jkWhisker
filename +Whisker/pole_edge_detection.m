function [pole_edge, varargout] = pole_edge_detection(video_fn)
% every value is in MATLAB convention (
    if isnumeric(video_fn)
        v = VideoReader([num2str(video_fn),'.mp4']);
    else
        v = VideoReader([video_fn,'.mp4']);
    end
    vv = read(v,1);
    width = size(vv,2);
    height = size(vv,1);
    vavg = zeros(height,width);
    for i = 1:v.NumberOfFrames
        vtemp = read(v,i);    
        vtemp = double(vtemp(:,:,1));
        vavg = vavg + vtemp/v.NumberOfFrames;
    end
    vavg = mat2gray(vavg);
%%    
    pole_edge_pic = edge(vavg,'Roberts');
    
    front_pole_edge = pole_edge_pic; front_pole_edge(floor(height*0.7):end,:) = 0; front_pole_edge(:, floor(width*0.5):end) = 0; front_pole_edge(:,1:10) = 0;    
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
    
    if length(result_j) > width*0.2 % linear fit of the front pole
        p = polyfit(result_j(1:floor(width*0.2)),result_i(1:floor(width*0.2)),1);
    else
        p = polyfit(result_j,result_i,1);
    end
    
    q = linspace(1,width*0.5);
%     origin_j = 1;
%     origin_i = polyval(p,origin_j);
%     end_j = width*0.5;
%     end_i = polyval(p,end_j);
   
    pole_edge{2} = [result_i, result_j];
    varargout{1}{2} = [polyval(p,q); q];
    
%%
    top_pole_edge = pole_edge_pic; top_pole_edge(floor(height*0.7):end,:) = 0; top_pole_edge(:, 1:floor(width*0.5)) = 0;
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
    
    if length(result_j) > width*0.2 % linear fit of the front pole
        p = polyfit(result_j(1:floor(width*0.2)),result_i(1:floor(width*0.2)),1);
    else
        p = polyfit(result_j,result_i,1);
    end
    
    q = linspace(width,width*0.5);
%     origin_j = width;
%     origin_i = polyval(p,origin_j);
%     end_j = width*0.5;
%     end_i = polyval(p,end_j);

    pole_edge{1} = [result_i, result_j];
    varargout{1}{1} = [polyval(p,q); q];
    
    varargout{2} = vavg; 
    
    %% Calculating pole available time (based on top view pole movement)
    
    
end