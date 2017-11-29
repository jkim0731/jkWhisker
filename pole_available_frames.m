base_d = 'D:\TPM\JK\tracked\';
mice = {'JK025', 'JK027', 'JK030'};
sessions = 0:22;
width_factor = 0.2;
height_factor = 0.7;

for mind = 2 : length(mice)
% for mind = 2
    if mind == 2
        sessions = 15:22;
    else
        sessions = 0:22;
    end
    for sind = 1 : length(sessions)
        mouse = mice{mind};
        session = sprintf('S%02d',sessions(sind));

        cd([base_d, mouse, session])

        vlist = dir('*.mp4');
        firsts = zeros(length(vlist),1);
        lasts = zeros(length(vlist),1);
        
        parfor vind = 1 : length(vlist)       
            tic
            [~,fn,~] = fileparts(vlist(vind).name);
            if exist([fn,'_WT.mat'],'file') % Let's run this only on those with WT files (otherwise there has been an error in the video file, e.g., interruption)
                v = VideoReader(vlist(vind).name);

                nof = fix(v.FrameRate*v.Duration); % number of frames
                vv = zeros(v.height,v.width,nof,'uint8');
                vedge = zeros(v.height,v.width,nof,'uint8');
                pole_shadow = zeros(nof,1);
                for i = 1 : nof
                    temp = readFrame(v);
                    vv(:,:,i) = temp(:,:,1);    
                    pole_shadow(i) = mean(mean(temp(1:round(v.height*height_factor),1:round(v.width*width_factor),1)));    
                end
                % figure, plot(1:length(pole_shadow),pole_shadow), hold on, plot(1:length(pole_shadow), ones(1,length(pole_shadow))*(mean(pole_shadow)), 'r-');

                pole_candidate = find(pole_shadow < mean(pole_shadow));
                bg_frames = find(pole_shadow > min([pole_shadow(1:100); pole_shadow(end-100:end)]));
                bg = imgaussfilt(mean(vv(:,:,bg_frames),3),2);

                % %%
                pole_im = zeros(round(v.height*height_factor),round(v.width*width_factor),length(pole_candidate));
                pole_bw = zeros(round(v.height*height_factor),round(v.width*width_factor),length(pole_candidate),'logical');
                edge_front = zeros(length(pole_candidate),1);
                for i = 1 : length(pole_candidate)
                    pole_im(:,:,i) = bg(1:round(v.height*height_factor),1:round(v.width*width_factor)) - imgaussfilt(double(vv(1:round(v.height*height_factor),1:round(v.width*width_factor),pole_candidate(i))),2);
                    pole_bw(:,:,i) = pole_im(:,:,i) > 3*std(std(pole_im(:,:,i)));
                    a = regionprops(pole_bw(:,:,i),'Extrema','Area');
                    maxArea = a(1).Area;
                    maxInd = 1;
                    if length(a) > 1
                        for j = 2 : length(a)
                            if a(j).Area > maxArea
                                maxArea = a(j).Area;
                                maxInd = j;
                            end
                        end
                    end    
                    edge_front(i) = a(maxInd).Extrema(5,2);
                end            

                [c, p] = histcounts(edge_front,min(edge_front):0.25:max(edge_front));
                [~,maxInd] = max(c);
                maxPos = p(maxInd-1);

                firsts(vind) = find(edge_front >= maxPos,1,'first') + pole_candidate(1) - 1;
                lasts(vind) = nof - find(edge_front >= maxPos,1,'last') - pole_candidate(1) + 1;
                %%
            %     figure, plot(1:length(edge_front),edge_front), hold on,
            %     plot(first,edge_front(first),'r.','MarkerSize',20)
            %     plot(last,edge_front(last),'r.','MarkerSize',20)
                %%
            %     implay(pole_bw)
                fprintf('%s %s %d/%d \n', mouse, session, vind, length(vlist));            
            else
                firsts(vind) = NaN; lasts(vind) = NaN;
                fprintf('Whisker tracking error: %s %s %d/%d \n', mouse, session, vind, length(vlist));
            end
        end
        save('firsts_n_lasts.mat','firsts','lasts')
        clear vlist v
    end
    
end
