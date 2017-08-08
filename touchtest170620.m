close all
clear
base_angle = 21;
behavior_base_dir = 'Z:\Data\2p\soloData\';
whisker_base_dir = 'Z:\Data\Video\JK\';
trial_types = {'rc', 'rf', 'lc', 'lf'};
mice = {'AH0648','AH0650','AH0651','AH0652','AH0653'};
% mice = {'AH0650'};
sessionNum = {[1:4,6:15],[1,2,4:6,8:10],[2,4:6,8:19],[6,8:13],[2,4:8,13:19]};
% sessionNum = {[21:24]};

theta_T = cell(length(mice),1);
theta_F = cell(length(mice),1);
for i = 1 : length(mice)    
    theta_T{i} = cell(length(sessionNum{i}),1);
    theta_F{i} = cell(length(sessionNum{i}),1);
    for j = 1 : length(sessionNum{i})
        theta_T{i}{j} = cell(length(trial_types),1); 
        theta_F{i}{j} = cell(length(trial_types),1); 
        for k = 1 : length(trial_types)
            theta_T{i}{j}{k} = cell(2,1); % {1} for free-whisking, {2} for touch
            theta_F{i}{j}{k} = cell(2,1);
            for l = 1 : 2
                theta_T{i}{j}{k}{l} = [];
                theta_F{i}{j}{k}{l} = [];
            end
        end
    end
end

for mi = 1 : length(mice)
    mouseName = mice{mi};
    for si = 1:length(sessionNum{mi})
        sessionName = sprintf('S%02d',sessionNum{mi}(si));
        whisker_d = [whisker_base_dir mouseName sessionName '\'];
        if ~exist('b','var') || ~iscell(b) || ~isfield(b{1},'mouseName') || ~strcmp(b{1}.mouseName,mouseName)
            load([behavior_base_dir mouseName filesep 'behavior.mat']) % load b
        end
        b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
        b_session = b{b_ind};

        tt_ind = cell(1,length(trial_types));

        for ti = 1 : length(trial_types)    
            tt_ind{ti} = find(cellfun(@(x) strcmp(x.trialType,trial_types{ti}),b_session.trials));
            temp_files = cell(length(tt_ind{ti}),1);
            for j = 1 : length(tt_ind{ti})
                temp_files{j} = num2str(tt_ind{ti}(j));
            end
            wl = Whisker.WhiskerTrialLiteArray_2pad(whisker_d,'include_files',temp_files);            
            temp_Tt = cellfun(@(x) x.thetaAtBase{1}(x.th_touch_frames),wl.trials, 'UniformOutput', false); % results in cell
            temp_Ft = cellfun(@(x) x.thetaAtBase{2}(x.th_touch_frames),wl.trials, 'UniformOutput', false);
            temp_Tw = cellfun(@(x) x.thetaAtBase{1}(x.time{1}<1.6),wl.trials, 'UniformOutput', false);
            temp_Fw = cellfun(@(x) x.thetaAtBase{2}(x.time{1}<1.6),wl.trials, 'UniformOutput', false);
            for i = 1 : length(temp_Tt)
                theta_T{mi}{si}{ti}{1} = [theta_T{mi}{si}{ti}{1}, temp_Tw{i}];
                theta_F{mi}{si}{ti}{1} = [theta_F{mi}{si}{ti}{1}, temp_Fw{i}];
                theta_T{mi}{si}{ti}{2} = [theta_T{mi}{si}{ti}{2}, temp_Tt{i}];
                theta_F{mi}{si}{ti}{2} = [theta_F{mi}{si}{ti}{2}, temp_Ft{i}];
            end            
        end
    end
end

% %%
for i = 1 : length(mice)    
    for j = 1 : length(sessionNum{i})
        for k = 1 : length(trial_types)
            for l = 1 : 2
                theta_T{i}{j}{k}{l} = theta_T{i}{j}{k}{l}+21;
            end
        end
    end
end

save('tethas.mat', 'theta_*')
%%
thetaTrange = -40:20;
mouse = 2;
session = 8;
thetaThist = cell(length(trial_types),2);
for wi = 1 : 2
    thetaThist{1,wi} = histcounts(theta_T{mouse}{session}{1}{wi},[thetaTrange,thetaTrange(end)+1]);
    thetaThist{2,wi} = histcounts(theta_T{mouse}{session}{2}{wi},[thetaTrange,thetaTrange(end)+1]);
    thetaThist{3,wi} = histcounts(theta_T{mouse}{session}{3}{wi},[thetaTrange,thetaTrange(end)+1]);
    thetaThist{4,wi} = histcounts(theta_T{mouse}{session}{4}{wi},[thetaTrange,thetaTrange(end)+1]);
end
t = zeros(size(thetaThist{1}));
for i = 1 : 4
    for j = 1 : 2
        t = t + thetaThist{i,j};
    end
end
figure, set(gcf,'Units','centimeters','Position',[8 8 11.7 9.5]);
plot(thetaTrange,thetaThist{1,2},'Color','red','LineWidth',2), hold on
plot(thetaTrange,thetaThist{2,2},'Color','magenta','LineWidth',2)
plot(thetaTrange,thetaThist{3,2},'Color','blue','LineWidth',2)
plot(thetaTrange,thetaThist{4,2},'Color','cyan','LineWidth',2)
plot(thetaTrange,thetaThist{1,1}+thetaThist{2,1}+thetaThist{3,1}+thetaThist{4,1},'Color','green','LineWidth',2)
plot(thetaTrange,thetaThist{1,2}+thetaThist{2,2}+thetaThist{3,2}+thetaThist{4,2},'Color',[0.5 0.5 0.5],'LineWidth',2)
plot(thetaTrange,t,'k-','LineWidth',2)
xlabel('Theta'), title({[mice{mouse},sprintf('S%02d',sessionNum{mouse}(session)), ' Top-view']}), 
legend('Close down','Far down','Close up', 'Far up','Free Whisking', 'All trial types')

%%
thetaFrange = -10:30;
mouse = 2;
session = 8;
thetaFhist = cell(length(trial_types),2);
for wi = 1 : 2
    thetaFhist{1,wi} = histcounts(theta_F{mouse}{session}{1}{wi},[thetaFrange,thetaFrange(end)+1]);
    thetaFhist{2,wi} = histcounts(theta_F{mouse}{session}{2}{wi},[thetaFrange,thetaFrange(end)+1]);
    thetaFhist{3,wi} = histcounts(theta_F{mouse}{session}{3}{wi},[thetaFrange,thetaFrange(end)+1]);
    thetaFhist{4,wi} = histcounts(theta_F{mouse}{session}{4}{wi},[thetaFrange,thetaFrange(end)+1]);
end
t = zeros(size(thetaFhist{1}));
for i = 1 : 4
    for j = 1 : 2
        t = t + thetaFhist{i,j};
    end
end
figure, set(gcf,'Units','centimeters','Position',[8 8 11.7 9.5]);
plot(thetaFrange,thetaFhist{1,2},'Color','red','LineWidth',2), hold on
plot(thetaFrange,thetaFhist{2,2},'Color','magenta','LineWidth',2)
plot(thetaFrange,thetaFhist{3,2},'Color','blue','LineWidth',2)
plot(thetaFrange,thetaFhist{4,2},'Color','cyan','LineWidth',2)
plot(thetaFrange,thetaFhist{1,1}+thetaFhist{2,1}+thetaFhist{3,1}+thetaFhist{4,1},'Color','green','LineWidth',2)
plot(thetaFrange,thetaFhist{1,2}+thetaFhist{2,2}+thetaFhist{3,2}+thetaFhist{4,2},'Color',[0.5 0.5 0.5],'LineWidth',2)
plot(thetaFrange,t,'k-','LineWidth',2)
xlabel('Theta'), title({[mice{mouse},sprintf('S%02d',sessionNum{mouse}(session)), ' Front-view']}) 
%%
thetaTrange = -40:10;
mouse = 1;
thetaThist = cell(length(sessionNum{mouse}),1);
for si = 1 : length(sessionNum{mouse})
    thetaThist{si} = zeros(1,length(thetaTrange));
    for i = 1 : 4
        for j = 1 : 2
            thetaThist{si} = thetaThist{si} + histcounts(theta_T{mouse}{si}{i}{j},[thetaTrange,thetaTrange(end)+1]);
        end
    end
end
figure, set(gcf,'Units','centimeters','Position',[8 8 11.7 9.5]);
alls = zeros(1,length(thetaTrange));
for si = 1 : length(sessionNum{mouse})
    cval = 0.5+0.5/length(sessionNum{mouse})*si;
    plot(thetaTrange,thetaThist{si}, 'Color', [cval cval cval], 'LineWidth', 2), hold on
    alls = alls + thetaThist{si};
end
plot(thetaTrange,alls/length(sessionNum{mouse}),'r-','LineWidth',2)
xlabel('Theta'), title('Telecentric lens top-view all sessions') 
%%
thetaFrange = -10:30;
mouse = 2;
thetaFhist = cell(length(sessionNum{mouse}),1);
for si = 1 : length(sessionNum{mouse})
    thetaFhist{si} = zeros(1,length(thetaFrange));
    for i = 1 : 4
        for j = 1 : 2
            thetaFhist{si} = thetaFhist{si} + histcounts(theta_F{mouse}{si}{i}{j},[thetaFrange,thetaFrange(end)+1]);
        end
    end
end
figure, set(gcf,'Units','centimeters','Position',[8 8 11.7 9.5]);
alls = zeros(1,length(thetaFrange));
for si = 1 : length(sessionNum{mouse})
    cval = 0.5+0.5/length(sessionNum{mouse})*si;
    plot(thetaFrange,thetaFhist{si}, 'Color', [cval cval cval], 'LineWidth', 2), hold on
    alls = alls + thetaFhist{si};
end
plot(thetaFrange,alls/length(sessionNum{mouse}),'r-','LineWidth',2)
xlabel('Theta'), title('Telecentric lens Front-view all sessions') 
%%
thetaTrange = -30:20;
thetaThist = cell(length(mice),1);
for mi = 1 : length(mice)
    thetaThist{mi} = zeros(1,length(thetaTrange));
    for si = 1 : length(sessionNum{mi})        
        for i = 1 : 4
            for j = 1 : 2
                thetaThist{mi} = thetaThist{mi} + histcounts(theta_T{mi}{si}{i}{j},[thetaTrange,thetaTrange(end)+1]);
            end
        end
    end
end
figure, set(gcf,'Units','centimeters','Position',[8 8 11.7 9.5]);
alls = zeros(1,length(thetaTrange));
for mi = 1 : length(mice)
    cval = 0.5+0.5/length(mice)*mi;
    plot(thetaTrange,thetaThist{mi}, 'Color', [cval cval cval], 'LineWidth', 2), hold on
    alls = alls + thetaThist{mi};
end
plot(thetaTrange,alls/length(mice),'r-','LineWidth',2)
xlabel('Theta'), title('Top-view all mice') 
%%
thetaFrange = 0:30;
thetaFhist = cell(length(mice),1);
for mi = 1 : length(mice)
    thetaFhist{mi} = zeros(1,length(thetaFrange));
    for si = 1 : length(sessionNum{mi})        
        for i = 1 : 4
            for j = 1 : 2
                thetaFhist{mi} = thetaFhist{mi} + histcounts(theta_F{mi}{si}{i}{j},[thetaFrange,thetaFrange(end)+1]);
            end
        end
    end
end
figure, set(gcf,'Units','centimeters','Position',[8 8 11.7 9.5]);
alls = zeros(1,length(thetaFrange));
for mi = 1 : length(mice)
    cval = 0.5+0.5/length(mice)*mi;
    plot(thetaFrange,thetaFhist{mi}, 'Color', [cval cval cval], 'LineWidth', 2), hold on
    alls = alls + thetaFhist{mi};
end
plot(thetaFrange,alls/length(mice),'r-','LineWidth',2)
xlabel('Theta'), title('Front-view all mice') 
%%
% mouseName = 'AH0653';
sessions{1} = {[4,6,7], [13:15]}; % 0648
sessions{2} = {[4:6], [10]}; % 0650
sessions{3} = {[4:6], [6,8,10]}; % 0651
sessions{4} = {[6,8], [10,11,13]}; % 0652
sessions{5} = {[4:6], [16:18]}; % 0653

theta_T = cell(length(mice),1);
kappa_T = cell(length(mice),1);
theta_F = cell(length(mice),1);
kappa_F = cell(length(mice),1);
for i = 1 : length(mice)
    theta_T{i} = cell(2,1); % before and after 
    kappa_T{i} = cell(2,1);
    theta_F{i} = cell(2,1);
    kappa_F{i} = cell(2,1);    
    for j = 1 : 2
        theta_T{i}{j} = cell(length(trial_types),1); 
        kappa_T{i}{j} = cell(length(trial_types),1);
        theta_F{i}{j} = cell(length(trial_types),1);
        kappa_F{i}{j} = cell(length(trial_types),1);
        for k = 1 : length(trial_types)
            theta_T{i}{j}{k} = cell(2,1); % pre-decision at the second cell
            kappa_T{i}{j}{k} = cell(2,1);
            theta_F{i}{j}{k} = cell(2,1);
            kappa_F{i}{j}{k} = cell(2,1);
%             for l = 1 : 2
%                 theta_T{i}{j}{k}{l} = {};
%                 kappa_T{i}{j}{k}{l} = {};
%                 theta_F{i}{j}{k}{l} = {};
%                 kappa_F{i}{j}{k}{l} = {};
%             end
        end
    end
end

%%
% 
% touch_chunks = cell(length(sessions),1);
% touch_chunks_pre = cell(length(sessions),1);
% for i = 1 : length(touch_chunks)
%     touch_chunks{i} = [];
%     touch_chunks_pre{i} = [];
% end
% % for i = 1 : length(sessions)
% %     touch_chunks{i} = cell(length(sessions{i}),1);
% %     for j = 1 : length(sessions{i})
% %         touch_chunks{i}{j} = cell(length(trial_types),1);
% %         for k = 1 : length(trial_types)
% %             touch_chunks{i}{j}{k} = [];
% %         end
% %     end
% % end
for mi = 1 : length(mice)
    mouseName = mice{mi};
for sind = 1 : length(sessions{mi})
    for s = 1 : length(sessions{mi}{sind})
        sessionName = sprintf('S%02d', sessions{mi}{sind}(s));
        whisker_d = [whisker_base_dir mouseName sessionName '\'];
        if ~exist('b','var') || ~iscell(b) || ~isfield(b{1},'mouseName') || ~strcmp(b{1}.mouseName,mouseName)
            load([behavior_base_dir mouseName filesep 'behavior.mat']) % load b
        end
        b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
        b_session = b{b_ind};

        tt_ind = cell(1,length(trial_types));
        wl_array = cell(1,length(trial_types));

        for trial_type_num = 1 : length(trial_types)    
            tt_ind{trial_type_num} = find(cellfun(@(x) strcmp(x.trialType,trial_types{trial_type_num}),b_session.trials));
            temp_files = cell(length(tt_ind{trial_type_num}),1);
            for j = 1 : length(tt_ind{trial_type_num})
                temp_files{j} = num2str(tt_ind{trial_type_num}(j));
            end
            wl = Whisker.WhiskerTrialLiteArray_2pad(whisker_d,'include_files',temp_files);
            wl_array{trial_type_num} = wl;
        end
        for wl_array_ind = 1 : length(trial_types)
            for wl_ind = 1 : length(wl_array{wl_array_ind}.trials)
                wt = wl_array{wl_array_ind}.trials{wl_ind};
                trial_temp_ind = find(cellfun(@(x) x.trialNum == wt.trialNum,b_session.trials));
                trial_temp = b_session.trials{trial_temp_ind};
                
                temp_touch = wt.th_touch_frames;
                if ~isempty(temp_touch)
                    touch_diff_inds = [0;find(diff(temp_touch) - 1);length(temp_touch)];
                    for cind = 1 : length(touch_diff_inds)-1
                        if length(wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * ...
                                length(wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) ...
                                > 0 % non of these are empty
                            theta_T{mi}{sind}{wl_array_ind}{1}{end+1} = wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) + base_angle;
                            kappa_T{mi}{sind}{wl_array_ind}{1}{end+1} = wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                            theta_F{mi}{sind}{wl_array_ind}{1}{end+1} = wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                            kappa_F{mi}{sind}{wl_array_ind}{1}{end+1} = wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                        end
                    end
                end                
                
                if ~isempty(trial_temp.beamBreakTimes)
                    first_lick_time = min(trial_temp.beamBreakTimes);                    
                else
                    first_lick_frame = 0;
                end                
                before_frames_max = find(wt.time{1} < first_lick_time,1,'last');
                before_frames = wt.th_touch_frames(wt.th_touch_frames < before_frames_max);                
                if ~isempty(before_frames)
                    touch_diff_inds = [0; find(diff(before_frames) - 1)];
                    for cind = 1 : length(touch_diff_inds)-1
                        if length(wt.thetaAtBase{1}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{1}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * ...
                                length(wt.thetaAtBase{2}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{2}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) ...
                                > 0 % non of these are empty
                            theta_T{mi}{sind}{wl_array_ind}{2}{end+1} = wt.thetaAtBase{1}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) + base_angle;
                            kappa_T{mi}{sind}{wl_array_ind}{2}{end+1} = wt.deltaKappa{1}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                            theta_F{mi}{sind}{wl_array_ind}{2}{end+1} = wt.thetaAtBase{2}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                            kappa_F{mi}{sind}{wl_array_ind}{2}{end+1} = wt.deltaKappa{2}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                        end
                    end
                end                
            end
        end
    end
end
end
%%
close all
maxhist = 601;
for i = 1 : length(mice)
    figure, title({mice{i}; 'Duration of touches'}), hold on,
    bl = []; al = []; blpd = []; alpd = [];
    for j = 1 : length(trial_types)        
        bl = [bl, cellfun(@(x) length(x), theta_T{i}{1}{j}{1})];
        al = [al, cellfun(@(x) length(x), theta_T{i}{2}{j}{1})];
        blpd = [blpd, cellfun(@(x) length(x), theta_T{i}{1}{j}{2})];
        alpd = [alpd, cellfun(@(x) length(x), theta_T{i}{2}{j}{2})];
    end
    blhist = histcounts(bl,30:maxhist);
    alhist = histcounts(al,30:maxhist);
    blpdhist = histcounts(blpd,30:maxhist);
    alpdhist = histcounts(alpd,30:maxhist);
    plot(30:maxhist-1,blhist,'m--','LineWidth',2), hold on,
    plot(30:maxhist-1,alhist,'c--','LineWidth',2)
    plot(30:maxhist-1,blpdhist,'r-','LineWidth',2)
    plot(30:maxhist-1,alpdhist,'b-','LineWidth',2),
    ylim([0 100]), xlabel('Duration (frame; 0.0032 sec)')    
    if i > 1 && i < 5
        legend('Before learning','After learning','Before learning pre-decision','After learning pre-decision');
    else        
        legend('Early phase','Late phase','Early phase pre-decision','Late phase pre-decision');
    end
end

%%
close all
minhist = -30; maxhist = 30;
for i = 1 : length(mice)
    figure,     
    for j = 1 : 2                
        for k = 1 : 2
            cd = cell2mat(theta_T{i}{j}{1}{k});
            fd = cell2mat(theta_T{i}{j}{2}{k});
            cu = cell2mat(theta_T{i}{j}{3}{k});
            fu = cell2mat(theta_T{i}{j}{4}{k});
            cdhist = histcounts(cd,minhist:maxhist+1);
            fdhist = histcounts(fd,minhist:maxhist+1);
            cuhist = histcounts(cu,minhist:maxhist+1);
            fuhist = histcounts(fu,minhist:maxhist+1);        
            subplot(2,2,(j-1)*2+k), plot(minhist:maxhist,cdhist,'r-','LineWidth',2), hold on,
            plot(minhist:maxhist,fdhist,'m-','LineWidth',2)
            plot(minhist:maxhist,cuhist,'b-','LineWidth',2)
            plot(minhist:maxhist,fuhist,'c-','LineWidth',2)
            switch (j-1)*2 + k
                case 1
                    titletext = 'Early phase';
                case 2
                    titletext = 'Early phase / pre-decision';
                case 3
                    titletext = 'Late phase';
                case 4
                    titletext = 'Late phase / pre-decision';
            end
            title({[mice{i} 'Touch angle span'];titletext}), xlim([minhist maxhist])
        end            
    end
end


%%

%%
close all
minhist = 0; maxhist = 30;
for i = 1 : length(mice)
    figure,     
    for j = 1 : 2                
        for k = 1 : 2
            cd = cell2mat(theta_F{i}{j}{1}{k});
            fd = cell2mat(theta_F{i}{j}{2}{k});
            cu = cell2mat(theta_F{i}{j}{3}{k});
            fu = cell2mat(theta_F{i}{j}{4}{k});
            cdhist = histcounts(cd,minhist:maxhist+1);
            fdhist = histcounts(fd,minhist:maxhist+1);
            cuhist = histcounts(cu,minhist:maxhist+1);
            fuhist = histcounts(fu,minhist:maxhist+1);        
            subplot(2,2,(j-1)*2+k), plot(minhist:maxhist,cdhist,'r-','LineWidth',2), hold on,
            plot(minhist:maxhist,fdhist,'m-','LineWidth',2)
            plot(minhist:maxhist,cuhist,'b-','LineWidth',2)
            plot(minhist:maxhist,fuhist,'c-','LineWidth',2)
            switch (j-1)*2 + k
                case 1
                    titletext = 'Early phase';
                case 2
                    titletext = 'Early phase / pre-decision';
                case 3
                    titletext = 'Late phase';
                case 4
                    titletext = 'Late phase / pre-decision';
            end
            title({[mice{i} 'Touch angle span'];titletext}), xlim([minhist maxhist])
        end            
    end
end

%%
% mouseName = 'AH0653';
sessions{1} = {[4,6,7], [13:15]}; % 0648
sessions{2} = {[4:6], [10]}; % 0650
sessions{3} = {[4:6], [6,8,10]}; % 0651
sessions{4} = {[6,8], [10,11,13]}; % 0652
sessions{5} = {[4:6], [16:18]}; % 0653

theta_fw_T = cell(length(mice),1);
kappa_fw_T = cell(length(mice),1);
theta_fw_F = cell(length(mice),1);
kappa_fw_F = cell(length(mice),1);
for i = 1 : length(mice)
    theta_fw_T{i} = cell(2,1); % before and after 
    kappa_fw_T{i} = cell(2,1);
    theta_fw_F{i} = cell(2,1);
    kappa_fw_F{i} = cell(2,1);    
    for j = 1 : 2
        theta_fw_T{i}{j} = cell(length(trial_types),1); 
        kappa_fw_T{i}{j} = cell(length(trial_types),1);
        theta_fw_F{i}{j} = cell(length(trial_types),1);
        kappa_fw_F{i}{j} = cell(length(trial_types),1);
        for k = 1 : length(trial_types)
            theta_fw_T{i}{j}{k} = cell(2,1); % pre-decision at the second cell
            kappa_fw_T{i}{j}{k} = cell(2,1);
            theta_fw_F{i}{j}{k} = cell(2,1);
            kappa_fw_F{i}{j}{k} = cell(2,1);
%             for l = 1 : 2
%                 theta_T{i}{j}{k}{l} = {};
%                 kappa_T{i}{j}{k}{l} = {};
%                 theta_F{i}{j}{k}{l} = {};
%                 kappa_F{i}{j}{k}{l} = {};
%             end
        end
    end
end

%%
% 
% touch_chunks = cell(length(sessions),1);
% touch_chunks_pre = cell(length(sessions),1);
% for i = 1 : length(touch_chunks)
%     touch_chunks{i} = [];
%     touch_chunks_pre{i} = [];
% end
% % for i = 1 : length(sessions)
% %     touch_chunks{i} = cell(length(sessions{i}),1);
% %     for j = 1 : length(sessions{i})
% %         touch_chunks{i}{j} = cell(length(trial_types),1);
% %         for k = 1 : length(trial_types)
% %             touch_chunks{i}{j}{k} = [];
% %         end
% %     end
% % end
for mi = 1 : length(mice)
    mouseName = mice{mi};
    for sind = 1 : length(sessions{mi})
        for s = 1 : length(sessions{mi}{sind})
            sessionName = sprintf('S%02d', sessions{mi}{sind}(s));
            whisker_d = [whisker_base_dir mouseName sessionName '\'];
            if ~exist('b','var') || ~iscell(b) || ~isfield(b{1},'mouseName') || ~strcmp(b{1}.mouseName,mouseName)
                load([behavior_base_dir mouseName filesep 'behavior.mat']) % load b
            end
            b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
            b_session = b{b_ind};

            tt_ind = cell(1,length(trial_types));
            wl_array = cell(1,length(trial_types));

            for trial_type_num = 1 : length(trial_types)    
                tt_ind{trial_type_num} = find(cellfun(@(x) strcmp(x.trialType,trial_types{trial_type_num}),b_session.trials));
                temp_files = cell(length(tt_ind{trial_type_num}),1);
                for j = 1 : length(tt_ind{trial_type_num})
                    temp_files{j} = num2str(tt_ind{trial_type_num}(j));
                end
                wl = Whisker.WhiskerTrialLiteArray_2pad(whisker_d,'include_files',temp_files);
                wl_array{trial_type_num} = wl;
            end
            for wl_array_ind = 1 : length(trial_types)
                for wl_ind = 1 : length(wl_array{wl_array_ind}.trials)
                    wt = wl_array{wl_array_ind}.trials{wl_ind};
                    trial_temp_ind = find(cellfun(@(x) x.trialNum == wt.trialNum,b_session.trials));
                    trial_temp = b_session.trials{trial_temp_ind};

                    temp_touch = wt.th_touch_frames;
                    if ~isempty(temp_touch)
                        touch_diff_inds = [0;find(diff(temp_touch) - 1);length(temp_touch)];
                        for cind = 1 : length(touch_diff_inds)-1
                            if length(wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * ...
                                    length(wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) ...
                                    > 0 % non of these are empty
                                theta_fw_T{mi}{sind}{wl_array_ind}{1}{end+1} = wt.thetaAtBase{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) + base_angle;
                                kappa_fw_T{mi}{sind}{wl_array_ind}{1}{end+1} = wt.deltaKappa{1}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                                theta_fw_F{mi}{sind}{wl_array_ind}{1}{end+1} = wt.thetaAtBase{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                                kappa_fw_F{mi}{sind}{wl_array_ind}{1}{end+1} = wt.deltaKappa{2}(temp_touch(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                            end
                        end
                    end                

                    if ~isempty(trial_temp.beamBreakTimes)
                        first_lick_time = min(trial_temp.beamBreakTimes);                    
                    else
                        first_lick_frame = 0;
                    end                
                    before_frames_max = find(wt.time{1} < first_lick_time,1,'last');
                    before_frames = wt.th_touch_frames(wt.th_touch_frames < before_frames_max);                
                    if ~isempty(before_frames)
                        touch_diff_inds = [0; find(diff(before_frames) - 1)];
                        for cind = 1 : length(touch_diff_inds)-1
                            if length(wt.thetaAtBase{1}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{1}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * ...
                                    length(wt.thetaAtBase{2}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) * length(wt.deltaKappa{2}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)))) ...
                                    > 0 % non of these are empty
                                theta_fw_T{mi}{sind}{wl_array_ind}{2}{end+1} = wt.thetaAtBase{1}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1))) + base_angle;
                                kappa_fw_T{mi}{sind}{wl_array_ind}{2}{end+1} = wt.deltaKappa{1}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                                theta_fw_F{mi}{sind}{wl_array_ind}{2}{end+1} = wt.thetaAtBase{2}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                                kappa_fw_F{mi}{sind}{wl_array_ind}{2}{end+1} = wt.deltaKappa{2}(before_frames(touch_diff_inds(cind)+1:touch_diff_inds(cind+1)));
                            end
                        end
                    end                
                end
            end
        end
    end
end

save('thetakappa.mat','theta_*','kappa_*')

%% 
% xlim_val = 40;
% h1 = histcounts(touch_chunks{1},1:xlim_val+1); 
% h2 = histcounts(touch_chunks{2},1:xlim_val+1); 
% h3 = histcounts(touch_chunks_pre{1},1:xlim_val+1);
% h4 = histcounts(touch_chunks_pre{2},1:xlim_val+1);
% figure, 
% plot(1:xlim_val,h1,'c--', 'LineWidth', 2), hold on,
% plot(1:xlim_val,h2,'m--', 'LineWidth', 2)
% plot(1:xlim_val,h3,'b-', 'LineWidth', 2)
% plot(1:xlim_val,h4,'r-', 'LineWidth', 2)
% title({mouseName;'Duration of touches'}), legend('Early phase','Late phase','Early phase pre-decision','Late phase pre-decision'), xlabel('Duration (frames; 0.0032sec)')
% xlim([01 80]), title({'AH0650';'Duration of touches'}), legend('Before Learning','After Learning'), xlabel('Duration (frames; 0.0032sec)')
% %% pre-decision vs after-result
% sgind = 2; % only in late phase
% spg = cell(2,1); % before licking vs after reward
% for i = 1 : 2
%     spg{i} = cell(4,1); % figures 1~4: rc rf lc lf
%     for j = 1 : 4
%         spg{i}{j} = cell(5,1); % subplot 1~5
%     end
% end
% for sind = 1 : length(sessions{sgind})
%     sessionName = sprintf('S%02d',sessions{sgind}(sind));
%     whisker_d = [whisker_base_dir mouseName sessionName '\'];
%     if ~exist('b','var') || ~iscell(b) || ~isfield(b{1},'mouseName') || ~strcmp(b{1}.mouseName,mouseName)
%         load([behavior_base_dir mouseName filesep 'behavior.mat']) % load b
%     end
%     b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
%     b_session = b{b_ind};
%     
%     tt_ind = cell(1,length(trial_types));
%     wl_array = cell(1,length(trial_types));
% 
%     for trial_type_num = 1 : length(trial_types)    
%         tt_ind{trial_type_num} = find(cellfun(@(x) strcmp(x.trialType,trial_types{trial_type_num}),b_session.trials));
%         temp_files = cell(length(tt_ind{trial_type_num}),1);
%         for j = 1 : length(tt_ind{trial_type_num})
%             temp_files{j} = num2str(tt_ind{trial_type_num}(j));
%         end
%         wl = Whisker.WhiskerTrialLiteArray_2pad(whisker_d,'include_files',temp_files);
%         wl_array{trial_type_num} = wl;
%     end 
%     
%     for wl_array_ind = 1 : length(trial_types)
%         for wl_ind = 1 : length(wl_array{wl_array_ind}.trials)
%             wt = wl_array{wl_array_ind}.trials{wl_ind};
%             trial_temp_ind = find(cellfun(@(x) x.trialNum == wt.trialNum,b_session.trials));
%             trial_temp = b_session.trials{trial_temp_ind};
%             
%             if ~isempty(trial_temp.beamBreakTimes)
%                 first_lick_time = min(trial_temp.beamBreakTimes);
%                 first_lick_frame = first_lick_time / wt.framePeriodInSec;
%             else
%                 first_lick_frame = 0;
%             end
%             
%             if trial_temp.trialCorrect == 1
%                 result_time = trial_temp.drinkingTime(1);
%             elseif trial_temp.trialCorrect == 0
%                 result_time = trial_temp.timeoutPeriodTimes{1}(1);
%             else
%                 result_time = max(wt.time{1});
%             end
%             result_frame = result_time / wt.framePeriodInSec;
%             
%             before_frames = wt.th_touch_frames(wt.th_touch_frames < first_lick_frame);
%             after_frames = wt.th_touch_frames(wt.th_touch_frames > result_frame);
%             
%             spg{1}{wl_array_ind}{1} = [spg{1}{wl_array_ind}{1}; wt.thetaAtBase{1}(before_frames)'+base_angle, wt.thetaAtBase{2}(before_frames)'];            
%             spg{1}{wl_array_ind}{2} = [spg{1}{wl_array_ind}{2}; wt.deltaKappa{1}(before_frames)', wt.deltaKappa{2}(before_frames)'];
%             spg{1}{wl_array_ind}{3} = [spg{1}{wl_array_ind}{3}; wt.thetaAtBase{1}(before_frames)'+base_angle, wt.deltaKappa{1}(before_frames)'];
%             spg{1}{wl_array_ind}{4} = [spg{1}{wl_array_ind}{4}; wt.thetaAtBase{2}(before_frames)', wt.deltaKappa{2}(before_frames)'];
%             spg{1}{wl_array_ind}{5} = [spg{1}{wl_array_ind}{5}; wt.thetaAtBase{1}(before_frames)'+base_angle, wt.deltaKappa{2}(before_frames)'];        
%             
%             spg{2}{wl_array_ind}{1} = [spg{2}{wl_array_ind}{1}; wt.thetaAtBase{1}(after_frames)'+base_angle, wt.thetaAtBase{2}(after_frames)'];            
%             spg{2}{wl_array_ind}{2} = [spg{2}{wl_array_ind}{2}; wt.deltaKappa{1}(after_frames)', wt.deltaKappa{2}(after_frames)'];
%             spg{2}{wl_array_ind}{3} = [spg{2}{wl_array_ind}{3}; wt.thetaAtBase{1}(after_frames)'+base_angle, wt.deltaKappa{1}(after_frames)'];
%             spg{2}{wl_array_ind}{4} = [spg{2}{wl_array_ind}{4}; wt.thetaAtBase{2}(after_frames)', wt.deltaKappa{2}(after_frames)'];
%             spg{2}{wl_array_ind}{5} = [spg{2}{wl_array_ind}{5}; wt.thetaAtBase{1}(after_frames)'+base_angle, wt.deltaKappa{2}(after_frames)'];     
%         end
%     end
% end
% 
% %% 
% close all
% for i = 1 : 4 % different figure
%     switch i
%         case 1
%             tt_text = 'Close Down';
%         case 2
%             tt_text = 'Far Down';
%         case 3
%             tt_text = 'Close Up';
%         case 4
%             tt_text = 'Far Up';
%     end
%     figure,set(gcf,'Units','centimeters','Position',[i+3 i+3 11.7 9.5]);     
%     for j = 1 : 5 % different subplot
%         subplot(3,2,j), hold on
%         scatter(spg{1}{i}{j}(:,1),spg{1}{i}{j}(:,2),0.5,'k.'); % , 'LineWidth', 0.001
%         scatter(spg{2}{i}{j}(:,1),spg{2}{i}{j}(:,2),0.5, 'b.', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha',0.2); % , 'LineWidth', 0.001, 
% 
%         switch j
%             case 1
%                 title({['                                                       ',tt_text,' touch frames'];'\theta_T vs \theta_F'});
%                 xlim([-30 30]), ylim([-10 30]);
%             case 2
%                 title('\kappa_T vs \kappa_F');
%                 xlim([-0.01 0.01]), ylim([-0.01 0.01]);
%             case 3
%                 title('\theta_T vs \kappa_T');
%                 xlim([-30 30]), ylim([-0.01 0.01]);
%             case 4
%                 title('\theta_F vs \kappa_F');
%                 xlim([-10 30]), ylim([-0.01 0.01]);
%             case 5
%                 title('\theta_T vs \kappa_F');
%                 xlim([-30 30]), ylim([-0.01 0.01]);
%         end
%     end
% end
    
    
% %% before learning vs after learning (or early phase vs late phase for not-learned animals)
% spg = cell(2,1);
% for sgind = 1 : length(sessions) % session group index
%     spg{sgind} = cell(length(sessions{sgind}),1);
%     for sind = 1 : length(sessions{sgind}) % session index
%         sessionName = sprintf('S%02d',sessions{sgind}(sind));
%         behavior_d = [behavior_base_dir mouseName '\'];
%         whisker_d = [whisker_base_dir mouseName sessionName '\'];
% 
%         load([behavior_d 'behavior.mat']) % loading b of the mouse (all the sessions)
% 
%         b_ind = find(cellfun(@(x) strcmp(x.sessionName,sessionName), b));
%         b_session = b{b_ind};
% 
%         load_fn = [mouseName sessionName '_post.mat'];
%         load([whisker_d load_fn]); % loading errorlist
% 
%         if ~isempty(b_ind) % try only the ones with behavior session
%             filelist=dir([whisker_d '*.measurements']);
%             dirTrialNums=zeros(1,size(filelist,1));
% 
%             for i=1:length(filelist)
%                 dirTrialNums(i)=str2double(filelist(i).name(1:end-13)); % extract out the trial number from each measurements file present in directory
%             end
%             dirTrialNums = setdiff(dirTrialNums,errorlist);
%             trialNums = sort(dirTrialNums);
%             trialNums = trialNums(~isnan(trialNums));
%             trialNums = intersect(trialNums,b{b_ind}.trialNums); % try only the ones with behavior trials
% 
%             includef=cell(size(trialNums,1),1);
%             for i = 1: length(trialNums)
%                 includef{i} = num2str(trialNums(i));
%             end
%         end
% 
%         %%
%         tt_ind = cell(1,length(trial_types));
%         wl_array = cell(1,length(trial_types));
% 
%         for trial_type_num = 1 : length(trial_types)    
%             tt_ind{trial_type_num} = find(cellfun(@(x) strcmp(x.trialType,trial_types{trial_type_num}),b_session.trials));
%             temp_files = cell(length(tt_ind{trial_type_num}),1);
%             for j = 1 : length(tt_ind{trial_type_num})
%                 temp_files{j} = num2str(tt_ind{trial_type_num}(j));
%             end
%             wl = Whisker.WhiskerTrialLiteArray_2pad(whisker_d,'include_files',temp_files);
%             wl_array{trial_type_num} = wl;
%         end
%         %%
%         close all
%         base_angle = 21;
%         sp = cell(5,1);% subplot {1~5}{1~5} / first: rc rf lc lf free-whisking / second: subplot 1 ~ 5
%         for i = 1 : 5
%             sp{i} = cell(5,1);
%         end
%         for wl_array_ind = 1 : length(trial_types)
% %             figure, set(gcf,'Units','centimeters','Position',[wl_array_ind+3 wl_array_ind+3 11.7 9.5]); hold all
%             switch wl_array_ind
%                 case 1
%                     tt_text = 'Close Down';
%                 case 2
%                     tt_text = 'Far Down';
%                 case 3
%                     tt_text = 'Close Up';
%                 case 4
%                     tt_text = 'Far Up';
%             end
% %             title([tt_text 'touch frames'])
%             
%             for wl_ind = 1 : length(wl_array{wl_array_ind}.trials)
%                 wt = wl_array{wl_array_ind}.trials{wl_ind};      
%                 fw_frames = find(wt.time{1} < 0.0032 * 500); % first 500 frames
%                 sp{wl_array_ind}{1} = [sp{wl_array_ind}{1}; wt.thetaAtBase{1}(wt.th_touch_frames)'+base_angle, wt.thetaAtBase{2}(wt.th_touch_frames)'];
%                 sp{wl_array_ind}{2} = [sp{wl_array_ind}{2}; wt.deltaKappa{1}(wt.th_touch_frames)', wt.deltaKappa{2}(wt.th_touch_frames)'];
%                 sp{wl_array_ind}{3} = [sp{wl_array_ind}{3}; wt.thetaAtBase{1}(wt.th_touch_frames)'+base_angle, wt.deltaKappa{1}(wt.th_touch_frames)'];
%                 sp{wl_array_ind}{4} = [sp{wl_array_ind}{4}; wt.thetaAtBase{2}(wt.th_touch_frames)', wt.deltaKappa{2}(wt.th_touch_frames)'];
%                 sp{wl_array_ind}{5} = [sp{wl_array_ind}{5}; wt.thetaAtBase{1}(wt.th_touch_frames)'+base_angle, wt.deltaKappa{2}(wt.th_touch_frames)'];        
%                 sp{5}{1} = [sp{5}{1}; wt.thetaAtBase{1}(fw_frames)'+base_angle, wt.thetaAtBase{2}(fw_frames)'];
%                 sp{5}{2} = [sp{5}{2}; wt.deltaKappa{1}(fw_frames)', wt.deltaKappa{2}(fw_frames)'];
%                 sp{5}{3} = [sp{5}{3}; wt.thetaAtBase{1}(fw_frames)'+base_angle, wt.deltaKappa{1}(fw_frames)'];
%                 sp{5}{4} = [sp{5}{4}; wt.thetaAtBase{2}(fw_frames)', wt.deltaKappa{2}(fw_frames)'];
%                 sp{5}{5} = [sp{5}{5}; wt.thetaAtBase{1}(fw_frames)'+base_angle, wt.deltaKappa{2}(fw_frames)'];        
%             end
% 
%         end
%         spg{sgind}{sind} = sp;        
%     end
% end
% %%
% close all
% for i = 1 : 5 % different figure
%     switch i
%         case 1
%             tt_text = 'Close Down';
%         case 2
%             tt_text = 'Far Down';
%         case 3
%             tt_text = 'Close Up';
%         case 4
%             tt_text = 'Far Up';
%         case 5
%             tt_text = 'Free Whisking';
%     end
%     figure,set(gcf,'Units','centimeters','Position',[i+3 i+3 11.7 9.5]);     
%     for j = 1 : 5 % different subplot
%         subplot(3,2,j), hold on
%         for k = 1 : length(sessions{1})
%             scatter(spg{1}{k}{i}{j}(:,1),spg{1}{k}{i}{j}(:,2),0.5,'k.', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha',0.2); % , 'LineWidth', 0.001
%         end
%         for k = 1 : length(sessions{2})
%             scatter(spg{2}{k}{i}{j}(:,1),spg{2}{k}{i}{j}(:,2),0.5, 'b.', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha',0.2); % , 'LineWidth', 0.001, 
%         end
%         switch j
%             case 1
%                 title({['                                                       ',tt_text,' touch frames'];'\theta_T vs \theta_F'});
%                 xlim([-30 30]), ylim([-10 30]);
%             case 2
%                 title('\kappa_T vs \kappa_F');
%                 xlim([-0.01 0.01]), ylim([-0.01 0.01]);
%             case 3
%                 title('\theta_T vs \kappa_T');
%                 xlim([-30 30]), ylim([-0.01 0.01]);
%             case 4
%                 title('\theta_F vs \kappa_F');
%                 xlim([-10 30]), ylim([-0.01 0.01]);
%             case 5
%                 title('\theta_T vs \kappa_F');
%                 xlim([-30 30]), ylim([-0.01 0.01]);
%         end
%     end
% end


        
        