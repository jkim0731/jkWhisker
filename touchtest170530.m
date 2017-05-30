%%
tid = [0 1]; % Set trajectory ID to view
wl = Whisker.WhiskerTrialLiteArray_2pad(whisker_d);
Whisker.viewdouble_WhiskerTrialLiteArray(wl,tid)
%%
ttf_ind = [];
for i = 1 : length(wl.trials)
    if ~isempty(wl.trials{i}.th_touch_frames)
        ttf_ind = [ttf_ind; i];
    end
end

%%
ttf_fn = cell(length(ttf_ind),1);
for i = 1 : length(ttf_ind)
    ttf_fn{i} = wl.trials{ttf_ind(i)}.trackerFileName;
end

%%
t_ind = 8;
confirmed_touch_frame = 551:556;
figure, 
% plot(wl.trials{t_ind}.intersect_coord(wl.trials{t_ind}.pole_available_timepoints,1), wl.trials{t_ind}.intersect_coord(wl.trials{t_ind}.pole_available_timepoints,2), 'b.'), hold on,
plot(wl.trials{t_ind}.th_polygon(:,1), wl.trials{t_ind}.th_polygon(:,2), 'k.'), hold on,
plot(wl.trials{t_ind}.intersect_coord(confirmed_touch_frame,1), wl.trials{t_ind}.intersect_coord(confirmed_touch_frame,2), 'b.')