function [diff1, diff2] = when_number_of_files_dont_match

%%
diff1 = [];
diff2 = [];

mp4 = dir('*.mp4');
whiskers = dir('*.whiskers');
measurements = dir('*.measurements');

mp4_num = zeros(length(mp4),1);
whiskers_num = zeros(length(whiskers),1);
measurements_num = zeros(length(measurements),1);

for i = 1 : length(mp4)
    mp4_num(i) = str2double(mp4(i).name(1:end-4));
end

for i = 1 : length(whiskers)
    whiskers_num(i) = str2double(whiskers(i).name(1:end-9));
end

for i = 1 : length(measurements)
    measurements_num(i) = str2double(measurements(i).name(1:end-13));
end

if length(mp4_num) ~= length(whiskers_num)
    diff1 = setdiff(mp4_num,whiskers_num);
    disp('mp4 and whiskers differences are in:')
    for i = 1 : length(diff1)
        disp(num2str(diff1(i)));
    end
end

if length(whiskers_num) ~= length(measurements_num)
    diff2 = setdiff(whiskers_num,measurements_num);
    disp('whiskers and measurements differences are in:')
    for i = 1 : length(diff2)
        disp(num2str(diff2(i)));
    end
end
