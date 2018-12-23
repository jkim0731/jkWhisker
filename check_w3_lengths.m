mice = {'JK025','JK027','JK030','JK036','JK037','JK038','JK039','JK041','JK052', 'JK053','JK054','JK056'};

videoloc = 'Y:\Whiskernas\JK\whisker\tracked\';
if strcmp(videoloc(end),filesep)
    whisker_d = videoloc;
else
    whisker_d = ([videoloc filesep]);
end

sessions = {    [1:19,22],  [1:22], [1:22], [1:18], [1:10,12:24],   [1:22,24:31],   [1:25], [1:19,21:30],   [1:29],     [1:3,5:21],     [1:26], [1:13]};
sessionsPre = { [1],        [1],    [1:2],  [1],    [1:2],          [1],            [1],    [1],            [1],        [1],            [1],    [1]};

errorFlag = zeros(length(mice), 31);

for mi = 10 : length(mice)
    parfor si = 1 : length(sessions{mi})
        w3a = Whisker.Whisker3D_2padArray(sprintf('%s%sS%02d',videoloc, mice{mi}, sessions{mi}(si)));
        errorFlag(mi,si) = sum(cellfun(@(x) ~isequal(length(x.time), length(x.kappaH), length(x.theta), length(x.kappaV), length(x.phi)), w3a.trials));
    end
%     parfor spi = 1 : length(sessionsPre{mi})
%         w3a = Whisker.Whisker3D_2padArray(sprintf('%s%spre%d',videoloc, mice{mi}, sessionsPre{mi}(spi)));
%         errorFlag(mi,32-spi) = sum(cellfun(@(x) ~isequal(length(x.time), length(x.kappaH), length(x.theta), length(x.kappaV), length(x.phi)), w3a.trials));
%     end
end
