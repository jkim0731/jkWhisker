function makeAllDirectory_WhiskerFinal_2pad(d)

cd(d)

wlflist = dir('*_WL_2pad.mat');
wltn = zeros(length(wlflist),1);
for i = 1 : length(wlflist)
    wltn(i) = str2double(strtok(wlflist(i).name,'_'));
end

w3flist = dir('*_W3_2pad.mat');
w3tn = zeros(length(w3flist),1);
for i = 1 : length(w3flist)
    w3tn(i) = str2double(strtok(w3flist(i).name,'_'));
end

tn = intersect(wltn, w3tn);

temp = strsplit(d,filesep);
sessionName = temp{end};

parfor i = 1 : length(tn)
    fprintf('Processing %d/%d from %s\n', i, length(tn), sessionName)
    wldat = load([num2str(tn(i)), '_WL_2pad.mat']);    
    w3dat = load([num2str(tn(i)), '_W3_2pad.mat']);
    wf = Whisker.WhiskerFinal_2pad(wldat.wl, w3dat.w3);
    savefn = [num2str(tn(i)), '_WF_2pad.mat'];
    pctsave(savefn,wf);
end

end

function pctsave(outfn,wf)
save(outfn,'wf');
end
