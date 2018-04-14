function [hh, amplitude, filteredSignal, midpoint, amplitudeS, midpointS, phase, phaseS] =  SAHWhiskerDecomposition(theta)
%% calculate whisker amplitude and midpoint


sampleRate=  310;

% make any nan thetaAtBase = mean of the surrounding points (10 on each
% side)
try
    theta(isnan(theta)) = nanmean(theta(repmat(find(isnan(theta)),21,1)+repmat([-10:10]',1,length(find(isnan(theta))))));
catch
    theta(isnan(theta)) = nanmean(theta);
end


BandPassCutOffsInHz = [8 30];  %%check filter parameters!!!
% From Sofroniew 2014, which sites Hill 2011
W1 = BandPassCutOffsInHz(1) / (sampleRate/2);
W2 = BandPassCutOffsInHz(2) / (sampleRate/2);
[b,a]=butter(2,[W1 W2]);
filteredSignal = filtfilt(b, a, theta);

[b,a]=butter(2, BandPassCutOffsInHz(1)/ (sampleRate/2),'low');
new_midpoint = filtfilt(b,a,theta-filteredSignal);
hh=hilbert(filteredSignal);

amplitude=abs(hh);
midpoint=new_midpoint;

amplitudeS=smooth(amplitude,10,'moving');
% amplitudeDS=amplitudeS(5:10:end);

midpointS=smooth(midpoint,10,'moving');
% midpointDS=midpointS(5:10:end);

phase = angle(hh);
phaseS = smooth(phase,10,'moving');
% phaseDS = phaseS(5:10:end);
 
%%


