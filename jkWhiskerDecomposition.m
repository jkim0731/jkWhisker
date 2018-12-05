function [hh, amplitude, filteredSignal, midpoint, amplitudeS, midpointS, phase, phaseS] =  jkWhiskerDecomposition(theta, varargin)
%% calculate whisker amplitude and midpoint
% varargin{1} = time;
% varargin{2} = sampleRate;

sampleRate = 311;
time = [0:length(theta)-1]/sampleRate;

if nargin == 2
    time = varargin{1};
    if length(time) ~= length(theta)
        error('Theta and time should be same length')
    end
    sampleRate = 1/min(diff(time));
elseif nargin == 3
    sampleRate = varargin{2};
end

if abs(max(diff(time)) - 1/sampleRate) > 0.5/sampleRate
    targetTime = 0:1/sampleRate:max(time);
    inputTheta = theta;
    theta = nan(length(targetTime),1);
    inds = zeros(length(time),1);
    for i = 1 : length(inds)
        [~, inds(i)] = min(abs(targetTime - time(i)));
    end
    theta(inds) = inputTheta;
end
% make any nan thetaAtBase = mean of the surrounding points (10 on each side)
try
    theta(isnan(theta)) = nanmean(theta(repmat(find(isnan(theta)),21,1)+repmat([-10:10]',1,length(find(isnan(theta))))));
catch
    theta(isnan(theta)) = nanmean(theta);
end

BandPassCutOffsInHz = [6 30];  %%check filter parameters!!!
% From Sofroniew 2014, which sites Hill 2011
W1 = BandPassCutOffsInHz(1) / (sampleRate/2);
W2 = BandPassCutOffsInHz(2) / (sampleRate/2);
[b,a]=butter(2,[W1 W2]);
filteredSignal = filtfilt(b, a, theta);

[b,a] = butter(2, BandPassCutOffsInHz(1) / (sampleRate/2),'low');
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


