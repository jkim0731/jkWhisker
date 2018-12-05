function [onsetFrame, amplitude, midpoint, whiskingAmp] = jkWhiskerOnsetNAmplitude(theta, varargin)

% calculate whisking onset frames calculated from phase and max amplitude
% of that whisking bout.
% input is assumed to be continuous frames with equal inter-frame interval.
% If not tracked, fill in NaN.
switch nargin
    case 1
        whiskingThreshold = 5; % in degrees
        sampleRate = 311;
    case 2        
        whiskingThreshold = varargin{1};
        sampleRate = 311;
    case 3
        whiskingThreshold = varargin{1};
        sampleRate = varargin{2};
    otherwise
        error('too much input arguments')
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
phase = angle(hh);

onsetCandid = [find(diff(phase) < -5) + 1; length(phase)];
whiskingAmp = zeros(length(onsetCandid)-1,1);
for i = 1 : length(onsetCandid)-1
    whiskingAmp(i) = max(amplitude(onsetCandid(i):onsetCandid(i+1)));
end
onsetFrame = onsetCandid((whiskingAmp > whiskingThreshold));
