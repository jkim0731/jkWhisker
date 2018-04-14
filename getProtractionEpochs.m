function epochs = getProtractionEpochs(theta)

if size(theta,1) > size(theta,2)
    theta = theta';
end

if length(theta) > 1
    diffTheta = diff(theta);
    if ~isempty(find(diffTheta > 0,1))
        [~, dTNZind, dTNZ] = find(diffTheta); % diff Theta Non Zero
        dTNZ = dTNZ ./ abs(dTNZ);
        troughsDTNZ = find(diff(dTNZ)>0)+1;
        if dTNZ(1) > 0 % starting with protraction
            troughsDTNZ = [1, troughsDTNZ];
        end
        peaksDTNZ = find(diff(dTNZ) < 0) + 1;

        epochs = cell(length(troughsDTNZ),1);
        for edi = 1 : length(epochs)-1
            epochs{edi} = dTNZind(troughsDTNZ(edi)) : dTNZind(peaksDTNZ(edi)); 
        end
        if length(peaksDTNZ) == length(troughsDTNZ)
            epochs{end} = dTNZind(troughsDTNZ(end)) : dTNZind(peaksDTNZ(end));
        else
            epochs{end} = dTNZind(troughsDTNZ(end)) : length(theta);
        end
    else
        epochs = {};
    end
else
    epochs = {};
end