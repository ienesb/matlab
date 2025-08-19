function [x_detected, threshold] = my_cfar(x)
    numGuardCells = 2;
    numTrainingCells = 200;
    numPadding = (numTrainingCells + numGuardCells) / 2;
    ProbabilityFalseAlarm = 1e-3;
    alpha = numTrainingCells * (ProbabilityFalseAlarm^(-1/numTrainingCells) - 1);
    x_padded = [x(end-numPadding+1:end) ; x; x(1: numPadding)];
    
    x_detected = zeros(size(x));

    threshold = zeros(size(x));

    for idx = 1:length(x)
        i1 = idx;
        i2 = idx+numTrainingCells/2-1;
        i3 = idx+numPadding+numGuardCells/2+1;
        i4 = idx+numPadding+numGuardCells/2+numTrainingCells/2;

        Pn = (sum(x_padded(i1:i2)) + sum(x_padded(i3:i4)))/numTrainingCells;
        T = alpha * Pn;
        x_detected(idx) = x(idx) >= T;
        threshold(idx) = T;
    end
end