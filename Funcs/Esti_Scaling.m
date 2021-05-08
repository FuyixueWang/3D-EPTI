function [scaling] = Esti_Scaling(kdata)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    tmp = dimnorm(ifft2c(kdata), 3);
    tmpnorm = dimnorm(tmp, 4);
    tmpnorm2 = sort(tmpnorm(:), 'ascend');
    % match convention used in BART
    p100 = tmpnorm2(end);
    p90 = tmpnorm2(round(.9 * length(tmpnorm2)));
    p50 = tmpnorm2(round(.5 * length(tmpnorm2)));
    if (p100 - p90) < 2 * (p90 - p50)
        scaling = p90;
    else
        scaling = p100;
    end
end

