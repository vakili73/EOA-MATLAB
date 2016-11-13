function [OPTIONS] = ComputeRandomShift(OPTIONS)
if OPTIONS.ShiftFlag > 0
    OPTIONS.OrderDependent = true; % the order of the independent variables does matter
    MaxShiftMagnitude = (OPTIONS.MaxDomain - OPTIONS.MinDomain) / 2;
    if OPTIONS.ShiftFlag == 1
        % use a uniform random number to shift the solution of the problem
        OPTIONS.ShiftAmount = MaxShiftMagnitude .* (2 * rand(1, OPTIONS.numVar) - 1);
    else
        % use a Gaussian random number to shift the solution of the problem
        sigma = (OPTIONS.MaxDomain - OPTIONS.MinDomain) / 6;
        OPTIONS.ShiftAmount = sigma .* randn(1, OPTIONS.numVar);
        OPTIONS.ShiftAmount(OPTIONS.ShiftAmount < -MaxShiftMagnitude) = -MaxShiftMagnitude(OPTIONS.ShiftAmount < -MaxShiftMagnitude);
        OPTIONS.ShiftAmount(OPTIONS.ShiftAmount > MaxShiftMagnitude) = MaxShiftMagnitude(OPTIONS.ShiftAmount > MaxShiftMagnitude);
    end
else
    OPTIONS.ShiftAmount = zeros(1, OPTIONS.numVar);
end
return