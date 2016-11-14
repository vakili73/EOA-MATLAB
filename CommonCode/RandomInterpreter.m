function [ fcn ] = RandomInterpreter( chaoticMethod, params )
% params = '1, OPTIONS.numVar'
% 'rand()-rand()' convert to "rand(1, OPTIONS.numVar)-rand(1,
% OPTIONS.numVar)"

fcn = strrep(chaoticMethod, '()', ['(', params, ')']);

end

