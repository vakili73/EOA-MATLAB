function [ Population ] = PopAppend(Population, Candidate)
%POPAPPEND Summary of this function goes here
%   Detailed explanation goes here

i = length(Population) + 1;
for j = 1:length(Candidate)
    Population(i) = Candidate(j);
    i = i + 1;
end

end