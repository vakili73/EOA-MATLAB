clc;
clear;
close all;

% Initializing
% use this function for generate random number from clock puls:
% fix(sum(100*clock));
Init__;

%% Start the program

chaoticPlot = dir('ChaoticPlot');

for i = 3 : length(chaoticPlot)
    temp = chaoticPlot(i).name;
    feval(temp(1, 1:length(temp)-2));
end

%% bernolli()
Bernoulli
Bernoulli(1, 5)
Bernoulli(5)