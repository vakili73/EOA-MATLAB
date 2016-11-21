function GADyn1

% Dynamic systems analysis of simple selection-only GA.
% Reeves & Rowe, p. 148.

f = [1 2 2 3 2 3 3 4]; % fitnesses
p = [.93 .01 .01 .01 .01 .01 .01 .01]'; % proportions of each individual in population
genMax = 15;
pSave = zeros(length(p), genMax+1);
pSave(:, 1) = p;
for gen = 1 : genMax
    p = diag(f) * p / (f * p);
    pSave(:, gen+1) = p;
end
close all; figure; set(gca, 'fontsize', 14); hold on; box on;
plot(0:genMax, pSave(1,:), 'r:', 'linewidth', 2);
plot(0:genMax, pSave(4,:), 'b-', 'linewidth', 2);
plot(0:genMax, pSave(8,:), 'k--', 'linewidth', 2);
xlabel('generation')
ylabel('proportion of population');
legend('p_1', 'p_4 = p_6 = p_7', 'p_8');