close all
SetPlotOptions
x = -2 : .1 : 4;
plot(x, (x-1).^2-4)
xlabel('{\it K}');
ylabel('$$\mbox{Discriminant } \Delta$$', 'Interpreter', 'latex');
text('Interpreter', 'latex', 'String', '$$K_1 = \frac{\phi_T+1-2\sqrt{\phi_T}}{(\phi_T-1)^2}$$', ...
    'Position', [-1, -3], 'FontSize', 14)
text('Interpreter', 'latex', 'String', '$$K_2 = \frac{\phi_T+1+2\sqrt{\phi_T}}{(\phi_T-1)^2}$$', ...
    'Position', [3, -3], 'FontSize', 14)
set(gca, 'XTick', [-1, 3])
set(gca, 'XTickLabel', []);
set(gca, 'YTick', 0)
grid