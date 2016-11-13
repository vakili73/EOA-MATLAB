function SphereShiftedPlot

xlimits = [-50 50];
bias = [10 -10];
dx = (xlimits(2) - xlimits(1)) / 100;
N = floor((xlimits(2) - xlimits(1)) / dx) + 1;
objShifted = zeros(N, N);
obj = zeros(N, N);
i = 0;
for x = xlimits(1) : dx : xlimits(2)
    i = i + 1;
    k = 0;
    for y = xlimits(1) : dx : xlimits(2)
        k = k + 1;
        objShifted(i,k) = (x-bias(1))^2 + (y-bias(2))^2;
        obj(i,k) = x^2 + y^2;
    end
end
close all
figure
[x, y] = meshgrid(xlimits(1) : dx : xlimits(2));
% Plot the unshifted function
subplot(1,2,1), meshc(x,y,obj);
zlabel('f(x)')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
temp = axis;
axis([xlimits, xlimits, temp(5), max(max(obj))])
% Plot the shifted function
subplot(1,2,2), meshc(x,y,objShifted);
zlabel('f(x)')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
temp = axis;
axis([xlimits, xlimits, temp(5), max(max(obj))])