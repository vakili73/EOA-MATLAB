function [x, y, obj] = AckleyContour(OPTIONS)
% Create a contour plot of the Ackley benchmark
xlimits = [OPTIONS.MinDomain(1) OPTIONS.MaxDomain(1)];
ylimits = [OPTIONS.MinDomain(2) OPTIONS.MaxDomain(2)];
dx = (xlimits(2) - xlimits(1)) / 100;
dy = (ylimits(2) - ylimits(1)) / 100;
Nx = floor((xlimits(2) - xlimits(1)) / dx) + 1;
Ny = floor((ylimits(2) - ylimits(1)) / dy) + 1;
obj = zeros(Nx, Ny);
i = 0;
for x = xlimits(1) : dx : xlimits(2)
    i = i + 1;
    k = 0;
    for y = ylimits(1) : dy : ylimits(2)
        k = k + 1;
        sum1 = x^2 + y^2;
        sum2 = cos(2*pi*x) + cos(2*pi*y);
        obj(i,k) = exp(1) - 20 * exp(-0.2*sqrt(sum1/2)) - exp(sum2/2);
    end
end
[x, y] = meshgrid(xlimits(1) : dx : xlimits(2), ylimits(1) : dy : ylimits(2));