function Penalty1Plot

xlimits = [-20 20];
dx = (xlimits(2) - xlimits(1)) / 100;
N = floor((xlimits(2) - xlimits(1)) / dx) + 1;
obj = zeros(N, N);
i = 0;
r = rand;
for x = xlimits(1) : dx : xlimits(2)
    i = i + 1;
    k = 0;
    for y = xlimits(1) : dx : xlimits(2)
        k = k + 1;
        if x > 10
            u1 = 100*(x-10)^4;
        elseif x < -10
            u1 = 100*(-x-10)^4;
        else
            u1 = 0;
        end
        if y > 10
            u2 = 100*(y-10)^4;
        elseif y < -10
            u2 = 100*(-y-10)^4;
        else
            u2 = 0;
        end
        y1 = 1 + (x + 1) / 4;
        y2 = 1 + (y + 1) / 4;
        obj(i,k) = pi/2 * (10*(sin(pi*y1))^2 + (y1-1)^2 * (1 + 10*(sin(pi*y2))^2) + (y2 - 1)^2) + u1 + u2;
    end
end
obj = obj / 1000000;
close all;
[x, y] = meshgrid(xlimits(1) : dx : xlimits(2));
meshc(x,y,obj);
zlabel('f(x)')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
temp = axis;
axis([xlimits, xlimits, temp(5), max(max(obj))])