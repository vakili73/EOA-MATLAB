function Penalty2Plot

xlimits = [-50 50];
dx = (xlimits(2) - xlimits(1)) / 100;
N = floor((xlimits(2) - xlimits(1)) / dx) + 1;
obj = zeros(N, N);
i = 0;
for x = xlimits(1) : dx : xlimits(2)
    i = i + 1;
    k = 0;
    for y = xlimits(1) : dx : xlimits(2)
        k = k + 1;
        if x > 5
            u1 = 100*(x-5)^4;
        elseif x < -5
            u1 = 100*(-x-5)^4;
        else
            u1 = 0;
        end
        if y > 5
            u2 = 100*(y-5)^4;
        elseif y < -5
            u2 = 100*(-y-5)^4;
        else
            u2 = 0;
        end
        obj(i,k) = 0.1 * ((sin(3*pi*x))^2 + (x-1)^2 * (1 + (sin(3*pi*y))^2) + (y - 1)^2 * (1 + (sin(2*pi*y))^2)) + u1 + u2;
    end
end
obj = obj / 100000000;
close all
[x, y] = meshgrid(xlimits(1) : dx : xlimits(2));
meshc(x,y,obj);
zlabel('f(x)')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
temp = axis;
axis([xlimits, xlimits, temp(5), max(max(obj))])