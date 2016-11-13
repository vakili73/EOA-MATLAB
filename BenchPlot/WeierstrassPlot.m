function WeierstrassPlot

xlimits = [-0.5 0.5];
MaxSumIndex = 20;
a = 0.5;
b = 3;
dx = (xlimits(2) - xlimits(1)) / 200;
N = floor((xlimits(2) - xlimits(1)) / dx) + 1;
obj = zeros(N, N);
i = 0;
for x = xlimits(1) : dx : xlimits(2)
    i = i + 1;
    k = 0;
    for y = xlimits(1) : dx : xlimits(2)
        k = k + 1;
        sumx1 = 0;
        sumy1 = 0;
        sum2 = 0;
        for m = 0 : MaxSumIndex
            sumx1 = sumx1 + a^m * cos(2 * pi * b^m * (x + 0.5));
            sumy1 = sumy1 + a^m * cos(2 * pi * b^m * (y + 0.5));
            sum2 = sum2 + a^m * cos(pi * b^m);
        end
        obj(i,k) = sumx1 + sumy1 - 2 * sum2;
    end
end
close all;
[x, y] = meshgrid(xlimits(1) : dx : xlimits(2));
meshc(x,y,obj);
zlabel('f(x)')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
temp = axis;
axis([xlimits, xlimits, temp(5), max(max(obj))])