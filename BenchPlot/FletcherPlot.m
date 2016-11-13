function FletcherPlot

a = 200 * rand(2,2) - 100;
b = 200 * rand(2,2) - 100;
alph = 2 * pi * rand(2,1) - pi;
A = zeros(2, 1);
for i = 1 : 2
    for j = 1 : 2
        A(i) = A(i) + a(i,j) * sin(alph(j)) + b(i,j) * cos(alph(j));
    end
end

xlimits = [-pi pi];
dx = (xlimits(2) - xlimits(1)) / 100;
N = floor((xlimits(2) - xlimits(1)) / dx) + 1;
obj = zeros(N, N);
i = 0;
for x = xlimits(1) : dx : xlimits(2)
    i = i + 1;
    k = 0;
    B = 0;
    for y = xlimits(1) : dx : xlimits(2)
        k = k + 1;
        B(1) = a(1,1) * sin(x) + a(1,2) * sin(y) + b(1,1) * cos(x) + b(1,2) * cos(y);
        B(2) = a(2,1) * sin(x) + a(2,2) * sin(y) + b(2,1) * cos(x) + b(2,2) * cos(y);
        obj(i,k) = (A(1) - B(1))^2 + (A(2) - B(2))^2;
    end
end
obj = obj / 10000;
close all;
[x, y] = meshgrid(xlimits(1) : dx : xlimits(2));
meshc(x,y,obj);
zlabel('f(x)')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
temp = axis;
diff = max(max(obj)) - min(min(obj));
axis([xlimits, xlimits, min(min(obj)) - diff/5, max(max(obj))])
