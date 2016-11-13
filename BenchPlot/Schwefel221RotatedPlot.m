function Schwefel221RotatedPlot

xlimits = [-3 3];
dx = (xlimits(2) - xlimits(1)) / 100;
N = floor((xlimits(2) - xlimits(1)) / dx) + 1;
MRot = createRotMatrix(2); % Random rotation matrix
% Create cell array of rotated independent coordinates
yCell = cell(1, N^2);
CellIndex = 0;
for i = 1 : N
    x1 = xlimits(1) + (i - 1) * dx;
    for j = 1 : N
        x2 = xlimits(1) + (j - 1) * dx;
        xRot = [x1, x2] * MRot;
        CellIndex = CellIndex + 1;
        yCell{CellIndex} = xRot;
    end
end
% Calculate rotated Schwefel 2.21 function
objRot = zeros(N, N);
i = 1;
k = 1;
for j = 1 : N^2
    x = yCell{j};
    objRot(i, k) = max(abs(x(1)), abs(x(2)));
    k = k + 1;
    if k > N
        k = 1;
        i = i + 1;
    end
end
% Calculated non-rotated Schwefel 2.21 function
obj = zeros(N, N);
i = 0;
for x = xlimits(1) : dx : xlimits(2)
    i = i + 1;
    k = 0;
    for y = xlimits(1) : dx : xlimits(2)
        k = k + 1;
        obj(i,k) = max(abs(x), abs(y));
    end
end
close all
figure
[x, y] = meshgrid(xlimits(1) : dx : xlimits(2));
% Plot non-rotated function
subplot(1, 2, 1), meshc(x, y, obj);
zlabel('f(x)')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
temp = axis;
axis([xlimits, xlimits, temp(5), max(max(obj))])
% Plot rotated function
subplot(1, 2, 2), meshc(x, y, objRot);
zlabel('f(x)')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
temp = axis;
axis([xlimits, xlimits, temp(5), max(max(objRot))])