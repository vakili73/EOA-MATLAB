function [m] = createRotMatrix(D)
% Create a random rotation matrix "m" with "D" dimensions
[m, ~] = qr(randn(D));
if det(m) < 1
    m = -m;
end
n = m*m' - eye(D);
f = norm(n, 'fro');
if f > 0.00001
    error('Error creating rotation matrix')
end
return