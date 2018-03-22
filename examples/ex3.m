p = .2;
N = 40000;
z = randn(2, N);
eta = (rand(2, N) > p);
x = z.*eta;
[X, Y] = meshgrid(-2:.01:2, -2:.01:2);
Z = X;
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        Z(i,j) = mean(mean( (x(1, :) - X(i,j)).^2 + (x(2, :) - Y(i,j)).^2 < .001));
    end
end
%%
mesh(X, Y, Z)
xlabel('x')
ylabel('y')
zlabel('density estimation using kernel')