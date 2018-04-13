function drawLineSeg(x,y)
% draw a line segment between x and y
hori = linspace(x(1), y(1), 100);
vert = linspace(x(2), y(2), 100);
hold on;
%plot(hori,vert,'color','blue');
hold off;
end