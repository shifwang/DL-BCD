function drawLineSeg(x,y)
% draw a line segment between x and y
hori = x(1):0.5:y(1);
vert = x(2):0.5:y(2);
hold on;
plot(hori,vert,'color','red');
hold off;
end