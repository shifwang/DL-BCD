N = 2000;
z = [[10 1; -10, 1] * randn(2, N/2), [10 1; 10, -1] * randn(2, N/2)];
p = .3;
eta = (rand(2, N) > p);
x = z.*eta;
trueDict = eye(2);
true_coef = x;
subplot(1,2,1);
scatter(x(1,:), x(2,:), 4, 'MarkerFaceColor', 'b',...
    'MarkerEdgeColor', 'b',...
    'MarkerEdgeAlpha', .5,...
    'MarkerFaceAlpha', .5)
axis([-30, 30, -30,30])
title('raw data points')
%%
% plot the surface of primal formulation for 2 by 2 dictionary
% stand alone, no dependencies
K = 2;
theta1 = -179:1:179;
theta2 = -179:1:0;
obj = nan(length(theta1),length(theta2));
s = 1;
for i = 1:length(theta1)
    for j = 1:length(theta2)
        if i ~= j && i + j ~= 0
            dict = [sin(theta(i)*pi/180) sin(theta(j)*pi/180);...
                  cos(theta(i)*pi/180) cos(theta(j)*pi/180)];  
                tmp = dict\x;
            obj(i,j) = 1/K*sum(sum(abs(tmp)))/N;
        end
    end
end
%% 3d mesh 
mesh(theta1,theta2,log10(min(obj', 10)));
%xlabel('Angle between two atoms');
%ylabel('Objective');
hold on;
view(-30,45)
plot3(0, -90, log10(mean(mean(abs(x))))+0.01,'r.','MarkerSize',30);
[~, ind] = min(obj(:));
xxx = floor(ind/size(obj,1));
yyy = ind - size(obj,1)*xxx;
plot3(theta1(yyy), theta2(xxx), log10(obj(ind))+.01,'g.','MarkerSize',30);
hold off;
%[~, trueInd] = min(abs(cos(theta*pi/180) - mu));
%semilogy(acos(mu)/pi*180,1,'r.','MarkerSize',30)
%hold off;
%saveas(gcf,strcat('contour/ddd',num2str(ii),'.jpg'))
%% 2d contour
subplot(1,2,2);
contour(theta1,theta2,min(log10(obj'),2),40);
title('Objective contour plot')
xlabel('Angle between 1st atom and x-axis (-180, 180)');
ylabel('Angle between 2nd atom and x-axis (-180, 0)');
hold on;
plot(0, -90,'r.','MarkerSize',30)
plot(90, -180,'r.','MarkerSize',30)
plot(180, -90,'r.','MarkerSize',30)
plot(-90, -180,'r.','MarkerSize',30)
plot(-90, 0,'r.','MarkerSize',30)
plot(-180, -90,'r.','MarkerSize',30)
text(5, -90, 'D^*', 'FontWeight', 'bold')
plot(theta1(yyy), theta2(xxx),'g.','MarkerSize',30)
plot(theta1(yyy)-90, theta2(xxx)-90,'g.','MarkerSize',30)
plot(theta1(yyy)+90, theta2(xxx)-90,'g.','MarkerSize',30)
plot(theta1(yyy)-180, theta2(xxx),'g.','MarkerSize',30)
text(theta1(yyy), theta2(xxx), '  global min', 'FontWeight', 'bold')
plot(acos(trueDict(1,2))*180/pi, acos(trueDict(2,2))*180/pi,'r.','MarkerSize',30)
hold off;
%saveas(gcf,strcat('dual/dual',num2str(ii),'.jpg'))


