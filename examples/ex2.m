N = 2000;
z = abs(randn(2, N));
p = .5;
eta = (rand(2, N) > p);
true_coef = z.*eta;
mu = .9;
trueDict = [1, mu; mu, 1]^.5;
x = trueDict * true_coef;
subplot(1,3,1);
scatter(x(1,:), x(2,:), 4, 'MarkerFaceColor', 'b',...
    'MarkerEdgeColor', 'b',...
    'MarkerEdgeAlpha', .5,...
    'MarkerFaceAlpha', .5)
%axis([-30, 30, -30,30])
title('raw data points')
%%
% plot the surface of primal formulation for 2 by 2 dictionary
% stand alone, no dependencies
K = 2;
theta1 = 0:.01:(pi/2 - .01);
theta2 = 0:.009:pi/2;
obj = nan(length(theta1),length(theta2));
s = 1;
for i = 1:length(theta1)
    for j = 1:length(theta2)
        if i ~= j && i + j ~= 0
            dict = [cos(theta1(i)) cos(theta2(j));...
                  sin(theta1(i)) sin(theta2(j))];  
                tmp = dict\x;
            obj(i,j) = 1/K*sum(sum(abs(tmp)))/N;
        end
    end
end
%% 2d contour
subplot(1,3,2);
contour(theta1,theta2,min(log10(obj'),2),40);
title('Objective contour plot')
xlabel('Angle between 1st atom and x-axis');
ylabel('Angle between 2nd atom and x-axis');
hold on;
plot(0, -90,'r.','MarkerSize',30)
plot(90, -180,'r.','MarkerSize',30)
plot(180, -90,'r.','MarkerSize',30)
plot(-90, -180,'r.','MarkerSize',30)
plot(-90, 0,'r.','MarkerSize',30)
plot(-180, -90,'r.','MarkerSize',30)
axis([-180, 180, -180, 0])
text(5, -90, 'D^*', 'FontWeight', 'bold')
plot(theta1(yyy), theta2(xxx),'g.','MarkerSize',30)
plot(theta1(yyy)+90, theta2(xxx)+90,'g.','MarkerSize',30)
plot(theta1(yyy)+180, theta2(xxx),'g.','MarkerSize',30)
plot(theta1(yyy)-90, theta2(xxx) + 90,'g.','MarkerSize',30)
text(theta1(yyy), theta2(xxx), '  global min', 'FontWeight', 'bold')
plot(acos(trueDict(1,2))*180/pi, acos(trueDict(2,2))*180/pi,'r.','MarkerSize',30)
hold off;
%saveas(gcf,strcat('dual/dual',num2str(ii),'.jpg'))
%% geometric curve
subplot(1,3,3)
n = size(obj, 2);
obj_line = diag(obj(91:(n+90), 1:n));
plot(theta2(30:119), obj_line(30:119))
axis([min(theta2(30:119)), max(theta2(30:119)),min(obj_line(30:119)) - .1, max(obj_line(30:119))+ .1])
xlabel('angle \theta')
ylabel('objective')
%hold on;
%plot(-90,obj_line(-90+180) ,'r.','MarkerSize',30)
text(-90, obj_line(-90+180), '  ref dict', 'FontWeight', 'bold')
%plot(-135,obj_line(-135 + 180) ,'g.','MarkerSize',30)
text(-135, obj_line(-135 + 180), '  global min', 'FontWeight', 'bold')
%% 3d mesh (outdated)
mesh(theta1,theta2,log10(min(obj', 10)));
%xlabel('Angle between two atoms');
%ylabel('Objective');
hold on;
view(-30,45)
%plot3(45 - acos(mu)/2*180/pi, 45 + acos(mu)/2*180/pi, log10(mean(mean(abs(x))))+0.01,'r.','MarkerSize',30);
%[~, ind] = min(obj(:));
%xxx = floor(ind/size(obj,1));
%yyy = ind - size(obj,1)*xxx;
%plot3(theta1(yyy), theta2(xxx), log10(obj(ind))+.01,'g.','MarkerSize',30);
hold off;
%[~, trueInd] = min(abs(cos(theta*pi/180) - mu));
%semilogy(acos(mu)/pi*180,1,'r.','MarkerSize',30)
%hold off;
%saveas(gcf,strcat('contour/ddd',num2str(ii),'.jpg'))
