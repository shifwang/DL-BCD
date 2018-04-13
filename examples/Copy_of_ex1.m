N = 1000;
z = [[10 1; -10, 1] * randn(2, N/2), [10 1; 10, -1] * randn(2, N/2)];
p = .3;
eta = (rand(2, N) > p);
x = z.*eta;
trueDict = eye(2);
true_coef = x;
subplot(1,3,1);
scatter(x(1,:), x(2,:), 4, 'MarkerFaceColor', 'b',...
    'MarkerEdgeColor', 'b',...
    'MarkerEdgeAlpha', .5,...
    'MarkerFaceAlpha', .5)
axis([-30, 30, -30,30])
title('raw data points')
hold on;
len = 20;
plot([-len, len], [0, 0], 'r', 'LineWidth', 3)
plot([0, 0], [-len, len], 'r', 'LineWidth', 3)
plot([-len, len], [-len, len], 'g', 'LineWidth', 3)
plot([-len, len], [len, -len], 'g', 'LineWidth', 3)

%%
% plot the surface of primal formulation for 2 by 2 dictionary
% stand alone, no dependencies
K = 2;
theta1 = -pi:(pi/21):pi;
theta2 = (-pi):(pi/22):0;
obj = nan(length(theta1),length(theta2));
s = 1;
for i = 1:length(theta1)
    for j = 1:length(theta2)
        if i ~= j && i + j ~= 0
            dict = [sin(theta1(i)) sin(theta2(j));...
                  cos(theta1(i)) cos(theta2(j))];  
                tmp = dict\x;
            obj(i,j) = 1/K*sum(sum(abs(tmp)))/N;
        end
    end
end
%% 3d mesh 
subplot(1, 3, 2)
surf(theta1,theta2,log10(min(obj', 9)));

%xlabel('Angle between two atoms');
%ylabel('Objective');
hold on;
view(-15,10)
plot3(0, -pi/2, log10(mean(mean(abs(x))))-0.0,'r.','MarkerSize',30);
text(0, -pi/2, log10(mean(mean(abs(x))))-0.0,'D^*');
[~, ind] = min(obj(:));
xxx = floor(ind/size(obj,1));
yyy = ind - size(obj,1)*xxx;
plot3(theta1(yyy), theta2(xxx), log10(obj(ind))-.0,'g.','MarkerSize',30);
text(theta1(yyy), theta2(xxx), log10(obj(ind))-0.0,'global min');
hold off;
%[~, trueInd] = min(abs(cos(theta*pi/180) - mu));
%semilogy(acos(mu)/pi*180,1,'r.','MarkerSize',30)
%hold off;
%saveas(gcf,strcat('contour/ddd',num2str(ii),'.jpg'))
xlabel('\theta_1');
ylabel('\theta_2');
%% 2d contour
subplot(1,3,2);
contour(theta1,theta2,min(log10(obj'),2),10);
title('Objective contour plot')

hold on;
plot(0, -pi/2,'r.','MarkerSize',30)
plot(pi/2, -pi,'r.','MarkerSize',30)
plot(pi, -pi/2,'r.','MarkerSize',30)
plot(-pi/2, -pi,'r.','MarkerSize',30)
plot(-pi/2, 0,'r.','MarkerSize',30)
plot(-pi, -pi,'r.','MarkerSize',30)
axis([-pi, pi, -pi, 0])
text(5, -pi/2, 'D^*', 'FontWeight', 'bold')
plot(theta1(yyy), theta2(xxx),'g.','MarkerSize',30)
plot(theta1(yyy)+pi/2, theta2(xxx)+pi/2,'g.','MarkerSize',30)
plot(theta1(yyy)+pi, theta2(xxx),'g.','MarkerSize',30)
plot(theta1(yyy)-pi/2, theta2(xxx) + pi/2,'g.','MarkerSize',30)
text(theta1(yyy), theta2(xxx), '  global min', 'FontWeight', 'bold')
plot(acos(trueDict(1,2)), acos(trueDict(2,2)),'r.','MarkerSize',30)
xlabel('\theta_1')
ylabel('\theta_2')
hold off;
%saveas(gcf,strcat('dual/dual',num2str(ii),'.jpg'))
%% geometric curve
subplot(1,3,3)
n = size(obj, 2);
obj_line = diag(obj(91:(n+90), 1:n));
plot(theta2(30:119), obj_line(30:119))
axis([min(theta2(30:119)), max(theta2(30:119)),min(obj_line(30:119)) - .1, max(obj_line(30:119))+ .1])
xlabel('\theta_1')
ylabel('objective')
hold on;
plot(-90/180*pi,obj_line(45 + 45) ,'r.','MarkerSize',30)
text(-90/180*pi, obj_line(45 + 45), '  ref dict', 'FontWeight', 'bold')
plot(-135/180*pi,obj_line(45) ,'g.','MarkerSize',30)
text(-135/180*pi, obj_line(45), '  global min', 'FontWeight', 'bold')
