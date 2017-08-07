function phase_plot(handle, result, title_name, ...
    x_label_name, x_tick, x_tick_label, ...
    y_label_name, y_tick, y_tick_label, ...
    CLim_value)
% Make phase pcolor plot for a matrix
if ~exist('CLim_value')
    CLim_value = [-120, 10];
end
pcolor(result);
%set(handle, 'CLim', CLim_value);
caxis(CLim_value); % works for matlab 2017
shading interp;
colormap();
colorbar();
title(title_name);
xlabel(x_label_name,'Interpreter','LaTex');
ylabel(y_label_name,'Interpreter','LaTex');
xticks(x_tick);
xticklabels(strread(num2str(x_tick_label),'%s'));
yticks(y_tick);
yticklabels(strread(num2str(y_tick_label),'%s'));
%set(handle, 'XTick', x_tick, 'XTickLabel', strread(num2str(x_tick_label),'%s'),...
%     'YTick', y_tick, 'YTickLabel', strread(num2str(y_tick_label),'%s'))
end