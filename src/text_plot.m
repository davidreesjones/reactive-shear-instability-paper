function tt = text_plot(ax,Start,str)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pos=get(ax,'Position');
xlim=get(ax,'XLim');
ylim=get(ax,'YLim');

Start
xlim
tt   = text((Start(1) + abs(min(xlim)))/diff(xlim) * pos(3) + pos(1),...
(Start(2) - min(ylim))/diff(ylim) * pos(4) + pos(2),...
str);
tt.Interpreter='latex';
end

