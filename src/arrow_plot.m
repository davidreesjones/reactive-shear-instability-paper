function ar = arrow_plot(ax,Start,End)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pos=get(ax,'Position');
xlim=get(ax,'XLim');
ylim=get(ax,'YLim');

ar   = annotation('arrow',...
[(Start(1) + abs(min(xlim)))/diff(xlim) * pos(3) + pos(1),...
(End(1) + abs(min(xlim)))/diff(xlim) * pos(3) + pos(1) ],...
[(Start(2) - min(ylim))/diff(ylim) * pos(4) + pos(2),...
(End(2) - min(ylim))/diff(ylim) * pos(4) + pos(2)]);
ar.HeadStyle='plain';
ar.HeadLength=6;
ar.HeadWidth=6;
%ar.Color='c'
ar.LineStyle='none';
end

