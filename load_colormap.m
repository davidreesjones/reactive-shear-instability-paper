map=[103,0,31
178,24,43
214,96,77
244,165,130
253,219,199
247,247,247
209,229,240
146,197,222
67,147,195
33,102,172
5,48,97]/255;
map=map(end:-1:1,:);
cmap=zeros(100,3);
for J=1:3
cmap(:,J)=interp1(linspace(0,1,numel(map(:,J))),map(:,J),linspace(0,1,numel(cmap(:,J))));
end
cmap_red=cmap(50:100,1:3);
cmap_blue=cmap(50:-1:1,1:3);
cmap_rev=cmap(end:-1:1,1:3);
newcolours={'#FAEB82','#E9B665','#DB844A','#C05736','#6B3F30'};
colororder(newcolours)

map_BrBG=[11,112,104
    88,176,166
    179,226,219
    244,244,244
    240,223,178
    207,162,85
    153,93,18]/255;
map_BrBG=map_BrBG(end:-1:1,:);
cmap_BrBG=zeros(100,3);
for J=1:3
cmap_BrBG(:,J)=interp1(linspace(0,1,numel(map_BrBG(:,J))),map_BrBG(:,J),linspace(0,1,numel(cmap_BrBG(:,J))));
end
cmap_BG=cmap_BrBG(50:100,1:3);
cmap_Br=cmap_BrBG(50:-1:1,1:3);
cmap_BrBG_rev=cmap_BrBG(end:-1:1,1:3);
