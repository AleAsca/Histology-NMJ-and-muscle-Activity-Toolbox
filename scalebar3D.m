function scalebar3D(xc,yc,zc, w, label)
%% scalebar3D.m
%
% Created by: Alessandro Ascani Orsini
%
% Date: 4/17/2024
%
% Version: 0.0.1
%
% *Description*: This function adds a 3D scalebar in the plot
hold on
x = [xc, xc-w, nan, xc, xc  , nan, xc, xc  ];
y = [yc, yc  , nan, yc, yc-w, nan, yc, yc  ];
z = [zc, zc  , nan, zc, zc  , nan, zc, zc+w];
hl = line(x,y,z,'HandleVisibility', 'off');
hl.LineWidth = 3;
hl.Color = [0,0,0];
ht = text(xc,yc,zc,[label]);
ht.FontSize = 18;
ht.Color = hl.Color;
ht.VerticalAlignment = 'bottom';
hold off
end