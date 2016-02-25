function plotOnCylinder(X,Y,Z,V,varargin)
% default options
options.radius = 1;
options.nLevels = 6;
options.camPos = [-2,-2,11.5];
options.plotCylinder = false;

% make gray-to-red colormap
cmap = colormap('lines');
red = cmap(2,:); % get the first red line
hsv = rgb2hsv(permute(red,[1,3,2]));
s = linspace(0,1,64);
h = repmat(hsv(1),1,64);
v = repmat(hsv(3),1,64);
rgb = hsv2rgb(cat(3,h,s,v));
options.cmap = permute(rgb,[2,3,1]);
% cmap = flipud(colormap('gray'));
% options.cmap = cmap(20:end,:);

% parse input options
options = utilSimpleInputParser(options,varargin);

% plot the contour
contour3d(X,Y,Z,V,...
  'nLevels',options.nLevels-1,...
  'colormap',options.cmap(ceil(size(cmap,1)/options.nLevels):end,:))

% plot the actual cylinder
if options.plotCylinder
	hold on
	plotUnitCylinderSolid([min(Z(:)),max(Z(:))],options.radius-0.015,options.cmap(1,:)) % 0.015
	hold off
end

% all the plotting options
axis tight
grid on
camproj('perspective')
campos(options.camPos)
alpha(0.5)
xlim([-1,1])
ylim([-1,1])

ax = gca;
ax.XTick = [-1,0,1];
ax.YTick = [-1,0,1];

%% fix order
dist = nan(length(ax.Children),1);
for i = 1:length(ax.Children)
  dist(i) = mean(sqrt((ax.Children(i).XData-options.camPos(1)).^2 ...
    +(ax.Children(i).YData-options.camPos(2)).^2 ...
    +(ax.Children(i).ZData-options.camPos(3)).^2));
end
[~,order] = sort(dist,1,'ascend');
ax.Children = ax.Children(order);