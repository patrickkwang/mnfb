function contour3d(X,Y,Z,V,varargin)
% contour3d(X,Y,Z,V,varargin)
%   plots a 3d contour with value Vi at position [Xi, Yi, Zi]
%
%   X, Y, Z, and V should be matrices of the same size

if nargin==0
	% test example
	disp('Testing mode...')
	[~,~,V] = peaks;
	[x,y] = meshgrid(1:size(V,2),1:size(V,1));
	X = cos(x/size(V,2)*pi);
	Y = sin(x/size(V,2)*pi);
	Z = y/size(V,1);
end

% default options
options.colormap = colormap('parula');
options.nLevels = 5;

% parse input options
options = prtUtilSimpleInputParser(options,varargin);

% construct contours from values
[x,y] = meshgrid(1:size(V,2),1:size(V,1));
C = contourc(1:size(V,2),1:size(V,1),V,options.nLevels);
[samples,levels] = deconstructContourMatrix(C);

% sort by bigness and punch out holes
[~,order] = sort(-abs(levels));
levels = levels(order);
samples = samples(order);
holysamples = breakOverlappingShapes(samples);

%% plot
% set up levels -> colors
levels = levels-min(levels);
levels = 1+ceil(levels/max(levels)*(size(options.colormap,1)-1));

% for each contour
for i = 1:length(holysamples)%-1
	% plot the outline
	xWarp = interp2(x,y,X,samples{i}(1,:),samples{i}(2,:));
	yWarp = interp2(x,y,Y,samples{i}(1,:),samples{i}(2,:));
	zWarp = interp2(x,y,Z,samples{i}(1,:),samples{i}(2,:));
	plot3(xWarp,yWarp,zWarp,'Color',options.colormap(levels(i),:))
	
	% fill in
	hold on
	xWarp = interp2(x,y,X,holysamples{i}(1,:),holysamples{i}(2,:));
	yWarp = interp2(x,y,Y,holysamples{i}(1,:),holysamples{i}(2,:));
	zWarp = interp2(x,y,Z,holysamples{i}(1,:),holysamples{i}(2,:));
	h = fill3(xWarp,yWarp,zWarp,options.colormap(levels(i),:));
	hold on
	h.EdgeColor = 'None';
end
hold off