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
options = utilSimpleInputParser(options,varargin);

%% shift to start at minimum
% to avoid contours at the edge
[~,vInd] = min(sum(V,2));
[~,hInd] = min(sum(V,1));
V = circshift(V,-[vInd,hInd]);
X = circshift(X,-[vInd,hInd]);
Y = circshift(Y,-[vInd,hInd]);
Z = circshift(Z,-[vInd,hInd]);
topNbottom = mod([0,size(V,1)-1]-vInd,size(V,1))+1;
while any(topNbottom<1)
  topNbottom(topNbottom<1) = topNbottom(topNbottom<1)+size(V,1);
end
top = topNbottom(1);
bottom = topNbottom(2);

%% construct contours from values
[x,y] = meshgrid(1:size(V,2),1:size(V,1));
C = contourc(1:size(V,2),1:size(V,1),V,options.nLevels);
[samples,levels] = deconstructContourMatrix(C);

% sort by bigness
[~,order] = sort(-abs(levels));
levels = levels(order);
samples = samples(order);

% simplify polygons
for i = 1:length(samples)
  samples{i} = reduce_poly(samples{i},20);
end

% punch out holes
[holysamples,uncontained] = breakOverlappingShapes(samples);

% cylinder granularity parameter
spacing = 10;

%% wrap cylinder around contours
for i = 1:length(uncontained)
  lowBound(i) = max(min(samples{uncontained(i)}(1,:))-1,1);
  highBound(i) = min(max(samples{uncontained(i)}(1,:))+1,size(V,2));
  xd = [lowBound(i):spacing:highBound(i),highBound(i)];
  sample = [xd,fliplr(xd),xd(1),...
    samples{uncontained(i)}(1,[1:end,1]),xd(1),xd(1);
    bottom*ones(1,length(xd)),top*ones(1,length(xd)),bottom,...
    samples{uncontained(i)}(2,[1:end,1]),bottom,bottom];
  holysamples = cat(2,holysamples,sample);
  levels = cat(2,levels,0);
end

%% fill in the rest of the cylinder
[~,order] = sort(lowBound);
lowBound = lowBound(order);
highBound = highBound(order);
if lowBound(1)>1
  xd = [1:spacing:lowBound(1),lowBound(1)];
  for i = 2:length(xd)
    sample = [xd(i-1),xd(i),xd(i),xd(i-1);
      top,top,bottom,bottom];
    holysamples = cat(2,holysamples,sample);
    levels = cat(2,levels,0);
  end
end
for i = 2:length(lowBound)
  xd = [highBound(i-1):spacing:lowBound(i),lowBound(i)];
  for j = 2:length(xd)
    sample = [xd(j-1),xd(j),xd(j),xd(j-1);
      top,top,bottom,bottom];
    holysamples = cat(2,holysamples,sample);
    levels = cat(2,levels,0);
  end
end
if highBound(end)<size(V,2)
  xd = [highBound(end):spacing:size(V,2),size(V,2)];
  for i = 2:length(xd)
    sample = [xd(i-1),xd(i),xd(i),xd(i-1);
      top,top,bottom,bottom];
    holysamples = cat(2,holysamples,sample);
    levels = cat(2,levels,0);
  end
end
% connect the two edges
sample = [1,size(V,2),size(V,2),1;
  top,top,bottom,bottom];
holysamples = cat(2,holysamples,sample);
levels = cat(2,levels,0);

samples{length(holysamples)}=[];

%% plot
% set up levels -> colors
levels = levels-min(levels);
levels = 1+ceil(levels/max(levels)*(size(options.colormap,1)-1));

% for each contour
for i = 1:length(holysamples)
  if ~isempty(samples{i})
    % plot the outline
    xWarp = interp2(x,y,X,samples{i}(1,:),samples{i}(2,:));
    yWarp = interp2(x,y,Y,samples{i}(1,:),samples{i}(2,:));
    zWarp = interp2(x,y,Z,samples{i}(1,:),samples{i}(2,:));
    plot3([xWarp,xWarp(1)],...
      [yWarp,yWarp(1)],...
      [zWarp,zWarp(1)],...
      'Color',options.colormap(levels(i),:))
  end
	
	% fill in
	hold on
	xWarp = interp2(x,y,X,holysamples{i}(1,:),holysamples{i}(2,:));
	yWarp = interp2(x,y,Y,holysamples{i}(1,:),holysamples{i}(2,:));
	zWarp = interp2(x,y,Z,holysamples{i}(1,:),holysamples{i}(2,:));
	h = fill3(xWarp,yWarp,zWarp,options.colormap(levels(i),:));
	hold on
	h.LineStyle = 'none';
end
hold off