function plotUnitCylinderSolid(zlims,radius,color)
if nargin==0
	% test example
	zlims = [-1,1];
	radius = 1;
end

if nargin<3
	color = 0.5*[1,1,1];
end

theta = linspace(0,pi,200)';
x = cos(theta)*radius;
y = sin(theta)*radius;
o = ones(size(theta));

barsX = [x,x,nan(size(x))]; barsX = barsX(1:10:end,:)';
barsY = [y,y,nan(size(x))]; barsY = barsY(1:10:end,:)';
barsZ = [o*zlims(1),o*zlims(2),nan(size(x))]; barsZ = barsZ(1:10:end,:)';


X = cat(1,x,flipud(x));
Y = cat(1,y,flipud(y));
Z = cat(1,o*zlims(1),o*zlims(2));
h = fill3(X,Y,Z,...
	color,...
	X,-Y,Z,...
	color);
h(1).EdgeColor = 'None';
h(2).EdgeColor = 'None';
hold on
plot3(x,y,o*zlims(1),...
	x,y,o*zlims(2),...
	x,-y,o*zlims(1),...
	x,-y,o*zlims(2),...
	'-','Color',color)
hold off
set(gca,'ColorOrderIndex',1)