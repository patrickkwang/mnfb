function plotUnitCylinder(zlims)
if nargin==0
	% test example
	zlims = [-1,1];
end

theta = linspace(0,2*pi,200)';
x = cos(theta);
y = sin(theta);
o = ones(size(theta));

barsX = [x,x,nan(size(x))]; barsX = barsX(1:10:end,:)';
barsY = [y,y,nan(size(x))]; barsY = barsY(1:10:end,:)';
barsZ = [o*zlims(1),o*zlims(2),nan(size(x))]; barsZ = barsZ(1:10:end,:)';

plot3(x,y,o*zlims(1),...
	x,y,o*zlims(2),...
	barsX(:),barsY(:),barsZ(:),...
	'--','Color',[0.8,0.8,0.8])