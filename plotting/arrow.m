function varargout = arrow(x,y,u,v,l,c,varargin)
% plot an arrow starting at [x,y] with components [u,v] and head length l
if nargout>1
  error('Too many output arguments requested.')
end

x = x(:);
y = y(:);
u = u(:);
v = v(:);
l = l(:);
c = reshape(c,[],3);
theta = pi/6;
R = [cos(theta),sin(theta);-sin(theta),cos(theta)];
uva = normc(R*[u';v']); uva = [uva(1,:).*l';uva(2,:).*l'];
uvb = normc(R'*[u';v']); uvb = [uvb(1,:).*l';uvb(2,:).*l'];

% center the arrow
x = x-u/2;
y = y-v/2;

xToPlot = cat(2,x,x+u,x+u-uva(1,:)',x+u,x+u-uvb(1,:)')';
% xToPlot = xToPlot(:);
yToPlot = cat(2,y,y+v,y+v-uva(2,:)',y+v,y+v-uvb(2,:)')';
% yToPlot = yToPlot(:);
h = plot(xToPlot,yToPlot,varargin{:});
for i=1:length(x)
  h(i).Color = c(i,:);
end
if nargout==1
  varargout{1} = h;
end