function varargout = arrow(x,y,u,v,l,varargin)
% plot an arrow starting at [x,y] with components [u,v] and head length l

x = x(:);
y = y(:);
u = u(:);
v = v(:);
theta = pi/6;
R = [cos(theta),sin(theta);-sin(theta),cos(theta)];
uva = R*[u';v']; uva = bsxfun(@rdivide,uva,sqrt(sum(uva.^2,1)))*l;
uvb = R'*[u';v']; uvb = bsxfun(@rdivide,uvb,sqrt(sum(uvb.^2,1)))*l;
xToPlot = [x,x+u,nan(size(x)),...
	x+u,x+u-uva(1,:)',nan(size(x)),...
	x+u,x+u-uvb(1,:)',nan(size(x))]';
xToPlot = xToPlot(:);
yToPlot = [y,y+v,nan(size(y)),...
	y+v,y+v-uva(2,:)',nan(size(y)),...
	y+v,y+v-uvb(2,:)',nan(size(y))]';
yToPlot = yToPlot(:);
if nargout==1
  varargout{1} = plot(xToPlot,yToPlot,varargin{:});
elseif nargout>1
  error('Too many output arguments requested.')
else
  plot(xToPlot,yToPlot,varargin{:})
end