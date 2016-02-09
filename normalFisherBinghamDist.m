classdef normalFisherBinghamDist
% this generalizes BinghamDist, FisherBinghamDist, and normalBinghamDist
  
properties
  % dimension vector d = [n,m1,m2,...]
  % first element is the number of linear dimensions (it may be zero)
  d = [1,2];
  
	% concentration matrix
	Z = [-1e1,-1e1,0];

	% orientation matrix
	V = eye(3);
  
  % Fisher term (encompasses the mean)
  a = [0,0,0];
	
	C
end

properties (Dependent)
  B
  mu
end  

methods
	function obj = normalFisherBinghamDist(varargin)
		obj = prtUtilAssignStringValuePairs(obj,varargin{:});
		obj = obj.sortEigenvectors();
		
		% compute normalization constant
% 		obj.C = exp(logNormConstSP(2,obj.a,obj.B));
    obj.C = 1;
	end
	
	function obj = sortEigenvectors(obj)
    for i=2:length(obj.d)
      inds = sum(obj.d(1:i-1))+1:sum(obj.d(1:i));
      [~,order] = sort(obj.Z(inds));
      obj.V(:,inds) = obj.V(:,inds(order));
      obj.Z(inds) = obj.Z(inds(order))-max(obj.Z(inds));
    end
	end
	
	function val = pdf(obj,points)
    % points should be an nxd matrix
    
    % compute exponent
    expt = obj.a'*points';
    for i = 1:sum(obj.d)
      expt = expt + obj.Z(i)*(obj.V(:,i)'*points').^2;
    end
    
    % finish up
		val = 1/obj.C*exp(expt);
  end
  
  function nfb = conditional(obj,d,points)
    % d should be an nx1 array of indices
    % points should be nx1 array of points
    %
    % We should check to make sure we're not trying to marginalize over a
    % subset of linked directional variables. Do this.
    
    % fix variable counts
    dCond = obj.d;
    for i = 1:length(d)
      j = find(d<cumsum(obj.d))-1;
      dCond(j) = dCond(j)-1;
    end

    % break into block matrix
    g = setdiff(1:sum(obj.d),d);
    aCond = obj.a(g) + 2*obj.B(g,d)*points';
    BCond = obj.B(g,g);
    
    nfb = normalFisherBinghamDist('a',aCond,'B',BCond,'d',dCond);
  end

	function plot2(obj)
		% for distributions on the unit circle, plot theta vs p(theta)
		
		% generate sample points
		n = 10000;
		t = linspace(0,2*pi,n+1)'; % t = t(1:end-1);
		x = [-sin(t),cos(t)];
		
		% evaluate pdf at sample points
		p = obj.pdf(x);
		
		% plot
		h = scatter(x(:,1),x(:,2),6,p);
		h.Parent.XTick = -1:1;
		h.Parent.YTick = -1:1;
		grid on
		axis equal tight
		
		% labels
		xlabel('x')
		ylabel('y')
	end
	
	function plot3(obj)
		% generate sample points
    nThetas = 1000;
    thetas = linspace(0,2*pi,nThetas+1);
    thetas = thetas(1:nThetas);
    
    % the follow allows for plotting (Fisher-)Bingham distributions, too
    if obj.d(1)
      nZs = 100;
      zs = 10+linspace(-3,3,nZs)/sqrt(-2*obj.Z(1)); % assumes mean of 10
      [T,X] = meshgrid(thetas,zs);
      x = [X(:),cos(T(:)),sin(T(:))];
    
      % compute likelihood
      like = obj.pdf(x);
      like = reshape(like,nZs,nThetas);

      % plot on cylinder
      plotOnCylinder(cos(T),sin(T),X,like)
      
      xlabel('q_1')
      ylabel('q_2')
      zlabel('x')
    else
      x = [cos(thetas(:)),sin(thetas(:))];
      
      % evaluate pdf at sample points
      p = obj.pdf(x);

      % plot
      plot3(x(:,1),x(:,2),p)

      % side bars
      hold on
      plot3(x(:,1),x(:,2),zeros(size(p)));
      stem3(x(round(1:nThetas/90:end),1),...
        x(round(1:nThetas/90:end),2),...
        p(round(1:nThetas/90:end)),'Marker','none')
      hold off

      % labels
      set(gca,'XTick',-0.9:0.3:0.9)
      set(gca,'YTick',-0.9:0.3:0.9)
      grid on
      xlabel('q_1')
      ylabel('q_2')
      zlabel('p([x,y]|a,B)')
    end
  end
  
  function val = get.B(obj)
    val = obj.V*diag(obj.Z)*obj.V';
  end
  
  function obj = set.B(obj,val)
    [v,z] = eig(val);
    obj.V = v;
    obj.Z = diag(z);
  end
  
  function obj = set.mu(obj,val)
    obj.a = -2*obj.B*val(:);
  end
  
  function obj = set.a(obj,val)
    obj.a = val(:);
  end
  
  function mode_hat = mode(obj)
    % find the mean of the Gaussian
    % we can only know that the mean lies somewhere on the major axis
    % let's say that it's on the x-axis...
    mu_hat = [-1/2*obj.a(1)/obj.V(1,1)^2/obj.Z(1);0];
    % to derive this, expand the expression a = -2 B mu
    % assuming that mu(2)==0
    
    % find the major axis of the Gaussian
    m = -obj.V(:,2);
    
    % find the intersection of the line (y-mu2)/(x-mu1)=m(2)/m(1) with the
    % unit circle x^2+y^2=1 (y = +-sqrt(1-x^2)):
    % 
    % here are the coefficients of the resulting quadratic equation
    a_ = 1+(m(2)/m(1)).^2;
    b_ = -2*(m(2)/m(1)).^2*mu_hat(1);
    c_ = (m(2)/m(1)).^2*mu_hat(1).^2-1;
    
    % find the four solutions to this quadratic equation
    mode_hatX = (-b_+sqrt(b_^2-4*a_*c_))/(2*a_);
    mode_hatY = sqrt(1-mode_hatX^2);
    mode_hat{1} = [mode_hatX;mode_hatY];
    mode_hatY = -sqrt(1-mode_hatX^2);
    mode_hat{3} = [mode_hatX;mode_hatY];
    mode_hatX = (-b_-sqrt(b_^2-4*a_*c_))/(2*a_);
    mode_hatY = sqrt(1-mode_hatX^2);
    mode_hat{2} = [mode_hatX;mode_hatY];
    mode_hatY = -sqrt(1-mode_hatX^2);
    mode_hat{4} = [mode_hatX;mode_hatY];
    
    % pick the intersection with the greatest likelihood
    llhs = obj.pdf(cat(2,mode_hat{:})');
    [~,ind] = max(llhs);
    mode_hat = mode_hat{ind};
  end
		
end

methods (Static)
  function test()
    Z = [-2, -10, 0];
    ry = -pi/4;
    rz = pi/6; % rotates about z-axis
    Ry = RotationMatrix.Rz(ry);
    Rz = RotationMatrix.Rx(rz);
    V = Rz*Ry;
    nfb = normalFisherBinghamDist('V',V,...
      'Z',Z,...
      'mu',[10,0,0]);
    figure(1), plot3(nfb)
    nfbCond = nfb.conditional(1,10.5);
    figure(2), plot3(nfbCond)
    nfbCond.mode
  end
end

end