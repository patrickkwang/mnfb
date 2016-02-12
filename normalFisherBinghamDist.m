classdef normalFisherBinghamDist
% this generalizes BinghamDist, FisherBinghamDist, and normalBinghamDist
  
properties
  % dimension vector d = [n,m1,m2,...]
  % first element is the number of linear dimensions (it may be zero)
  d = [1,2];

	% orientation matrix
	V = eye(3);
  
  % Fisher term (encompasses the mean)
  a = [0,0,0];
	
  logC
end

properties (Dependent)
  B
  mu
	Z
end

properties (Hidden, SetAccess=private)
	% concentration matrix
	Z_ = [-1e1,-1e1,0];
end

methods
	function obj = normalFisherBinghamDist(varargin)
		obj = prtUtilAssignStringValuePairs(obj,varargin{:});
		obj = obj.sortEigenvectors();
		
		% compute normalization constant
		obj.logC = obj.logNormConst;
  end
	
  function val = logNormConst(obj)
    d = 1:obj.d(1);
    g = obj.d(1)+1:sum(obj.d);
    normalPart = obj.d(1)/2*log(2*pi)-1/2*log(det(-2*obj.B(d,d)));
    if sum(obj.d)>obj.d(1)
%       obj.C = obj.C*exp(logNormConstSP(sum(obj.d(2:end)),...
%         zeros(sum(obj.d(2:end)),1),...
%         obj.B(g,g)-obj.B(g,d)/obj.B(d,d)*obj.B(d,g)));
      b = obj.B(g,g)-obj.B(g,d)/obj.B(d,d)*obj.B(d,g);
      [~,z] = eig(b);
      binghamPart = log(bingham_F(z(1)));
      if any(obj.a)
        fisherPart = -1/4*obj.a'*pinv(obj.B)*obj.a;
      else
        fisherPart = 0;
      end
    else
      binghamPart = 0;
      fisherPart = 0;
    end
    val = normalPart+binghamPart+fisherPart;
  end
  
	function obj = sortEigenvectors(obj)
    for i=2:length(obj.d)
      inds = sum(obj.d(1:i-1))+1:sum(obj.d(1:i));
      [~,order] = sort(obj.Z_(inds));
      obj.V(:,inds) = obj.V(:,inds(order));
      obj.Z_(inds) = obj.Z_(inds(order));
			if obj.d(1)==0 % only if there are no linear components
				obj.Z_(inds) = obj.Z_(inds)-max(obj.Z_(inds));
			end
    end
	end
	
	function val = pdf(obj,points)
    % points should be an nxd matrix
    
    % compute exponent
    expt = points*obj.a;
    expt = expt + (points*obj.V).^2*obj.Z;
    
    % finish up
		val = exp(expt-obj.logC);
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
      j = find(d(i)<=cumsum(obj.d),1,'first');
      dCond(j) = dCond(j)-1;
    end

    % break into block matrix
    g = setdiff(1:sum(obj.d),d);
    aCond = obj.a(g) + 2*obj.B(g,d)*points';
    BCond = obj.B(g,g);
    
    if strcmp(class(obj),'normalFisherBinghamDist')
      nfb = normalFisherBinghamDist('d',dCond,'a',aCond,'B',BCond);
    elseif strcmp(class(obj),'mirroredNormalFisherBinghamDist')
      nfb = mirroredNormalFisherBinghamDist('d',dCond,'a',aCond,'B',BCond);
    end
  end

	function plot2(obj)
		% for distributions on the unit circle, plot theta vs p(theta)
		assert(obj.d(1)==0 && obj.d(2)==2, 'should be a Bingham distribution on S_1')
    
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
	
	function plot(obj)
    clf
		% generate sample points
    nThetas = 1000;
    thetas = linspace(0,2*pi,nThetas+1);
    thetas = thetas(1:nThetas);
    
    % the follow allows for plotting (Fisher-)Bingham distributions, too
    if obj.d(1)==2 && obj.d(2)==2
      n = 50;
      mnfbMarg = obj.marginal([3,4]);
      x = obj.mu(1)+1/sqrt(-mnfbMarg.B(1,1))*linspace(-3,3,n)';
      y = obj.mu(2)+1/sqrt(-mnfbMarg.B(2,2))*linspace(-3,3,n)';
      [X,Y] = meshgrid(x,y);
      samples = [X(:),Y(:)];
      v = mnfbMarg.pdf(samples);
      h = nan(n^2,1);
      s = nan(n^2,1);
      modes = nan(n^2,2);
      for i = 1:n^2
        mnfbCond = obj.conditional([1,2],samples(i,:));
        mnfbAppx = mnfbCond.approximate;
        modes(i,:) = mnfbAppx.V(:,2);
  %       figure(1), plot(mnfbAppx)
  %       figure(2), plot(mnfbCond)
        ang = atan2(mnfbAppx.V(2,2),mnfbAppx.V(1,2))*2;
        while ang>=2*pi
          ang = ang-2*pi;
        end
        while ang<0
          ang = ang+2*pi;
        end
    %     px = [ang/(2*pi), 1-exp(obj.Z(1))];
        h(i,:) = ang;
        s(i,:) = -mnfbAppx.Z(1);
        if mod(i,n)==0
          fprintf('plotting: %d of %d\n',i,size(samples,1));
        end
      end
      hsv = cat(2,h/(2*pi),min(s,quantile(s,0.9))/quantile(s,0.9),v/max(v));
      rgb = hsv2rgb(reshape(hsv,[n,n,3]));
      im = imagesc(x,y,rgb);
      xlabel('x')
      ylabel('y')
      
      % fix axis
      ax(1) = im.Parent;
      ax(1).Position = [ax(1).Position(1),ax(1).Position(2),ax(1).Position(3)-0.05,ax(1).Position(4)-0.1];
      axis equal tight xy
      
      % add hue colorbar
      pos = plotboxpos(im.Parent);
      ax(2) = axes('Position',[pos(1)+pos(3),pos(2),0.05,pos(4)]);
      imagesc(permute(colormap('hsv'),[1,3,2]))
      ax(2).YAxisLocation='right';
      ax(2).XTick=[];
      ax(2).YTick = pi/3*(0:6)/(2*pi)*64;
      ax(2).YTickLabel = {'0','\pi/3','2\pi/3','\pi','4\pi/3','5\pi/3','2\pi'};
      ylabel('angle (hue)')
      
      % add saturation colorbar
      ax(3) = axes('Position',[pos(1),pos(2)+pos(4),pos(3),0.05]);
      imagesc(hsv2rgb(bsxfun(@times,cat(3,ones(1,64),linspace(0,1,64),ones(1,64)),cat(3,linspace(0,1,64)',ones(64,1),ones(64,1)))))
      ax(3).XAxisLocation='top';
      ax(3).YTick=[];
      ax(3).XTick = linspace(0,64,6);
      ax(3).XTickLabel = cellfun(@(a){num2str(a,'%.2f')},num2cell(ax(3).XTick/64*quantile(s,0.9)));
      xlabel('precision (saturation)')
      
      % focus on main axis
      axes(ax(1))
%       downsamp = floor(n/12.5);
%       angles = atan2(modes(1:downsamp:end,2),modes(1:downsamp:end,1))*2;
%       figure, quiver(X(1:downsamp:end)',Y(1:downsamp:end)',cos(angles).*v(1:downsamp:end)/max(v),sin(angles).*v(1:downsamp:end)/max(v))
		elseif obj.d(1)==1 && obj.d(2)==2 % cylinder
      nZs = 100;
      zs = obj.mu(1)+1/sqrt(-obj.B(1,1))*linspace(-3,3,nZs); % assumes mean of 10
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
		elseif obj.d(1)==1 && obj.d(2)==0 % 1-D Gaussian
			% generate sample points
			x = -1/2/obj.B*obj.a+1/sqrt(-obj.B)*linspace(-3,3,1000)';

			% evaluate pdf at sample points
			p = obj.pdf(x);

			% plot
			plot(x,p)
			xlabel('x')
		elseif obj.d(1)==2 && obj.d(2)==0 % 2-D Gaussian
			% generate sample points
			x = obj.mu(1)+1/sqrt(-obj.B(1,1))*linspace(-3,3,1000)';
			y = obj.mu(2)+1/sqrt(-obj.B(2,2))*linspace(-3,3,1000)';
			[X,Y] = meshgrid(x,y);

			% evaluate pdf at sample points
			p = obj.pdf([X(:),Y(:)]);
			
			% display
			imagesc(x,y,reshape(p,size(X)))
      axis equal tight xy
      xlabel('x')
      ylabel('y')
      colormap gray
			
    elseif obj.d(1)==0 && obj.d(2)==2 % classic Bingham
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
  
  function val = get.mu(obj)
%     val = -1/2*(obj.B\obj.a);
    val = -1/2*pinv(obj.B)*obj.a;
  end
  
  function obj = set.a(obj,val)
    obj.a = val(:);
	end
	
	function obj = set.Z(obj,val)
		obj.Z_ = val(:);
		obj = obj.sortEigenvectors();
	end
	
	function val = get.Z(obj)
		val = obj.Z_;
	end
  
  function mode_hat = mode(obj)
		% find the major axis of the Gaussian
		m = obj.V(:,2);
		
		assert(m(1)~=0,'this does not work')
			
    % find the mean of the Gaussian
    mu_hat = (-2*obj.B*[1;obj.V(2,1)/obj.V(1,1)])\obj.a*[1;obj.V(2,1)/obj.V(1,1)];
    % to derive this, expand the expression a = -2 B mu
    % assuming that mu(2)==m2/m2*mu(1) (it's on the Gaussian's minor axis)

    % find the intersection of the line (y-mu2)/(x-mu1)=m(2)/m(1) with the
    % unit circle x^2+y^2=1 (y = +-sqrt(1-x^2)): http://mathworld.wolfram.com/Circle-LineIntersection.html
    x1 = mu_hat(1);
    x2 = mu_hat(1)+m(1);
    y1 = mu_hat(2);
    y2 = mu_hat(2)+m(2);
    dx = x2-x1;
    dy = y2-y1;
    dr = sqrt(dx^2+dy^2);
    D = x1*y2-x2*y1;

    mode_hat = {};
    if dr^2-D^2>=0 % at least one intersection
      mode_hat = cat(1,mode_hat,{[(D*dy+sign(dy)*dx*sqrt(dr^2-D^2))/dr^2;...
        (-D*dx+abs(dy)*sqrt(dr^2-D^2))/dr^2]});
    end
    if dr^2-D^2>0 % two intersections
      mode_hat = cat(1,mode_hat,{[(D*dy-sign(dy)*dx*sqrt(dr^2-D^2))/dr^2;...
        (-D*dx-abs(dy)*sqrt(dr^2-D^2))/dr^2]});
    end

    % and the minor axis
    m = obj.V(:,1);
    x2 = mu_hat(1)+m(1);
    y2 = mu_hat(2)+m(2);
    dx = x2-x1;
    dy = y2-y1;
    dr = sqrt(dx^2+dy^2);
    D = x1*y2-x2*y1;
    if dr^2-D^2>=0 % at least one intersection
      mode_hat = cat(1,mode_hat,{[(D*dy+sign(dy)*dx*sqrt(dr^2-D^2))/dr^2;...
        (-D*dx+abs(dy)*sqrt(dr^2-D^2))/dr^2]});
    end
    if dr^2-D^2>0 % two intersections
      mode_hat = cat(1,mode_hat,{[(D*dy-sign(dy)*dx*sqrt(dr^2-D^2))/dr^2;...
        (-D*dx-abs(dy)*sqrt(dr^2-D^2))/dr^2]});
    end
    
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