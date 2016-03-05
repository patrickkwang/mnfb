classdef normalFisherBinghamDist
% this generalizes the Bingham, Fisher-Bingham, and normal-Bingham
% distributions
  
properties
  % dimension vector d = [n,m1,m2,...]
  % first element is the number of linear dimensions (it may be zero)
  d = [1,2];

	% orientation matrix
	V = eye(3);
  
  % Fisher term (encompasses the mean)
  a = [0,0,0];
	
  %  log normalizing constant
  logC
end

properties (Dependent)
  B
  mu
	Z
  mode
end

properties (Hidden, SetAccess=private)
	% concentration matrix
	Z_ = [-1e1,-1e1,0];
end

methods
	function obj = normalFisherBinghamDist(varargin)
		obj = utilAssignStringValuePairs(obj,varargin{:});
		obj = obj.sortEigenvectors();
		
    if isempty(obj.logC)
      % compute normalization constant
      obj.logC = obj.logNormConst;
    end
    
    if nargout==0
      clear obj
      normalFisherBinghamDist.test
    end
  end
	
  function val = logNormConst(obj)
    d = 1:obj.d(1);
    g = obj.d(1)+1:sum(obj.d);
    normalPart = obj.d(1)/2*log(2*pi)-1/2*log(det(-2*obj.B(d,d)));
    if sum(obj.d)>obj.d(1)
      b = obj.B(g,g)-obj.B(g,d)/obj.B(d,d)*obj.B(d,g);
%       [~,z] = eig(b);
%       binghamPart = log(bingham_F(z(1))); % only works for a subset of parameters
      binghamPart = logNormConstSP(sum(obj.d(2:end)),...
        zeros(sum(obj.d(2:end)),1),...
        b);
    else
      binghamPart = 0;
    end
    if any(obj.a)
      fisherPart = -1/4*obj.a'*pinv(obj.B)*obj.a;
    else
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
    val = exp(obj.logPdf(points));
  end
  
  function val = logPdf(obj,points)
    % points should be an nxd matrix
    
    % compute exponent
    expt = points*obj.a + (points*obj.V).^2*obj.Z;
    
    % finish up
		val = expt-obj.logC;
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
    
    if isa(obj,'mirroredNormalFisherBinghamDist')
      nfb = mirroredNormalFisherBinghamDist('d',dCond,'a',aCond,'B',BCond,'logC',1);
    else
      nfb = normalFisherBinghamDist('d',dCond,'a',aCond,'B',BCond,'logC',1);
    end
  end
	
	function mnfb = marginal(obj,d)
    % d should be an nx1 array of indices
    %
    % We should check to make sure we're not trying to marginalize over a
    % subset of linked directional variables. Do this.
    
    % this is a bit of a hack
    % - these components are very small but not below eps
		if any(abs(obj.mu(obj.d(1)+1:end)/obj.mu(1))>10*eps)
      warning('There should be no Fisher components.')
    end
		
    % fix variable counts
    dMarg = obj.d;
    for i = 1:length(d)
      j = find(d(i)<=cumsum(obj.d),1,'first');
      dMarg(j) = dMarg(j)-1;
    end
		
    % break into block matrix
		d = sort(d);
    g = setdiff(1:sum(obj.d),d);
		mode = [-1/2*(obj.B(1:obj.d(1),1:obj.d(1))\obj.a(1:obj.d(1)));zeros(sum(obj.d(2:end)),1)];
		modeMarg = mode(g);
		
		% compute second moment
		n = 180;
    
    % this has been standard
% 		t = linspace(0,2*pi,n+1)'; t = t(1:end-1);
%     x = [cos(t),sin(t)];

    % this may work better for high precisions (and worse for low)
    [V,Z] = eig(obj.B(d,d));
    z = diag(Z);
    [~,minInd] = min(z);
    [~,maxInd] = max(z);
    t = linspace(-4,4,n)'*sqrt(-2/z(minInd))+atan2(V(2,maxInd),V(1,maxInd));
		x = [cos(t),sin(t)]; x = cat(1,x,-x);
        
    if isa(obj,'mirroredNormalFisherBinghamDist')
      b = mirroredNormalFisherBinghamDist('d',[0,2],'B',obj.B(d,d),'a',[0;0]);
    else
      b = normalFisherBinghamDist('d',[0,2],'B',obj.B(d,d),'a',[0;0]);
    end
		lp = b.logPdf(x); p = exp(lp); %exp(lp-prtUtilSumExp(lp));
		m = x'*bsxfun(@times,x,p)*(t(2)-t(1))/(2*trapz(t,p(1:180))); % with extra normalization
    
    BMarg = obj.B(g,g)+2*obj.B(g,d)*m*obj.B(d,g);
		
    if isa(obj,'mirroredNormalFisherBinghamDist')
      mnfb = mirroredNormalFisherBinghamDist('d',dMarg,'B',BMarg,'mu',modeMarg);
    else
      mnfb = normalFisherBinghamDist('d',dMarg,'B',BMarg,'mu',modeMarg);
    end
  end
	
	function plot(obj)
		% generate sample points
    nThetas = 180;
    thetas = linspace(0,2*pi,nThetas+1)';
    thetas = thetas(1:nThetas);
    
    % the follow allows for plotting (Fisher-)Bingham distributions, too
    if obj.d(1)==2 && obj.d(2)==2
      n = 50;
      mnfbMarg = obj.marginal([3,4]);
      x = obj.mu(1)+1/sqrt(-mnfbMarg.B(1,1))*linspace(-3,3,n)';
      y = obj.mu(2)+1/sqrt(-mnfbMarg.B(2,2))*linspace(-3,3,n)';
      [X,Y] = meshgrid(x,y);
      samples = [X(:),Y(:)];
      v = reshape(mnfbMarg.pdf(samples),size(X));
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
      hsv = cat(2,h/(2*pi),min(s,quantile(s,0.9))/quantile(s,0.9),v(:)./max(v(:)));
      rgb = hsv2rgb(reshape(hsv,[n,n,3]));
      h = imagesc(x,y,rgb);
%       downsamp = 2;
%       arrowLength = (Y(2)-Y(1))*downsamp;
%       angles = reshape(atan2(modes(:,2),modes(:,1))*2,size(X));
%       h = arrow(X(1:downsamp:end,1:downsamp:end)',...
%         Y(1:downsamp:end,1:downsamp:end)',...
%         arrowLength*v(1:downsamp:end,1:downsamp:end)/max(v(:)).*cos(angles(1:downsamp:end,1:downsamp:end)),...
%         arrowLength*v(1:downsamp:end,1:downsamp:end)/max(v(:)).*sin(angles(1:downsamp:end,1:downsamp:end)),...
%         arrowLength/4*v(1:downsamp:end,1:downsamp:end)/max(v(:)),...
%         rgb(1:downsamp:end,1:downsamp:end,:));
%       h(1).Parent.Color = [0,0,0];
      
      xlabel('$x$')
      ylabel('$y$')
      
      % fix axis
      ax(1) = h.Parent;
      ax(1).Position = [ax(1).Position(1),ax(1).Position(2),ax(1).Position(3)-0.05,ax(1).Position(4)-0.1];
      axis xy% equal tight
      
      ax(1).Units = 'pixels';
      pos = plotboxpos(ax(1));
      
      % add hue colorbar
      ax(2) = axes('Units','pixels','Position',[pos(1)+pos(3),pos(2),30,pos(4)]);
      imagesc([],linspace(0,2*pi,64),permute(colormap('hsv'),[1,3,2]))
      ax(2).YAxisLocation='right';
      ax(2).XTick=[];
      tickLabelsToPiFractions(ax(2),'y',3)
%       ax(2).YTick = pi/3*(0:6)/(2*pi)*64;
%       ax(2).YTickLabel = {'0','\pi/3','2\pi/3','\pi','4\pi/3','5\pi/3','2\pi'};
      ylabel('angle (hue)')
      
      % add saturation colorbar
      ax(3) = axes('Units','pixels','Position',[pos(1),pos(2)+pos(4),pos(3),30]);
      imagesc(hsv2rgb(bsxfun(@times,cat(3,ones(1,64),linspace(0,1,64),ones(1,64)),cat(3,linspace(0,1,64)',ones(64,1),ones(64,1)))))
      ax(3).XAxisLocation='top';
      ax(3).YTick=[];
      ax(3).XTick = linspace(0,64,6);
      ax(3).XTickLabel = cellfun(@(a){num2str(a,'%.2f')},num2cell(ax(3).XTick/64*quantile(s,0.9)));
      xlabel('precision (saturation)')
      
      % focus on main axis
      axes(ax(1))
		elseif obj.d(1)==1 && obj.d(2)==2 % cylinder
      nZs = 100;
      zs = obj.mu(1)+1/sqrt(-obj.B(1,1))*linspace(-3,3,nZs); % assumes mean of 10
      [T,X] = meshgrid(thetas,zs);
      x = [X(:),cos(T(:)),sin(T(:))];
    
      % compute likelihood
      like = obj.pdf(x);
      like = reshape(like,nZs,nThetas);

      % plot on cylinder
      zStd = 1/sqrt(-obj.B(1,1));
      options = struct('camPos',obj.mode([2,3,1]).*[3,3,1]'+[0,0,8*zStd]',...
        'cmap',cat(2,linspace(0.8,0,64)',...
          linspace(0.8,0,64)',...
          linspace(0.8,0.8,64)'));
      plotOnCylinder(cos(T),sin(T),X,like,options)
      
      xlabel('$q_1$')
      ylabel('$q_2$')
      zlabel('$x$')
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
      axis xy% equal tight
      xlabel('$x$')
      ylabel('$y$')
      colormap('gray')
			
    elseif obj.d(1)==0 && obj.d(2)==2 % classic Bingham
      x = [cos(thetas([1:end,1])),sin(thetas([1:end,1]))];
      
      % evaluate pdf at sample points
      p = obj.pdf(x);

      % plot
      h = plot3([x(round(1:nThetas/90:end),1),x(round(1:nThetas/90:end),1)]',...
        [x(round(1:nThetas/90:end),2),x(round(1:nThetas/90:end),2)]',...
        [zeros(size(p(round(1:nThetas/90:end)))),p(round(1:nThetas/90:end))]',...
        x(:,1),x(:,2),zeros(size(p)),...
        x(:,1),x(:,2),p,...
        'Marker','none');
      colors = get(groot,'DefaultAxesColorOrder');
      for i = 1:length(h)-2
        h(i).Color = colors(3,:);
      end
      h(end-1).Color = colors(2,:);
      h(end).Color = colors(1,:);

      % perspective
      littleSpin = -pi/6;
      R = [cos(littleSpin), sin(littleSpin); -sin(littleSpin), cos(littleSpin)];
      h(1).Parent.CameraPosition = 3*[R*obj.mode*2;0.5]';
      h(1).Parent.Projection = 'perspective';
      % axis equal
      h(1).Parent.DataAspectRatio(2) = h(1).Parent.DataAspectRatio(1);
      xlim([-1.1,1.1])
      ylim([-1.1,1.1])
      
      % labels
      set(gca,'XTick',-1:0.5:1)
      set(gca,'YTick',-1:0.5:1)
      grid on
      xlabel('$q_1$')
      ylabel('$q_2$')
      zlabel('$p(q|M,Z)$')
    end
  end

	function plot1(obj,varargin)
		% for distributions on the unit circle, plot theta vs p(theta)
		assert(obj.d(1)==0 && obj.d(2)==2, 'should be a Bingham distribution on S_1')
    
		% generate sample points
		n = 1000;
		t = linspace(0,2*pi,n+1)'; % t = t(1:end-1);
		x = [-sin(t),cos(t)];
		
		% evaluate pdf at sample points
		p = obj.pdf(x);
		
		% plot
		plot(2*t,p,varargin{:})
    xlim([0,4*pi])
    tickLabelsToPiFractions(gca,'x',2)
		
		% labels
		xlabel('$\theta$')
		ylabel('$p(\theta|M,Z)$')
  end

	function plot2(obj)
		% for distributions on the unit circle, plot theta vs p(theta)
		assert(obj.d(1)==0 && obj.d(2)==2, 'should be a Bingham distribution on S_1')
    
    x = linspace(-1.5,1.5,100);
    y = linspace(-1.5,1.5,100);
    [X,Y] = meshgrid(x,y);
    im = exp(([X(:),Y(:)]*obj.V).^2*obj.Z);
    h = imagesc(x,y,reshape(im,size(X)));
    h.AlphaData = 0.2;
    
		% generate sample points
		n = 250;
		t = linspace(0,2*pi,n+1)'; % t = t(1:end-1);
		x = [-sin(t),cos(t)];
		
		% evaluate pdf at sample points
		p = exp((x*obj.V).^2*obj.Z);
    parula = colormap('parula');
    c = parula(ceil(p./max(p)*64),:);
		
		% plot
    hold on
		h = scatter(x(:,1),x(:,2),6,c,...
      'MarkerFaceColor','flat');
    hold off
		h.Parent.XTick = -1:1;
		h.Parent.YTick = -1:1;
		grid on
		axis equal tight
		
		% labels
		xlabel('$x$')
		ylabel('$y$')
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
  
  function val = get.mode(obj)
    val = nan(sum(obj.d),1);
    val(1:obj.d(1)) = obj.mu(1:obj.d(1));
    for i = 2:length(obj.d)
      inds = obj.d(i-1)+(1:obj.d(i))';
      assert(length(inds)==2,'Mode currently works only for 1D rotations.')

      % find the mean of the Gaussian
      % to derive this, expand the expression a = -2 B mu
      % assuming that mu(2)==m2/m1*mu(1) (it's on the Gaussian's minor axis)
      a_demeaned = -2*obj.B*(obj.mu.*[zeros(obj.d(1),1);ones(sum(obj.d(2:end)),1)]);
      mu_hat = (-2*obj.B(inds,inds)*[1;obj.V(inds(2),inds(1))/obj.V(inds(1),inds(1))])\a_demeaned(inds)*[1;obj.V(inds(2),inds(1))/obj.V(inds(1),inds(1))];

      mode_hat = {};
      x1 = mu_hat(1);
      y1 = mu_hat(2);

      % find the intersection of the line (y-mu2)/(x-mu1)=m(2)/m(1) with the
      % unit circle x^2+y^2=1:
      % http://mathworld.wolfram.com/Circle-LineIntersection.html
      % first with the major axis of the Gaussian
      m = obj.V(inds,inds(2));
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

      % then the minor axis
      m = obj.V(inds,inds(1));
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

      % force real modes
      mode_hat = cellfun(@(a){real(a)},mode_hat);
      
      % pick the intersection with the greatest likelihood
      llhs = obj.pdf(cat(2,repmat(val(1:obj.d(1))',length(mode_hat),1),cat(2,mode_hat{:})'));
      [~,ind] = max(llhs);
      val(inds) = mode_hat{ind};
    end
  end
  
  function x = sample(obj,n,burnin)
    if nargin<3
      burnin = 100;
    end
    
    addnoise = @(x,stdv)x+randn(size(x)).*stdv(:);
    prop = @(x)addnoise(x(1:obj.d(1)),1./sqrt(-1/2*obj.Z(1:obj.d(1))));
    for i = 2:length(obj.d)
      prop = @(x)cat(1,prop(x),normc(addnoise(x(sum(obj.d(1:i-1))+(1:obj.d(i))),3)));
    end

    m = 0;
    x = zeros(n,sum(obj.d));

    sample = obj.mu;
    prevLike = -Inf;
    while m<n+burnin
      candidate = prop(sample);
      %logLike = sum(bpLikelihood(candidate,mode,Y));
      logLike = log(obj.pdf(candidate'));
      r = logLike-prevLike;
      if r>0 || exp(r)>rand
        m = m+1;
        sample = candidate;
        if m>burnin
          x(m-burnin,:) = sample;
        end
        prevLike = logLike;
      else

      end
    end
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
    figure(1), plot(nfb), axis equal
    nfbCond = nfb.conditional(1,10.5);
    figure(2), plot(nfbCond)
  end
end

end