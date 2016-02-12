classdef mirroredNormalFisherBinghamDist < normalFisherBinghamDist
% this generalizes mirroredFisherBinghamDist

methods
	function obj = mirroredNormalFisherBinghamDist(varargin)
    obj = obj@normalFisherBinghamDist(varargin{:});
  end
	
	function val = pdf(obj,points)
    % points should be an nxd matrix
    
    for i = 2:length(obj.d)
      inds = sum(obj.d(1:i-1))+1:sum(obj.d(1:i));
      points(:,inds) = bsxfun(@times,normr(points(:,inds)),sign(points(:,inds)*-obj.V(inds,sum(obj.d(1:i)))));
    end
    
    val = pdf@normalFisherBinghamDist(obj,points);
  end
	
	function mnfb = marginal(obj,d)
    % d should be an nx1 array of indices
    %
    % We should check to make sure we're not trying to marginalize over a
    % subset of linked directional variables. Do this.
    
		warning('this only works if there are no Fisher components')
		
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
		n = 1000;
		t = linspace(0,2*pi,n+1)'; t = t(1:end-1);
		x = [cos(t),sin(t)];
		b = mirroredNormalFisherBinghamDist('d',[0,2],'B',obj.B(d,d),'a',[0;0]);
		p = b.pdf(x);
		m = x'*bsxfun(@times,x,p)*(t(2)-t(1));
		
    BMarg = obj.B(g,g)+2*obj.B(g,d)*m*obj.B(d,g);
		
    mnfb = mirroredNormalFisherBinghamDist('d',dMarg,'B',BMarg,'mu',modeMarg);
		
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
        m = m+1
        sample = candidate;
        if m>burnin
          x(m-burnin,:) = sample;
        end
        prevLike = logLike;
      else

      end
    end
  end
  
  function bingham = approximate(obj)
    x = obj.mode;
    
    % if mode is along the minor axis, enforce zero precision
    if all(abs(x*sign(x(1))-obj.V(:,1)*sign(obj.V(1,1)))<(2*eps))
      B = zeros(2);
    else
      B = (-obj.a*x' - x*obj.a' - obj.a'*x*(eye(2) - 3*(x*x'))...
        +2*obj.B - 4*obj.B*(x*x') - 4*(x*x')*obj.B - 2*x'*obj.B*x*(eye(2) - 4*(x*x')))/2;
    end
    
    bingham = mirroredNormalFisherBinghamDist('d',[0,2],'B',B,'a',[0;0]);
  end
		
end

methods (Static)
  function test()
%     Z = [-2,-2];
%     theta = pi/4;
%     V = [cos(theta),-sin(theta);sin(theta),cos(theta)];
%     mnfb = mirroredNormalFisherBinghamDist('d',[2,0],...
%       'V',V,...
%       'Z',Z,...
%       'mu',[0,0]);
%     plot(mnfb)
%     samples = mnfb.sample(100);
%     hold on, scatter(samples(:,1),samples(:,2)), hold off

%     Z = [-7, 0];
%     theta = -pi/4;
%     V = [cos(theta),-sin(theta);sin(theta),cos(theta)];
%     mnfb = mirroredNormalFisherBinghamDist('d',[0,2],...
%       'V',V,...
%       'Z',Z,...
%       'mu',[0,0]);
%     plot(mnfb)
%     samples = mnfb.sample(100);
%     hold on, scatter(samples(:,1),samples(:,2)), hold off
    
%     Z = [-2, -10, 0];
%     ry = -pi/4;
%     rz = pi/6; % rotates about z-axis
%     Ry = RotationMatrix.Rz(ry);
%     Rz = RotationMatrix.Rx(rz);
%     V = Rz*Ry;
%     mnfb = mirroredNormalFisherBinghamDist('d',[1,2],...
%       'V',V,...
%       'Z',Z,...
%       'mu',[10,0,0]);
%     figure(1), plot(mnfb)
%     samples = mnfb.sample(100);
%     hold on, scatter3(samples(:,2),samples(:,3),samples(:,1)), hold off
%     mnfbCond = mnfb.conditional(1,9);
%     figure(2), plot(mnfbCond)
%     mnfbAppx = mnfbCond.approximate;
%     figure(3), plot(mnfbAppx)
%     mnfbMarg = mnfb.marginal([2,3]);
%     figure(4), plot(mnfbMarg)
		
    Z = [-2, -2, -30, 0];
		Ry = RotationMatrix('euler',[2.1847 0.9425 3.6825]); %rand(1,3)); %
		Ry = blkdiag(Ry.R,1);
		Rz = RotationMatrix.Rx(pi/3);
		Rz = blkdiag(1,Rz);
		V = Rz*Ry;
    mnfb = mirroredNormalFisherBinghamDist('V',V,...
      'Z',Z,...
      'mu',[10,10,0,0],...
			'd',[2,2]);
    figure(1), plot(mnfb)
    samples = mnfb.sample(100);
    scale = max(max(samples(:,1))-min(samples(:,1)),max(samples(:,2))-min(samples(:,2)))/20;
    samples(:,3:4) = [samples(:,3).^2-samples(:,4).^2,2*samples(:,3).*samples(:,4)];
    hold on
    arrow(samples(:,1),samples(:,2),samples(:,3)*scale,samples(:,4)*scale,scale/3,'Color',0.9*ones(1,3))
    hold off
    mnfbCond = mnfb.conditional(1,10);
    figure(2), plot(mnfbCond)
    mnfbMarg = mnfb.marginal([3,4]);
    figure(3), plot(mnfbMarg)
  end
end

end