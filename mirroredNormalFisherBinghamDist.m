classdef mirroredNormalFisherBinghamDist < normalFisherBinghamDist
% this generalizes mirroredFisherBinghamDist

methods
	function obj = mirroredNormalFisherBinghamDist(varargin)
		obj = utilAssignStringValuePairs(obj,varargin{:});
		obj = obj.sortEigenvectors();
		
		% compute normalization constant
% 		obj.C = exp(logNormConstSP(2,obj.a,obj.B));
    obj.C = 1;
  end
	
	function val = pdf(obj,points)
    % points should be an nxd matrix
    
    for i = 2:length(obj.d)
      inds = sum(obj.d(1:i-1))+1:sum(obj.d(1:i));
      points(:,inds) = bsxfun(@times,normr(points(:,inds)),sign(points(:,inds)*-obj.V(inds,sum(obj.d(1:i)))));
    end
    
    % compute exponent
    expt = obj.a'*points';
    for i = 1:sum(obj.d)
      expt = expt + obj.Z(i)*(obj.V(:,i)'*points').^2;
    end
    
    % finish up
		val = 1/obj.C*exp(expt)';
  end
  
  function mnfb = conditional(obj,d,points)
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
    
    mnfb = utilAssignStringValuePairs(obj,'d',dCond,'a',aCond,'B',BCond);
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
		b = normalFisherBinghamDist('d',[0,2],'B',obj.B(d,d),'a',[0;0]);
		p = b.pdf(x);
		m = x'*bsxfun(@times,x,p)*(t(2)-t(1));
		
    BMarg = obj.B(g,g)+2*obj.B(g,d)*m*obj.B(d,g);
		
    mnfb = utilAssignStringValuePairs(obj,'d',dMarg,'B',BMarg,'mu',modeMarg);
		
	end
  
  function bingham = approximate(obj)
    x = obj.mode;
    B = (-obj.a*x' - x*obj.a' - obj.a'*x*(eye(2) - 3*(x*x'))...
      +2*obj.B - 4*obj.B*(x*x') - 4*(x*x')*obj.B - 2*x'*obj.B*x*(eye(2) - 4*(x*x')))/2;
    [V,Z] = eig(B);
    bingham = utilAssignStringValuePairs(obj,'d',[0,2],'B',B,'a',[0;0]);
  end
  
  function px = hs(obj)
    % this is designed for a normal-Bingham distribution (d=[1,2])
    % infer the angle from the eigenvectors (hue)
    % take the first eigenvalue to represent the concentration (saturation)
    
		ang = atan2(obj.V(2,2),obj.V(1,2))*2;
		while ang>=2*pi
			ang = ang-2*pi;
		end
		while ang<0
			ang = ang+2*pi;
		end
    px(1) = ang/(2*pi);
    px(2) = 1-exp(obj.Z(1));
    
  end
		
end

methods (Static)
  function test()
%     Z = [-2, -10, 0];
%     ry = -pi/4;
%     rz = pi/6; % rotates about z-axis
%     Ry = RotationMatrix.Rz(ry);
%     Rz = RotationMatrix.Rx(rz);
%     V = Rz*Ry;
%     mnfb = mirroredNormalFisherBinghamDist('V',V,...
%       'Z',Z,...
%       'mu',[10,0,0]);
%     figure(1), plot3(mnfb)
%     mnfbCond = mnfb.conditional(1,9);
%     figure(2), plot3(mnfbCond)
%     mnfbAppx = mnfbCond.approximate;
%     figure(3), plot3(mnfbAppx)
%     mnfbMarg = mnfb.marginal([2,3]);
%     figure(4), plot3(mnfbMarg)
%     % plot with colors
% 		mu = -1/2*(mnfbMarg.B\mnfbMarg.a);
% 		zs = mu+1/sqrt(-mnfbMarg.B)*linspace(-3,3,100)';
% 		v = mnfbMarg.pdf(zs);
% 		for i = 1:100
% 			mnfbCond = mnfb.conditional(1,zs(i));
% 			mnfbAppx = mnfbCond.approximate;
% 			hs(i,:) = mnfbAppx.hs;
% 		end
% 		hsv = cat(2,hs,v/max(v));
% 		rgb = hsv2rgb(permute(hsv,[1,3,2]));
% 		figure(5), imagesc(rgb)
		
    Z = [-2, -2, -10, 0];
		Ry = RotationMatrix('euler',[2.1847 0.9425 3.6825]); %rand(1,3)
		Ry = blkdiag(Ry.R,1);
		Rz = RotationMatrix.Rx(pi/3);
		Rz = blkdiag(1,Rz);
		V = Rz*Ry;
    mnfb = mirroredNormalFisherBinghamDist('V',V,...
      'Z',Z,...
      'mu',[10,10,0,0],...
			'd',[2,2]);
    mnfbCond = mnfb.conditional(1,10);
    figure(2), plot3(mnfbCond)
    mnfbMarg = mnfb.marginal([3,4]);
    figure(4), plot3(mnfbMarg)
    % plot with colors
		mu = -1/2*(mnfbMarg.B\mnfbMarg.a);
		x = mu(1)+1/sqrt(-mnfbMarg.B(1,1))*linspace(3,-3,100)';
		y = mu(2)+1/sqrt(-mnfbMarg.B(2,2))*linspace(-3,3,100)';
		[X,Y] = meshgrid(x,y);
		samples = [X(:),Y(:)];
		v = mnfbMarg.pdf(samples);
		hs = nan(size(samples,1),2);
		for i = 1:size(samples,1)
			mnfbCond = mnfb.conditional([1,2],samples(i,:));
			mnfbAppx = mnfbCond.approximate;
			hs(i,:) = mnfbAppx.hs;
			if mod(i,10)==0
				fprintf('%d of %d\n',i,size(samples,1));
			end
		end
		hsv = cat(2,hs,v/max(v));
		rgb = hsv2rgb(reshape(hsv,[100,100,3]));
		figure(5), imagesc(rgb)
  end
end

end