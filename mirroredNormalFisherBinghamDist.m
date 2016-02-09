classdef mirroredNormalFisherBinghamDist < normalFisherBinghamDist
% this generalizes mirroredFisherBinghamDist

methods
	function obj = mirroredNormalFisherBinghamDist(varargin)
		obj = prtUtilAssignStringValuePairs(obj,varargin{:});
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
		val = 1/obj.C*exp(expt);
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
      j = find(d<cumsum(obj.d))-1;
      dCond(j) = dCond(j)-1;
    end

    % break into block matrix
    g = setdiff(1:sum(obj.d),d);
    aCond = obj.a(g) + 2*obj.B(g,d)*points';
    BCond = obj.B(g,g);
    
    mnfb = mirroredNormalFisherBinghamDist('a',aCond,'B',BCond,'d',dCond);
  end
  
  function bingham = approximate(obj)
    x = obj.mode;
    B = (-obj.a*x' - x*obj.a' - obj.a'*x*(eye(2) - 3*(x*x'))...
      +2*obj.B - 4*obj.B*(x*x') - 4*(x*x')*obj.B - 2*x'*obj.B*x*(eye(2) - 4*(x*x')))/2;
    [V,Z] = eig(B);
    bingham = mirroredNormalFisherBinghamDist('Z',diag(Z),'V',V,'d',[0,2],'a',[0;0]);
  end
  
  function px = hs(obj)
    % this is designed for a normal-Bingham distribution (d=[1,2])
    % infer the angle from the eigenvectors (hue)
    % take the first eigenvalue to represent the concentration (saturation)
    
    hs(1) = atan2(obj.V(2,2),obj.V(1,2))*2/(2*pi);
    hs(2) = 1-exp(obj.Z(1));
    
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
    mnfb = mirroredNormalFisherBinghamDist('V',V,...
      'Z',Z,...
      'mu',[10,0,0]);
    figure(1), plot3(mnfb)
    mnfbCond = mnfb.conditional(1,10.5);
    figure(2), plot3(mnfbCond)
    mnfbAppx = mnfbCond.approximate;
    figure(3), plot3(mnfbAppx)
    keyboard
  end
end

end