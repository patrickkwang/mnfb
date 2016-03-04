classdef mirroredNormalFisherBinghamDist < normalFisherBinghamDist
% this generalizes mirroredFisherBinghamDist

methods
	function obj = mirroredNormalFisherBinghamDist(varargin)
    if nargout==0
      mirroredNormalFisherBinghamDist.test
      return
    end
    obj = obj@normalFisherBinghamDist(varargin{:});
  end
	
	function val = pdf(obj,points)
    val = exp(obj.logPdf(points));
  end
  
  function val = logPdf(obj,points)
    % points should be an nxd matrix
    
    for i = 2:length(obj.d)
      inds = sum(obj.d(1:i-1))+1:sum(obj.d(1:i));
      signs = sign(points(:,inds)*-obj.V(inds,sum(obj.d(1:i))));
      signs(signs==0) = 1;
      points(:,inds) = bsxfun(@times,normr(points(:,inds)),signs);
    end
    
    val = logPdf@normalFisherBinghamDist(obj,points);
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
    close all
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

%     Z = [-3, 0];
%     theta = -pi/3;
%     V = [cos(theta),-sin(theta);sin(theta),cos(theta)];
%     mnfb = mirroredNormalFisherBinghamDist('d',[0,2],...
%       'V',V,...
%       'Z',Z,...
%       'mu',[0,0]);
%     figure(1), plot(mnfb), zlim([0,0.5])
%     figure(2), plot1(mnfb,'k-')
%     figure(3), plot2(mnfb)
    
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
%     figure(1), plot(mnfb)%, axis equal
% %     samples = mnfb.sample(100);
% %     hold on, scatter3(samples(:,2),samples(:,3),samples(:,1)), hold off
%     mnfbCond = mnfb.conditional(1,9.3);
%     mnfbAppx = mnfbCond.approximate;
%     figure(2)
%     subplot(2,1,1), plot(mnfbCond), axis equal
%     subplot(2,1,2), plot(mnfbAppx), axis equal
%     mnfbMarg = mnfb.marginal([2,3]);
%     figure(3), plot(mnfbMarg)
		
    Z = [-2, -2, -20, 0];
		Ry = RotationMatrix('euler',[2.1847 0.9425 3.6825]); %rand(1,3)); %
		Ry = blkdiag(Ry.R,1);
		Rz = RotationMatrix.Rx(pi/3);
		Rz = blkdiag(1,Rz);
		V = Rz*Ry;
    mnfb = mirroredNormalFisherBinghamDist('V',V,...
      'Z',Z,...
      'mu',[10,10,0,0],...
			'd',[2,2]);
%     figure(1), plot(mnfb)
%     samples = mnfb.sample(100);
%     scale = max(max(samples(:,1))-min(samples(:,1)),max(samples(:,2))-min(samples(:,2)))/20;
%     samples(:,3:4) = [samples(:,3).^2-samples(:,4).^2,2*samples(:,3).*samples(:,4)];
%     hold on
%     arrow(samples(:,1),samples(:,2),samples(:,3)*scale,samples(:,4)*scale,scale/3,...
%       1*ones(100,3))
%     hold off
    mnfbCond = mnfb.conditional(1,10);
    figure(2), plot(mnfbCond), axis equal
    title('Conditioned on one Euclidean component')
    mnfbMarg = mnfb.marginal([3,4]);
    figure(3), plot(mnfbMarg)
    title('Marginalized over directional components')
  end
end

end