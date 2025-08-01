classdef fd < handle
    %FD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh % this is another class of mesh
    end
    
    methods
        function [self] = fd(mesh)
            self.mesh = mesh;
        end
        function y=averagemp(self,x)
            % the function of solve the average of x
            
            y=(x(1:end-1)+x(2:end))/2;
        end
        function y = differ(self)
            y = 1./self.mesh.dx;
        end
        function y=diff1(self,x)
            y=x(1:end-1)-x(2:end);
        end
        function y=diffmp(self,x)
            % the funcation give the numerical difference
            y=(x(2:end)-x(1:end-1))./self.mesh.dx;
        end
        function y=diff3p(self,x)
           y=(x(3:end)-x(1:end-2))./(self.mesh.dx(1:end-1)+self.mesh.dx(2:end)); 
        end

        function y=diffmidp(self,x)
           % the funcation give the numerical difference for the mid point
           y=(x(2:end)-x(1:end-1))./self.mesh.dx1;
        end
        function y=diff2p(self,x)
            % second derivative
            temp=(x(2:end)-x(1:end-1))./self.mesh.dx;
            y=(temp(2:end)-temp(1:end-1))./self.mesh.dx1;
        end
        function M = CreateZOM(self)
            % create linear term jacob matrix
            dx = self.mesh.dx;
            N=length(dx)+1;
            temp = ones(1,N-2);
            M = spdiags(temp',0,N-2,N-2);
        end
        function M = CreateFMOM(self)
            % create first order difference(Mid) jacob matrix
            % this is used as when density do not have the value on the
            % mid point, so we use the average to represent the mid point.
            % so the first order difference is as dp(i)/dx=p(i+1)-p(i-1)/(2*dx)
            dx = self.mesh.dx;
            N=length(dx)+1;
            % Note we pad a at the start and b at the end by one
            % as the spdiags function ignores these values
            a = ones(N-2,1)./[(dx(1:N-3)+dx(2:N-2));1];% be careful
            b = ones(N-2,1)./[1;(dx(2:N-2)+dx(3:N-1))];% be careful
            %e = zeros(N,1);
            M = spdiags([a b],-1:2:1,N-2,N-2);
        end
        function M = CreateSOM_re(self,Eb)
            % create second order difference matrix
            % N is the total nodes, besides the boundary condition(two
            % nodes), there are N-2 nodes, so the jacob matrix is N-2*N-2
            dx = self.mesh.dx;
            dx1 = self.mesh.dx1;
            N=length(dx)+1;
            %Eb is relatively permittivity to InGaAs
            e = -(Eb(2:N-1)./dx(2:N-1)+Eb(1:N-2)./dx(1:N-2))./(dx1); 
            % Note we pad a at the start and b at the end by one
            % as the spdiags function ignores these values
            a = ones(N-2,1)./[(dx(2:N-2)./Eb(2:N-2).*dx1(2:N-2));1]; % be care
            b = ones(N-2,1)./[1;(dx(2:N-2)./Eb(2:N-2).*dx1(1:N-3))];
            M = spdiags([a e b],-1:1,N-2,N-2);
        end
        function M = CreateSOM(self)
            % create second order difference matrix
            % N is the total nodes, besides the boundary condition(two
            % nodes), there are N-2 nodes, so the jacob matrix is N-2*N-2
            dx = self.mesh.dx;
            dx1 = self.mesh.dx1;
            N=length(dx)+1;
            e = -(1./dx(2:N-1)+1./dx(1:N-2))./(dx1); 
            % Note we pad a at the start and b at the end by one
            % as the spdiags function ignores these values
            a = ones(N-2,1)./[(dx(2:N-2).*dx1(2:N-2));1]; % be care
            b = ones(N-2,1)./[1;(dx(2:N-2).*dx1(1:N-3))];
            M = spdiags([a e b],-1:1,N-2,N-2);
        end
        function M = CreateDiffAve(self)
            dx = self.mesh.dx;
            dx1 = self.mesh.dx1;
            N=length(dx)+1;
            % this is the function to create the Jacobian Matrix for use
            % the average to represent the mid point.ave(p)
            a = -ones(N-2,1)./[dx1(2:N-2);1]; % be care
            
            b = ones(N-2,1)./[1;dx1(1:N-3)];
            M = spdiags([a  b],-1:2:1,N-2,N-2);
        end
            
    end
    
end

