classdef Variablemesh < handle
    %   Create a non-uniform mesh
    properties
        Dmax
        Length
        NX% normliazed x
        dx % the spacing between main point
        dx1 % the spacing between the middle point
        Lx% the position of the nodes
        zis % interface locations
        interfaces % p-i and i-n interface numbers
        Node % the layer node information
        N
        N1 % p-i interface
        N2 % i-n interface
    end
    properties
        Lx_half
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = Variablemesh(Dmax,Length,NX,zis,interfaces)
            self.Length = Length;
            self.Dmax = Dmax;
            self.zis= zis;
            self.NX = NX;
            self.interfaces = interfaces;
            [self.dx,self.dx1,self.Lx,self.N,self.Node]=self.createmesh(zis);
            self.dx=self.dx'/NX;
            self.dx1=self.dx1'/NX;
            self.Lx=self.Lx'/NX;
            self.Lx_half = (self.Lx(1:end-1)+self.Lx(2:end))/2;
             self.N1=self.Node(interfaces(1));
             self.N2=self.Node(interfaces(2));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [y,N] = create_nonuniform_mesh(self,X0,X1,h)
            n1 = round((X1-X0)/h);
            [p,t]=distmeshnd(2/n1,[-1;1],[]);
            p(1) = -1;
            p(end) = 1;
            y = X0+(X1-X0)/2*(p.'+1);
            N = length(y);
            function [p,t]=distmeshnd(h,box,fix,varargin)
                %DISTMESHND N-D Mesh Generator using Distance Functions.
                %   [P,T]=DISTMESHND(FDIST,FH,H,BOX,FIX,FDISTPARAMS)
                %
                %      P:           Node positions (NxNDIM)
                %      T:           Triangle indices (NTx(NDIM+1))
                %      FDIST:       Distance function
                %      FH:          Edge length function
                %      H:           Smallest edge length
                %      BOX:         Bounding box [xmin,xmax;ymin,ymax; ...] (NDIMx2)
                %      FIX:         Fixed node positions (NFIXxNDIM)
                %      FDISTPARAMS: Additional parameters passed to FDIST
                %
                %   Example: Unit ball
                %      dim=3;
                %      d=inline('sqrt(sum(p.^2,2))-1','p');
                %      [p,t]=distmeshnd(d,@huniform,0.2,[-ones(1,dim);ones(1,dim)],[]);
                %
                %   See also: DISTMESH2D, DELAUNAYN, TRIMESH, MESHDEMOND.
                
                %   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
                fdist = @(p) sqrt(sum(p.^2,2))-1;
                fh = @(p) sqrt(sum(p.^2,2))-2;
                dim=size(box,2);
                ptol=.001; ttol=.1; L0mult=1+.4/2^(dim-1); deltat=.1; geps=1e-1*h; deps=sqrt(eps)*h;
                
                % 1. Create initial distribution in bounding box
                p=(box(1):h:box(2))';
                
                % 2. Remove points outside the region, apply the rejection method
                p=p(feval(fdist,p,varargin{:})<geps,:);
                r0=feval(fh,p);
                p=[fix; p(rand(size(p,1),1)<min(r0)^dim./r0.^dim,:)];
                N=size(p,1);
                
                count=0;
                p0=inf;
                while 1
                    % 3. Retriangulation by Delaunay
                    if max(sqrt(sum((p-p0).^2,2)))>ttol*h
                        p0=p;
                        t=delaunayn(p);
                        pmid=zeros(size(t,1),dim);
                        for ii=1:dim+1
                            pmid=pmid+p(t(:,ii),:)/(dim+1);
                        end
                        t=t(feval(fdist,pmid,varargin{:})<-geps,:);
                        % 4. Describe each edge by a unique pair of nodes
                        pair=zeros(0,2);
                        localpairs=nchoosek(1:dim+1,2);
                        for ii=1:size(localpairs,1)
                            pair=[pair;t(:,localpairs(ii,:))];
                        end
                        pair=unique(sort(pair,2),'rows');
                        % 5. Graphical output of the current mesh
                        count=count+1;
                    end
                    
                    % 6. Move mesh points based on edge lengths L and forces F
                    bars=p(pair(:,1),:)-p(pair(:,2),:);
                    L=sqrt(sum(bars.^2,2));
                    L0=feval(fh,(p(pair(:,1),:)+p(pair(:,2),:))/2);
                    L0=L0*L0mult*(sum(L.^dim)/sum(L0.^dim))^(1/dim);
                    F=max(L0-L,0);
                    Fbar=[bars,-bars].*repmat(F./L,1,2*dim);
                    dp=full(sparse(pair(:,[ones(1,dim),2*ones(1,dim)]), ...
                        ones(size(pair,1),1)*[1:dim,1:dim], ...
                        Fbar,N,dim));
                    dp(1:size(fix,1),:)=0;
                    p=p+deltat*dp;
                    
                    % 7. Bring outside points back to the boundary
                    d=feval(fdist,p,varargin{:}); ix=d>0;
                    gradd=zeros(sum(ix),dim);
                    for ii=1:dim
                        a=zeros(1,dim);
                        a(ii)=deps;
                        d1x=feval(fdist,p(ix,:)+ones(sum(ix),1)*a,varargin{:});
                        gradd(:,ii)=(d1x-d(ix))/deps;
                    end
                    p(ix,:)=p(ix,:)-d(ix)*ones(1,dim).*gradd;
                    
                    % 8. Termination criterion
                    maxdp=max(deltat*sqrt(sum(dp(d<-geps,:).^2,2)));
                    if maxdp<ptol*h, break; end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dx,dx1,Lx,N,Node]=createmesh(self,zis)
            d2=self.Dmax;
            LX0 = [0 cumsum(zis)];
            ni = zeros(1, length(zis));
            % keyboard
            if LX0(2)-LX0(1)>60
                [y,ni(1)] = self.create_nonuniform_mesh(LX0(1),LX0(2),d2);
            else
                [y,ni(1)] = self.create_nonuniform_mesh(LX0(1),LX0(2),d2/2);
            end
            for ioi = 2:length(zis)
                if LX0(ioi+1)-LX0(ioi)>60
                   [yi,ni(ioi)] = self.create_nonuniform_mesh(LX0(ioi),LX0(ioi+1),d2);
                else
                   [yi,ni(ioi)] = self.create_nonuniform_mesh(LX0(ioi),LX0(ioi+1),d2/2); 
                end
                y = [y yi(2:end)];
            end
            y = y*1e-7;
            Lx=y;
            dx=diff(y);
            dx1=0.5*(dx(1:end-1)+dx(2:end));
            N=length(y);
            Node = cumsum(ni);
            Node(2:end) = Node(2:end)-(1:length(Node(2:end)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
end