function A = A3d_EtaMinusDelta(eta,xmin,xmax,ymin,ymax,zmin,zmax,Jx,Jy,Jz,BCtype,pl,pri,pb,pt,pre,pf)
% Computes A, a finite difference approximation of the three
% dimensional operator eta-Delta on the domain
% Omega=(xmin,xmax)x(ymin,ymax)x(zmin,zmax) using Jx x Jy x Jz points.

    % Check that all arguments are valid   
    if (BCtype ~= 0) && (BCtype ~= 1)
        error('Wrong argument for the boundary conditions');
    end
    
    % Check that the mesh sizes in x and y coincide
    if ( (abs((xmax-xmin)/(Jx-1)-(ymax-ymin)/(Jy-1)) > 1e-14) || (abs((xmax-xmin)/(Jx-1)-(zmax-zmin)/(Jz-1)) > 1e-14) )
        error('The mesh sizes in x, y and z do not coincide');
    end
    
    % Initialization
    h = (xmax-xmin)/(Jx-1);
    % eta is replaced by eta/3 to ensure the validity of tensor formulas
    eta3d = eta/3;
    
    % Dirichlet case
    % --------------
    if (BCtype == 0)
        
        Dxx = A1d_EtaMinusDelta(eta3d,xmin,xmax,Jx,0,0,0);
        Dyy = A1d_EtaMinusDelta(eta3d,ymin,ymax,Jy,0,0,0);
        Dzz = A1d_EtaMinusDelta(eta3d,zmin,zmax,Jz,0,0,0);
        A = kron(speye(size(Dzz)),kron(speye(size(Dxx)),Dyy)) + kron(speye(size(Dzz)),kron(Dxx,speye(size(Dyy)))) + kron(Dzz,kron(speye(size(Dxx)),speye(size(Dyy))));
        
    % Robin case
    % ----------
    else
        
        % Initialize standard FD matrices of size (Jx x Jy x Jz)
        Dxx = A1d_EtaMinusDelta(eta3d,xmin-h,xmax+h,Jx+2,0,0,0);
        Dyy = A1d_EtaMinusDelta(eta3d,ymin-h,ymax+h,Jy+2,0,0,0);
        Dzz = A1d_EtaMinusDelta(eta3d,zmin-h,zmax+h,Jz+2,0,0,0);
        
        % Modify those matrices at boundary indices
        Dxx(1,1) = 0.5*Dxx(1,1);
        Dxx(end,end) = 0.5*Dxx(end,end);
        Dyy(1,1) = 0.5*Dyy(1,1);
        Dyy(end,end) = 0.5*Dyy(end,end);
        Dzz(1,1) = 0.5*Dzz(1,1);
        Dzz(end,end) = 0.5*Dzz(end,end);  

        % Define useful quantities
        Ix = speye(Jx);
        Ix(1,1) = 0.5;
        Ix(end,end) = 0.5;
        Ox = sparse(Jx,Jx);
        Iy = speye(Jy);
        Iy(1,1) = 0.5;
        Iy(end,end) = 0.5;
        Oy = sparse(Jy,Jy);
        Iz = speye(Jz);
        Iz(1,1) = 0.5;
        Iz(end,end) = 0.5;
        Oz = sparse(Jz,Jz);

        % Left and right sides
        Hxl = sparse(1,1,1/h,Jx,Jx);
        Hxr = sparse(Jx,Jx,1/h,Jx,Jx);
        Axl = kron(Oz,kron(Ox,Oy));
        Axr = Axl;
        % Bottom and top sides
        Hyb = sparse(1,1,1/h,Jy,Jy);
        Hyt = sparse(Jy,Jy,1/h,Jy,Jy);
        Ayb = kron(Oz,kron(Ox,Oy));
        Ayt = Ayb;
        % Loop on the indices k to fill these 4 matrices
        for k = 1:Jz
            % Left
            Plk = spdiags(pl(:,k),0,Jy,Jy);
            Axl(1+(k-1)*Jx*Jy:k*Jx*Jy,1+(k-1)*Jx*Jy:k*Jx*Jy) = kron(Hxl,Plk);
            % Right
            Prk = spdiags(pri(:,k),0,Jy,Jy);
            Axr(1+(k-1)*Jx*Jy:k*Jx*Jy,1+(k-1)*Jx*Jy:k*Jx*Jy) = kron(Hxr,Prk);
            % Bottom
            Pbk = spdiags(pb(:,k),0,Jx,Jx);
            Ayb(1+(k-1)*Jx*Jy:k*Jx*Jy,1+(k-1)*Jx*Jy:k*Jx*Jy) = kron(Pbk,Hyb);
            % Top
            Ptk = spdiags(pt(:,k),0,Jx,Jx);
            Ayt(1+(k-1)*Jx*Jy:k*Jx*Jy,1+(k-1)*Jx*Jy:k*Jx*Jy) = kron(Ptk,Hyt);
        end

        % Rear and front sides
        Hzr = sparse(1,1,1/h,Jz,Jz);
        Hzf = sparse(Jz,Jz,1/h,Jz,Jz);
        pre_t = pre';
        Pre = spdiags(pre_t(:),0,Jx*Jy,Jx*Jy);
        pf_t = pf';
        Pf = spdiags(pf_t(:),0,Jx*Jy,Jx*Jy);
        Azr = kron(Hzr,Pre);
        Azf = kron(Hzf,Pf);

        % Build A
        A =  kron(Iz,kron(Ix,Dyy)) + kron(Iz,kron(Dxx,Iy)) + kron(Dzz,kron(Ix,Iy)) + Axl + Axr + Ayb + Ayt + Azr + Azf; 

        % Correct edge indices
        % Edges along the z-axis
        for k = 2:Jz-1
            % Bottom-left
            A(1+(k-1)*Jx*Jy,1+(k-1)*Jx*Jy) = 0.25*(eta+6/h^2) + 0.5*(pl(1,k)+pb(1,k))/h;
            % Bottom-right
            A(1+(Jx-1)*Jy+(k-1)*Jx*Jy,1+(Jx-1)*Jy+(k-1)*Jx*Jy) = 0.25*(eta+6/h^2) + 0.5*(pri(1,k)+pb(end,k))/h;
            % Top-left
            A(Jy+(k-1)*Jx*Jy,Jy+(k-1)*Jx*Jy) = 0.25*(eta+6/h^2) + 0.5*(pl(end,k)+pt(1,k))/h;
            % Top-right
            A(k*Jx*Jy,k*Jx*Jy) = 0.25*(eta+6/h^2) + 0.5*(pri(end,k)+pt(end,k))/h;
        end
        % Edges along the x-axis
        for j = 2:Jx-1
            % Bottom-rear
            A(1+(j-1)*Jy,1+(j-1)*Jy) = 0.25*(eta+6/h^2) + 0.5*(pb(j,1)+pre(1,j))/h;
            % Bottom-front
            A(1+(Jz-1)*Jx*Jy+(j-1)*Jy,1+(Jz-1)*Jx*Jy+(j-1)*Jy) = 0.25*(eta+6/h^2) + 0.5*(pb(j,end)+pf(1,j))/h;
            % Top-rear
            A(j*Jy,j*Jy) = 0.25*(eta+6/h^2) + 0.5*(pt(j,1)+pre(end,j))/h;
            % Top-front
            A((Jz-1)*Jx*Jy+j*Jy,(Jz-1)*Jx*Jy+j*Jy) = 0.25*(eta+6/h^2) + 0.5*(pt(j,end)+pf(end,j))/h;
        end
        % Edges along the y-axis
        for i = 2:Jy-1
            % Left-rear
            A(i,i) = 0.25*(eta+6/h^2) + 0.5*(pl(i,1)+pre(i,1))/h; 
            % Left-front
            A((Jz-1)*Jx*Jy+i,(Jz-1)*Jx*Jy+i) = 0.25*(eta+6/h^2) + 0.5*(pl(i,end)+pf(i,1))/h;
            % Right-rear
            A((Jx-1)*Jy+i,(Jx-1)*Jy+i) = 0.25*(eta+6/h^2) + 0.5*(pri(i,1)+pre(i,end))/h;
            % Right-front
            A((Jz-1)*Jx*Jy+(Jx-1)*Jy+i,(Jz-1)*Jx*Jy+(Jx-1)*Jy+i) = 0.25*(eta+6/h^2) + 0.5*(pri(i,end)+pf(i,end))/h;
        end

        % Correct corner indices
        % Bottom-left-rear
        A(1,1) = 0.125*(eta+6/h^2) + 0.25*(pb(1,1)+pl(1,1)+pre(1,1))/h;
        % Bottom-left-front
        A(1+(Jz-1)*Jx*Jy,1+(Jz-1)*Jx*Jy) = 0.125*(eta+6/h^2) + 0.25*(pb(1,end)+pl(1,end)+pf(1,1))/h;
        % Bottom-right-rear
        A(1+(Jx-1)*Jy,1+(Jx-1)*Jy) = 0.125*(eta+6/h^2) + 0.25*(pb(end,1)+pri(1,1)+pre(1,end))/h;
        % Bottom-right-front
        A(1+(Jz-1)*Jx*Jy+(Jx-1)*Jy,1+(Jz-1)*Jx*Jy+(Jx-1)*Jy) = 0.125*(eta+6/h^2) + 0.25*(pb(end,end)+pri(1,end)+pf(1,end))/h;
        % Top-left-rear
        A(Jy,Jy) = 0.125*(eta+6/h^2) + 0.25*(pt(1,1)+pl(end,1)+pre(end,1))/h;
        % Top-left-front
        A((Jz-1)*Jx*Jy+Jy,(Jz-1)*Jx*Jy+Jy) = 0.125*(eta+6/h^2) + 0.25*(pt(1,end)+pl(end,end)+pf(end,1))/h;
        % Top-right-rear
        A(Jx*Jy,Jx*Jy) = 0.125*(eta+6/h^2) + 0.25*(pt(end,1)+pri(end,1)+pre(end,end))/h;
        % Top-right-front
        A(end,end) = 0.125*(eta+6/h^2) + 0.25*(pt(end,end)+pri(end,end)+pf(end,end))/h;

    end


end