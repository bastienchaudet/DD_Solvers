function A = A2d_Helmholtz(k,xmin,xmax,ymin,ymax,Jx,Jy,BCtype,pl,pr,pb,pt)
% Computes A, a finite difference approximation of the two
% dimensional Helmholtz operator on the domain
% Omega=(xmin,xmax)x(ymin,ymax)using Jx x Jy points.
% BCtype : specifies the type of boundary conditions (0:Dirichlet or 1:Robin)
% Note that pl, pr, pb and pt are only required for Robin BC

    % Check that all arguments are valid
    if (BCtype ~= 0) && (BCtype ~= 1)
        error('Wrong argument for the boundary conditions');
    end
    
    % Check that the mesh sizes in x and y coincide
    if (abs((xmax-xmin)/(Jx-1)-(ymax-ymin)/(Jy-1)) > 1e-14)
        error('The mesh sizes in x and y do not coincide');
    end
    
    % Initialization
    h = (xmax-xmin)/(Jx-1);
    % k is replaced by k/sqrt(2) to ensure the validity of tensor formulas
    k2d = k/sqrt(2);
    
    % Dirichlet case
    if (BCtype == 0)
        Dxx = A1d_Helmholtz(k2d,xmin,xmax,Jx,0,0,0,0);
        Dyy = A1d_Helmholtz(k2d,ymin,ymax,Jy,0,0,0,0);
        A = kron(speye(size(Dxx)),Dyy)+kron(Dxx,speye(size(Dyy)));
            
    % Robin case
    else
        % Initialize standard FD matrices of size (Jx x Jy)
        Dxx = A1d_Helmholtz(k2d,xmin-h,xmax+h,Jx+2,0,0,0);
        Dyy = A1d_Helmholtz(k2d,ymin-h,ymax+h,Jy+2,0,0,0);
        % Modify those matrices at boundary indices
        Dxx(1,1) = 0.5*Dxx(1,1);
        Dxx(end,end) = 0.5*Dxx(end,end);
        Dyy(1,1) = 0.5*Dyy(1,1);
        Dyy(end,end) = 0.5*Dyy(end,end);
        
        % Vertical parts of the boundary
        Ix = speye(Jx);
        Ix(1,1) = 0.5;
        Ix(end,end) = 0.5;
        % Boundary condition on the left side
        Hxl = sparse(1,1,1/h,Jx,Jx);
        Pl = sparse(Jy,Jy);
        Pl = spdiags(pl,0,Pl);
        % Boundary condition on the right side
        Hxr = sparse(Jx,Jx,1/h,Jx,Jx);
        Pr = sparse(Jy,Jy);
        Pr = spdiags(pr,0,Pr);
       
        % Horizontal parts of the boundary
        Iy = speye(Jy);
        Iy(1,1) = 0.5;
        Iy(end,end) = 0.5;
        % Boundary condition on the bottom side
        Hxb = sparse(1,1,1/h,Jy,Jy);
        Pb = sparse(Jx,Jx);
        Pb = spdiags(pb,0,Pb);
        % Boundary condition on the top side
        Hxt = sparse(Jy,Jy,1/h,Jy,Jy);
        Pt = sparse(Jx,Jx);
        Pt = spdiags(pt,0,Pt);
        
        % Build A
        A = kron(Ix,Dyy)-kron(Hxl,Pl)-kron(Hxr,Pr)+kron(Dxx,Iy)-kron(Pb,Hxb)-kron(Pt,Hxt);
        
        % Correct corner indices
        % Bottom-left corner
        A(1,1) = 0.25*(k^2-4/h^2)-0.5*(pl(1)+pb(1))/h;
        % Bottom-right corner
        A(1+(Jx-1)*Jy,1+(Jx-1)*Jy) = 0.25*(k^2-4/h^2)-0.5*(pr(1)+pb(end))/h;
        % Top-left corner
        A(Jy,Jy) = 0.25*(k^2-4/h^2)-0.5*(pl(end)+pt(1))/h;
        % Top-right corner
        A(end,end) = 0.25*(k^2-4/h^2)-0.5*(pr(end)+pt(end))/h;

    end
    
end