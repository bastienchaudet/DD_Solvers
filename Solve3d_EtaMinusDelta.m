function u = Solve3d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,zmin,zmax,gl,gri,gb,gt,gre,gf,BCtype,pl,pri,pb,pt,pre,pf)
% Solves the tridimensional equation (eta-Delta)u=f on the domain
% Omega=(xmin,xmax)x(ymin,ymax)x(zmin,zmax) with Dirichlet boundary conditions u=gl at
% x=xmin, u=gri at x=xmax, u=gb at y=ymin, u=gt at y=ymax, u=gre at z=zmin and
% g=gf at z=zmax or Robin boundary conditions (dn+pl)u=gl at x=xmin, (dn+pri)u=gri
% at x=xmax, (dn+pb)u=gb at y=ymin, (dn+pt)u=gt at y=ymax, (dn+pre)u=gre at z=zmin
% and (dn+pf)u=gf at z=zmax using a finite difference approximation with
% size(f) grid points.
% BCtype : specifies the type of boundary conditions (0:Dirichlet or 1:Robin)
% Note that pl, pri, pb, pt, gre, gf are only required for Robin BC

    % Initialization
    Jx = size(f,2);
    Jy = size(f,1);
    Jz = size(f,3);
     
    h = (xmax-xmin)/(Jx-1);
    
    % Check the sizes of the arguments
    if (size(gl,1) ~= Jy) || (size(gl,2) ~= Jz) || (size(gri,1) ~= Jy) || (size(gri,2) ~= Jz)
        error('The dimensions of gl or gri are not compatible with the dimensions of f.');
    end
    if (size(gb,1) ~= Jx) || (size(gb,2) ~= Jz) || (size(gt,1) ~= Jx) || (size(gt,2) ~= Jz)
        error('The dimensions of gb or gt are not compatible with the dimensions of f.');
    end
    if (size(gre,1) ~= Jy) || (size(gre,2) ~= Jx) || (size(gf,1) ~= Jy) || (size(gf,2) ~= Jx)
        error('The dimensions of gre or gf are not compatible with the dimensions of f.');
    end
    
    % Compatibility test between Dirichlet conditions
    if (BCtype == 0)
        % Side edges
        if (norm(gb(1,:)-gl(1,:)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (bottom-left edge).');
        elseif (norm(gb(end,:)-gri(1,:)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (bottom-right edge).');
        elseif (norm(gt(1,:)-gl(end,:)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (top-left edge).');
        elseif (norm(gt(end,:)-gri(end,:)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (top-right edge).');
        % Front edges
        elseif (norm(gf(:,1)-gl(:,end)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (left-front edge).');
        elseif (norm(gf(:,end)-gri(:,end)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (right-front edge).');
        elseif (norm(gf(1,:)-gb(:,end)') > 1e-14)
            error('Compatibility problem in the Dirichlet data (bottom-front edge).');
        elseif (norm(gf(end,:)-gt(:,end)') > 1e-14)
            error('Compatibility problem in the Dirichlet data (top-front edge).');
        % Rear edges
        elseif (norm(gre(:,1)-gl(:,1)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (left-rear edge).');
        elseif (norm(gre(:,end)-gri(:,1)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (right-rear edge).');
        elseif (norm(gre(1,:)-gb(:,1)') > 1e-14)
            error('Compatibility problem in the Dirichlet data (bottom-rear edge).');
        elseif (norm(gre(end,:)-gt(:,1)') > 1e-14)
            error('Compatibility problem in the Dirichlet data (top-rear edge).'); 
        end
    end
    
    % Build A
    A = A3d_EtaMinusDelta(eta,xmin,xmax,ymin,ymax,zmin,zmax,Jx,Jy,Jz,BCtype,pl,pri,pb,pt,pre,pf);
        
    % Dirichlet case
    % --------------
    if (BCtype == 0)
        
        % Add boundary conditions to the rhs
        % Rear and front faces
        f(2:end-1,2:end-1,2) = f(2:end-1,2:end-1,2)+gre(2:end-1,2:end-1)/h^2;
        f(2:end-1,2:end-1,end-1) = f(2:end-1,2:end-1,end-1)+gf(2:end-1,2:end-1)/h^2;
        % Left and right faces
        fp = permute(f, [1 3 2]);
        fp(2:end-1,2:end-1,2) = fp(2:end-1,2:end-1,2)+gl(2:end-1,2:end-1)/h^2;
        fp(2:end-1,2:end-1,end-1) = fp(2:end-1,2:end-1,end-1)+gri(2:end-1,2:end-1)/h^2;
        f = ipermute(fp, [1 3 2]);
        % Bottom and top faces
        fp = permute(f, [2 3 1]);
        fp(2:end-1,2:end-1,2) = fp(2:end-1,2:end-1,2)+gb(2:end-1,2:end-1)/h^2;
        fp(2:end-1,2:end-1,end-1) = fp(2:end-1,2:end-1,end-1)+gt(2:end-1,2:end-1)/h^2;
        f = ipermute(fp, [2 3 1]);
        
        % Solve
        fint = f(2:end-1,2:end-1,2:end-1);
        uint = A\fint(:);
        uint = reshape(uint,Jy-2,Jx-2,Jz-2);
        
        % Add Dirichlet boundary conditions
        uintz = [gre(2:end-1,2:end-1) reshape(uint,Jy-2,(Jx-2)*(Jz-2)) gf(2:end-1,2:end-1)];
        uint = reshape(uintz,Jy-2,Jx-2,Jz);
        uintpx = permute(uint, [1 3 2]);
        uintx = [gl(2:end-1,:) reshape(uintpx,Jy-2,(Jx-2)*Jz) gri(2:end-1,:)];
        uint = ipermute(reshape(uintx,Jy-2,Jz,Jx),[1 3 2]);
        uintpy = permute(uint, [2 3 1]);
        uinty = [gb reshape(uintpy,Jx,(Jy-2)*Jz) gt];
        u = ipermute(reshape(uinty,Jx,Jz,Jy),[2 3 1]);
     
    % Robin case
    % ----------
    else
        
        % Add boundary conditions to the rhs
        % For the corners
        % Bottom-left-rear
        f(1,1,1) = 0.125*f(1,1,1)+0.25*(gl(1,1)+gb(1,1)+gre(1,1))/h;
        % Bottom-left-front
        f(1,1,Jz) = 0.125*f(1,1,Jz)+0.25*(gb(1,end)+gl(1,end)+gf(1,1))/h;
        % Bottom-right-rear
        f(1,Jx,1) = 0.125*f(1,Jx,1)+0.25*(gb(end,1)+gri(1,1)+gre(1,end))/h;
        % Bottom-right-front
        f(1,Jx,Jz) = 0.125*f(1,Jx,Jz)+0.25*(gb(end,end)+gri(1,end)+gf(1,end))/h;
        % Top-left-rear
        f(Jy,1,1) = 0.125*f(Jy,1,1)+0.25*(gt(1,1)+gl(end,1)+gre(end,1))/h;
        % Top-left-front
        f(Jy,1,Jz) = 0.125*f(Jy,1,Jz)+0.25*(gt(1,end)+gl(end,end)+gf(end,1))/h;
        % Top-right-rear
        f(Jy,Jx,1) = 0.125*f(Jy,Jx,1)+0.25*(gt(end,1)+gri(end,1)+gre(end,end))/h;
        % Top-right-front
        f(Jy,Jx,Jz) = 0.125*f(Jy,Jx,Jz)+0.25*(gt(end,end)+gri(end,end)+gf(end,end))/h;
        
        % For the edges
        % Along the z-axis
        fp = permute(f,[3 1 2]);
        fp(2:end-1,1,1) = 0.25*fp(2:end-1,1,1)+0.5*(gl(1,2:end-1)'+gb(1,2:end-1)')/h;
        fp(2:end-1,1,end) = 0.25*fp(2:end-1,1,end)+0.5*(gri(1,2:end-1)'+gb(end,2:end-1)')/h;
        fp(2:end-1,end,1) = 0.25*fp(2:end-1,end,1)+0.5*(gl(end,2:end-1)'+gt(1,2:end-1)')/h;
        fp(2:end-1,end,end) = 0.25*fp(2:end-1,end,end)+0.5*(gri(end,2:end-1)'+gt(end,2:end-1)')/h;
        f = ipermute(fp,[3 1 2]);
        % Along the x-axis
        fp = permute(f,[2 1 3]);
        fp(2:end-1,1,1) = 0.25*fp(2:end-1,1,1)+0.5*(gb(2:end-1,1)+gre(1,2:end-1)')/h;
        fp(2:end-1,end,1) = 0.25*fp(2:end-1,end,1)+0.5*(gt(2:end-1,1)+gre(end,2:end-1)')/h;
        fp(2:end-1,1,end) = 0.25*fp(2:end-1,1,end)+0.5*(gb(2:end-1,end)+gf(1,2:end-1)')/h;
        fp(2:end-1,end,end) = 0.25*fp(2:end-1,end,end)+0.5*(gt(2:end-1,end)+gf(end,2:end-1)')/h;
        f = ipermute(fp,[2 1 3]);
        % Along the y-axis
        f(2:end-1,1,1) = 0.25*f(2:end-1,1,1)+0.5*(gl(2:end-1,1)+gre(2:end-1,1))/h;
        f(2:end-1,end,1) = 0.25*f(2:end-1,end,1)+0.5*(gri(2:end-1,1)+gre(2:end-1,end))/h;
        f(2:end-1,1,end) = 0.25*f(2:end-1,1,end)+0.5*(gl(2:end-1,end)+gf(2:end-1,1))/h;
        f(2:end-1,end,end) = 0.25*f(2:end-1,end,end)+0.5*(gri(2:end-1,end)+gf(2:end-1,end))/h;
        
        % For the rest of the boundary (interior of the faces)
        % Rear and front faces
        f(2:end-1,2:end-1,1) = 0.5*f(2:end-1,2:end-1,1)+gre(2:end-1,2:end-1)/h;
        f(2:end-1,2:end-1,end) = 0.5*f(2:end-1,2:end-1,end)+gf(2:end-1,2:end-1)/h;
        % Left and right faces
        fp = permute(f, [1 3 2]);
        fp(2:end-1,2:end-1,1) = 0.5*fp(2:end-1,2:end-1,1)+gl(2:end-1,2:end-1)/h;
        fp(2:end-1,2:end-1,end) = 0.5*fp(2:end-1,2:end-1,end)+gri(2:end-1,2:end-1)/h;
        f = ipermute(fp, [1 3 2]);
        % Bottom and top faces
        fp = permute(f, [2 3 1]);
        fp(2:end-1,2:end-1,1) = 0.5*fp(2:end-1,2:end-1,1)+gb(2:end-1,2:end-1)/h;
        fp(2:end-1,2:end-1,end) = 0.5*fp(2:end-1,2:end-1,end)+gt(2:end-1,2:end-1)/h;
        f = ipermute(fp, [2 3 1]);

        % Solve
        u = A\f(:);
        u = reshape(u,Jy,Jx,Jz);
        
    end

end