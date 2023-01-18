function u = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt)
% Solves the bidimensional equation (eta-Delta)u=f on the domain
% Omega=(xmin,xmax)x(ymin,ymax)with Dirichlet boundary conditions u=gl at
% x=xmin, u=gr at x=xmax, u=gb at y=ymin, u=gt at y=ymax or Robin boundary 
% conditions (dn+pl)u=gl at x=xmin, (dn+pr)u=gr at x=xmax, (dn+pb)u=gb at 
% y=ymin and (dn+pt)u=gt at y=ymax using a finite difference approximation
% with size(f) grid points.
% BCtype : specifies the type of boundary conditions (0:Dirichlet or 1:Robin)
% Note that pl, pr, pb, pt are only required for Robin BC

    % Initialization
    Jx = size(f,2);
    Jy = size(f,1);
    h = (xmax-xmin)/(Jx-1);
    
    % Check the sizes of the arguments
    if (length(gl) ~= Jy) || (length(gr) ~= Jy)
        error('The dimensions of gl or gr are not compatible with the dimensions of f.');
    end
    if (length(gb) ~= Jx) || (length(gt) ~= Jx)
        error('The dimensions of gb or gt are not compatible with the dimensions of f.');
    end
    
    % Compatibility test between Dirichlet conditions
    if (BCtype == 0)
        if (abs(gb(1)-gl(1)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (bottom-left corner).');
        elseif (abs(gb(end)-gr(1)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (bottom-right corner).');
        elseif (abs(gt(1)-gl(end)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (top-left corner).');
        elseif (abs(gt(end)-gr(end)) > 1e-14)
            error('Compatibility problem in the Dirichlet data (top-right corner).');
        end
    end
    
    % Build A
    A = A2d_EtaMinusDelta(eta,xmin,xmax,ymin,ymax,Jx,Jy,BCtype,pl,pr,pb,pt);
        
    % Dirichlet case
    if (BCtype == 0)   
        % Add boundary conditions to the rhs
        f(2,2:end-1) = f(2,2:end-1)+gb(2:end-1).'/h^2;
        f(end-1,2:end-1) = f(end-1,2:end-1)+gt(2:end-1).'/h^2;
        f(2:end-1,2) = f(2:end-1,2)+gl(2:end-1)/h^2;
        f(2:end-1,end-1) = f(2:end-1,end-1)+gr(2:end-1)/h^2;
        % Solve
        fint = f(2:end-1,2:end-1);
        uint = A\fint(:);
        uint = reshape(uint,Jy-2,Jx-2);
        uintx = [gl(2:end-1) uint gr(2:end-1)];
        u = [gb.'; uintx; gt.'];
     
    % Robin case
    else
        % Add boundary conditions to the rhs
        % For the four corners
        f(1,1) = 0.25*f(1,1)+0.5*(gl(1)+gb(1))/h;
        f(1,end) = 0.25*f(1,end)+0.5*(gr(1)+gb(end))/h;
        f(end,1) = 0.25*f(end,1)+0.5*(gl(end)+gt(1))/h;
        f(end,end) = 0.25*f(end,end)+0.5*(gr(end)+gt(end))/h;
        % For the rest of the boundary
        f(2:end-1,1) = 0.5*f(2:end-1,1)+gl(2:end-1)/h;
        f(2:end-1,end) = 0.5*f(2:end-1,end)+gr(2:end-1)/h;
        f(1,2:end-1) = 0.5*f(1,2:end-1)+gb(2:end-1).'/h;
        f(end,2:end-1) = 0.5*f(end,2:end-1)+gt(2:end-1).'/h;
        % Solve
        u = A\f(:);
        u = reshape(u,Jy,Jx);
        
    end

end