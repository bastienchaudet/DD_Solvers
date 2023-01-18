function u = Solve1d_Helmholtz(f,k,a,b,gl,gr,FDscheme,BCtype,pl,pr)
% Solves the one dimensional equation (k^2+Delta)u=f on the domain
% Omega=(a,b)with Dirichlet boundary conditions u=gl at x=a and u=gr at x=b
% or Sommerfeld boundary conditions (dn+ik)u=0 at x=a and at x=b or Robin
% boundary conditions (dn+pl)u=gl at x=a and (dn+pr)u=gr at x=b using a
% finite difference approximation with length(f) grid points.
% FDscheme : specifies the finite difference scheme (0:standard or 1:new)
% BCtype : specifies the type of boundary conditions (0:Dirichlet or 1:Sommerfled or 2:Robin)
% Note that pg and pd are only required for Robin BC
    
    % Initialization
    J = length(f);
    A = A1d_Helmholtz(k,a,b,J,FDscheme,BCtype,pl,pr);
    h = (b-a)/(J-1);
        
    % Dirichlet case
    if (BCtype == 0)
        % Take the values of f only at interior points
        fint = f(2:end-1);
        % Add boundary conditions to the rhs
        fint(1) = fint(1)-gl/h^2;
        fint(end) = fint(end)-gr/h^2;
        % Solve
        u = A\fint;
        % Add values at the boundaries
        u = [gl; u; gr];
    
    % Sommerfeld case
    elseif (BCtype == 1)
        % Add boundary conditions to the rhs
        f(1) = 0.5*f(1);
        f(end) = 0.5*f(end);
        % Solve
        u = A\f; 
    
    % Robin case  
    else
        % Add boundary conditions to the rhs
        f(1) = 0.5*f(1)-gl/h;
        f(end) = 0.5*f(end)-gr/h;
        % Solve
        u = A\f;        
    end

end