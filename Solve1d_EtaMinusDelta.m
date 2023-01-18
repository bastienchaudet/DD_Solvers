function u = Solve1d_EtaMinusDelta(f,eta,a,b,gl,gr,BCtype,pl,pr)
% Solves the one dimensional equation (eta-Delta)u=f on the domain
% Omega=(a,b)with Dirichlet boundary conditions u=gl at x=a and u=gr at x=b
% or Robin boundary conditions (dn+pl)u=gl at x=a and (dn+pr)u=gr at x=b
% using a finite difference approximation with length(f) grid points.
% BCtype : specifies the type of boundary conditions (0:Dirichlet or 1:Robin)
% Note that pl and pr are only required for Robin BC

    % Initialization
    J = length(f);
    A = A1d_EtaMinusDelta(eta,a,b,J,BCtype,pl,pr);
    h = (b-a)/(J-1);
    
    % Dirichlet case
    if (BCtype == 0)
        % Take the values of f only at interior points
        fint = f(2:end-1);
        % Add boundary conditions to the rhs
        fint(1) = fint(1)+gl/h^2;
        fint(end) = fint(end)+gr/h^2;
        % Solve
        u = A\fint;
        % Add values at the boundaries
        u = [gl; u; gr];
    
    % Robin case  
    else
        % Add boundary conditions to the rhs
        f(1) = 0.5*f(1)+gl/h;
        f(end) = 0.5*f(end)+gr/h;
        % Solve
        u = A\f;        
    end

end