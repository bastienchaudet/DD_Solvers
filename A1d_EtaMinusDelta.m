function A = A1d_EtaMinusDelta(eta,a,b,J,BCtype,pl,pr)
% Computes A, a finite difference approximation of the one
% dimensional operator eta-Delta on the domain Omega=(a,b) using J points.
% BCtype : specifies the type of boundary conditions (0:Dirichlet or 1:Robin)
% Note that pl and pr are only required for Robin BC

    % Check that all arguments are valid   
    if (BCtype ~= 0) && (BCtype ~= 1)
        error('Wrong argument for the boundary conditions');
    end
    
    % Initialization
    h = (b-a)/(J-1);
    
    % Dirichlet case 
    if (BCtype == 0)
        e = ones(J-2,1);
        A = spdiags([-e/h^2 (eta+2/h^2)*e -e/h^2],[-1 0 1],J-2,J-2); 
     
    % Robin case
    else
        e = ones(J,1);
        A = spdiags([-e/h^2 (eta+2/h^2)*e -e/h^2],[-1 0 1],J,J);
        A(1,1) = 0.5*A(1,1)+pl/h;
        A(end,end) = 0.5*A(end,end)+pr/h;
    end

end