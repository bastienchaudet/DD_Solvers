function A = A1d_Helmholtz(k,a,b,J,FDscheme,BCtype,pl,pr)
% Computes A, a finite difference approximation of the one
% dimensional Helmholtz operator on the domain Omega=(a,b) using J points.
% FDscheme : specifies the finite difference scheme (0:standard or 1:new)
% BCtype : specifies the type of boundary conditions (0:Dirichlet or 1:Sommerfeld or 2:Robin)
% Note that pl and pr are only required for Robin BC
    
    % Check that all arguments are valid
    if (FDscheme ~= 0) && (FDscheme ~= 1)
        error('Wrong argument for the finite difference scheme');
    end    
    if (BCtype ~= 0) && (BCtype ~= 1) && (BCtype ~= 2)
        error('Wrong argument for the boundary conditions');
    end
    if (FDscheme == 1) && (BCtype == 2)
        error('The new FD scheme does not support Robin boundary conditions');
    end
        
    % Initialization
    h = (b-a)/(J-1);
    
    % Build A
    
    % Dirichlet case
    if (BCtype == 0)
        e = ones(J-2,1);
        if (FDscheme == 0)
            A = spdiags([e/h^2 (k^2-2/h^2)*e e/h^2],[-1 0 1],J-2,J-2); 
        else
            w = 2*cos(k*h)+(k*h)^2;
            A = spdiags([e/h^2 (k^2-w/h^2)*e e/h^2],[-1 0 1],J-2,J-2);
        end
        
    % Sommerfeld case
    elseif (BCtype == 1)
        e = ones(J,1);
        if (FDscheme == 0)
            A = spdiags([e/h^2 (k^2-2/h^2)*e e/h^2],[-1 0 1],J,J);
            A(1,1) = 0.5*A(1,1)-1i*k/h;
            A(end,end) = 0.5*A(end,end)-1i*k/h;
        else
            A = spdiags([e/h^2 (k^2-2/h^2)*e e/h^2],[-1 0 1],J,J);
            A(1,1) = 1;
            A(1,2) = -2+1i*sin(k*h)+cos(k*h);
            A(end,end-1) = -2+1i*sin(k*h)+cos(k*h);
            A(end,end) = 1;
        end
        
    % Robin case
    else
        e = ones(J,1);
        if (FDscheme == 0)
            A = spdiags([e/h^2 (k^2-2/h^2)*e e/h^2],[-1 0 1],J,J);
            A(1,1) = 0.5*A(1,1)-pl/h;
            A(end,end) = 0.5*A(end,end)-pr/h;
        end
    end    

end