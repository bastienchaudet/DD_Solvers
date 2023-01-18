% Script to solve the Eta Minus Delta equation in three dimensions on the domain
% Omega = (xmin,xmax) x (ymin,ymax) with Jx x Jy x Jz grid points.


% --------------- %
% Initialization  %
% --------------- %

% Choice of the test case
Test = 4;
% Initialize all data
[f,eta,x,y,z,xmin,xmax,ymin,ymax,zmin,zmax,h,Jx,Jy,Jz,gl,gri,gb,gt,gre,gf,BCtype,pl,pri,pb,pt,pre,pf,uex] = InitializeData3d_EtaMinusDelta(Test);

% -------------------------------------------------------------------------------------------------------------------------- %
% Solve the problem using different solvers                                                                                  %
% Direct:0                                                                                                                   %
% Two subdomains DD  --- Dirichlet-Neumann:1, Neumann-Neumann:2, (Robin-Robin:3)                                             %
% Four subdomains DD --- Dirichlet-Neumann(Seq):4, Dirichlet-Neumann(Par):5, Neumann-Neumann:6,                              %                           
%                        Dirichlet-Neumann(Seq,Mixed):7, Dirichlet-Neumann(Par,Mixed):8,
%                        Neumann-Neumann(Mixed):9                                                                            %                                    
% Multilevel DD      --- Dirichlet-Neumann(New):10, Dirichlet-Neumann(NewBis):11,                                            %
% -------------------------------------------------------------------------------------------------------------------------- %

Solvertype = 5;


% -------------- %
% Direct solver  %
% -------------- %
if (Solvertype == 0)
    u = Solve3d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,zmin,zmax,gl,gri,gb,gt,gre,gf,BCtype,pl,pri,pb,pt,pre,pf);
    uex = u;
    

% -------------------------- %
% Two subdomains DD solvers  %
% -------------------------- %

% Dirichlet-Neumann solver (sequential version)
% ---------------------------------------------
elseif (Solvertype == 1)
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve3d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,zmin,zmax,gl,gri,gb,gt,gre,gf,BCtype,pl,pri,pb,pt,pre,pf);
    end
    % Parameters for the DN method
    a = 0.5;                   % a in (0,1) : position of the vertical interface
    th = 0.5;                 % relaxation parameter
    Nit = 2;                   % number of iterations
    % Solve
    u = DNSeq3d_TwoSubdomains_EtaMinusDelta(f,eta,x,y,z,gl,gri,gb,gt,gre,gf,BCtype,pl,pri,pb,pt,pre,pf,uex,a,th,Nit,true,true);
 
    
% --------------------------- %
% Four subdomains DD solvers  %
% --------------------------- %

% Dirichlet-Neumann solver for even symmetric problems
% Choice of transmission conditions:
%   - Dirichlet for left bottom and top right 
%   - Neumann for right bottom and top left
% --------------------------------------------
elseif (Solvertype == 4)
    
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve3d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,zmin,zmax,gl,gri,gb,gt,gre,gf,BCtype,pl,pri,pb,pt,pre,pf);
    end    
    % Solve
    u = DNEven3d_FourSubdomains_EtaMinusDelta(f,eta,x,y,z,gl,gb,gre,gf,BCtype,pl,pb,pre,pf);
    
% Dirichlet-Neumann solver for odd symmetric problems
% Choice of mixed transmission conditions:
% --------------------------------------------
elseif (Solvertype == 5)
    
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve3d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,zmin,zmax,gl,gri,gb,gt,gre,gf,BCtype,pl,pri,pb,pt,pre,pf);
    end    
    % Solve
    u = DNOdd3d_FourSubdomains_EtaMinusDelta(f,eta,x,y,z,gl,gb,gre,gf,BCtype,pl,pb,pre,pf);    
    
end
 
% ------------------------ %
% Analysis of the results  %
% ------------------------ %

% Display useful information
if (Solvertype == 0)
    fprintf('Direct solver used. \n');
end

% Error (if uex known or approximated)
if (true)%(uex ~= false)
    err = (u-uex)/norm(uex(:),'inf');
    fprintf('Infinity norm of the relative error = %d \n', norm(err(:),Inf));
    %fprintf('L2 norm of the relative error = %d \n', norm(err(:))/length(err(:)));
end

% Plot the results
% Define the planes to display the solution
xs0 = xmin+(xmax-xmin)/2; 
xs1 = xmin+(xmax-xmin)/3;
xs2 = xmin+2*(xmax-xmin)/3;
ys0 = ymin+(ymax-ymin)/2;
ys1 = ymin+(ymax-ymin)/3;
ys2 = ymin+2*(ymax-ymin)/3;
zs0 = zmin+(zmax-zmin)/2;
% Define the slices
% xslice = [xs1,xs2];   
% yslice = [ys1,ys2];
% zslice = zs0;
xslice = [xmin,xmax];
yslice = [ymin,ymax];
zslice = [zmin,zmax];
% The solution
figure('Name','Solution');
set(gcf, 'Color', 'w');
slice(x,y,z,u,xslice,yslice,zslice);
% The exact solution (if known) and the error
if (true)%(uex ~= false)
    % Exact solution
    figure('Name','Exact solution');
    slice(x,y,z,uex,xslice,yslice,zslice);
    % Exact error
    figure('Name','Relative error');
    set(gcf, 'Color', 'w');
    slice(x,y,z,abs(err),xslice,yslice,zslice);
end
