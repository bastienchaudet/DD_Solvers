% Script to solve the Eta Minus Delta equation in two dimensions on the domain
% Omega = (xmin,xmax) x (ymin,ymax) with Jx x Jy grid points.


% --------------- %
% Initialization  %
% --------------- %

% Choice of the test case
Test = 0;
% Initialize all data
[f,eta,x,y,xmin,xmax,ymin,ymax,h,Jx,Jy,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex] = InitializeData2d_EtaMinusDelta(Test);


% -------------------------------------------------------------------------------------------------------------------------- %
% Solve the problem using different solvers                                                                                  %
% Direct:0                                                                                                                   %
% Two subdomains DD  --- Dirichlet-Neumann:1, Neumann-Neumann:2, (Robin-Robin:3)                                             %
% Four subdomains DD --- Dirichlet-Neumann(Seq):4, Dirichlet-Neumann(Par):5, Neumann-Neumann:6,                              %                           
%                        Dirichlet-Neumann(Seq,Mixed):7, Dirichlet-Neumann(Par,Mixed):8,
%                        Neumann-Neumann(Mixed):9                                                                            %                                    
% Multilevel DD      --- Dirichlet-Neumann(New):10, Dirichlet-Neumann(NewBis):11,                                            %
% -------------------------------------------------------------------------------------------------------------------------- %

Solvertype = 9;


% -------------- %
% Direct solver  %
% -------------- %
if (Solvertype == 0)
    u = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    uex = u;

    
% -------------------------- %
% Two subdomains DD solvers  %
% -------------------------- %

% Dirichlet-Neumann solver (sequential version)
% ---------------------------------------------
elseif (Solvertype == 1)
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end
    % Parameters for the DN method
    a = 0.5;                   % a in (0,1) : position of the vertical interface
    th = 0.5;                  % relaxation parameter
    Nit = 2;                   % number of iterations
    % Solve
    u = DNSeq2d_TwoSubdomains_EtaMinusDelta(f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,th,Nit,true,true);
   
% Neumann-Neumann solver
% -----------------------
elseif (Solvertype == 2)    
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end
    % Parameters for the NN method
    a = 0.5;                    % a in (0,1) : position of the vertical interface
    th = 0.25;                   % relaxation parameter
    Nit = 2;                    % number of iterations
    % Solve
    u = NN2d_TwoSubdomains_EtaMinusDelta(f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,th,Nit,false,true);    
   

% --------------------------- %
% Four subdomains DD solvers  %
% --------------------------- %

% Dirichlet-Neumann solver (sequential version)
% Choice of transmission conditions:
%   - Dirichlet for left bottom and top right 
%   - Neumann for right bottom and top left
% --------------------------------------------
elseif (Solvertype == 4)
    
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end    
    % Parameters for the DN method
    a = 0.5;                    % a in (0,1) : position of the vertical interface
    b = 0.5;                    % b in (0,1) : position of the horizontal interface
    th = 0.45;                  % relaxation parameter
    Nit = 2;                    % number of iterations
    % Solve
    u = DNSeq2d_FourSubdomains_EtaMinusDelta(f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,b,th,Nit,true,true);

% Dirichlet-Neumann solver (parallel version)
% Choice of transmission conditions:
%   - Dirichlet for left bottom and top right 
%   - Neumann for right bottom and top left
% --------------------------------------------
elseif (Solvertype == 5)
    
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end    
    % Parameters for the DN method
    a = 0.5;                    % a in (0,1) : position of the vertical interface
    b = 0.5;                    % b in (0,1) : position of the horizontal interface
    th = 0.5;                  % relaxation parameter
    Nit = 40;                    % number of iterations
    % Solve
    u = DNPar2d_FourSubdomains_EtaMinusDelta(f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,b,th,Nit,false,true);
    
% Neumann-Neumann solver
% -----------------------
elseif (Solvertype == 6)
    
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end
    % Parameters for the NN method
    a = 0.5;                    % a in (0,1) : position of the vertical interface
    b = 0.5;                    % b in (0,1) : position of the horizontal interface
    th = 0.25;                   % relaxation parameter
    Nit = 2;                    % number of iterations
    % Solve
    u = NN2d_FourSubdomains_EtaMinusDelta(f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,b,th,Nit,true,true);        
       
% Dirichlet-Neumann solver (sequential version)
% Choice of mixed transmission conditions for each subdomain
% ----------------------------------------------------------
elseif (Solvertype == 7)
    
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end
    % Parameters for the DN method
    a = 0.5;                    % a in (0,1) : position of the vertical interface
    b = 0.5;                    % b in (0,1) : position of the horizontal interface
    th = 0.5;                   % relaxation parameter
    Nit = 3;                    % number of iterations
    % Solve
    u = DNSeqMixed2d_FourSubdomains_EtaMinusDelta(f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,b,th,Nit,true,true); 
 
% Dirichlet-Neumann solver (parallel version)
% Choice of mixed transmission conditions for each subdomain
% ----------------------------------------------------------
elseif (Solvertype == 8)
    
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end
    % Parameters for the DN method
    a = 0.5;                    % a in (0,1) : position of the vertical interface
    b = 0.5;                    % b in (0,1) : position of the horizontal interface
    th = 0.5;                   % relaxation parameter
    Nit = 5;                    % number of iterations
    % Solve
    u = DNParMixed2d_FourSubdomains_EtaMinusDelta(f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,b,th,Nit,true,true);
    
% Neumann-Neumann solver
% Choice of mixed transmission conditions for each subdomain
% ----------------------------------------------------------
elseif (Solvertype == 9)
    
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end
    % Parameters for the NN method
    a = 0.5;                    % a in (0,1) : position of the vertical interface
    b = 0.5;                    % b in (0,1) : position of the horizontal interface
    th = 0.23;                   % relaxation parameter
    Nit = 8;                    % number of iterations
    % Solve
    u = NNMixed2d_FourSubdomains_EtaMinusDelta(f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,b,th,Nit,true,true);
    
    
% ---------------------- %
% Multilevel DD solvers  %
% ---------------------- %
 
% New Dirichlet-Neumann solver (sequential version)
% Based on a multilevel four-subdomains decomposition
% ---------------------------------------------------
elseif (Solvertype == 10)
    
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end
    % Parameters
    plot_evenodd = true;
    nb_level = 1;
    % Solve
    fprintf('------------------ Beginning of the DDM ------------------ \n');
    u = NewDN2d_EtaMinusDelta(nb_level,f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,plot_evenodd);
    fprintf('--------------------- End of the DDM --------------------- \n');
    
% NewBis Dirichlet-Neumann solver (sequential version)
% Based on a multilevel two-domain decomposition
% ----------------------------------------------------
elseif (Solvertype == 11)
    
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_EtaMinusDelta(f,eta,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end
    % Parameters
    nb_level = 2;
    % Solve
    [u,nb_solves] = NewDNBis2d_EtaMinusDelta(nb_level,f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,0);
    fprintf('Total number of solves on subdomains of smallest size : %d . \n', nb_solves);     
    
end


% ------------------------ %
% Analysis of the results  %
% ------------------------ %

% Display useful information
if (Solvertype == 0)
    fprintf('Direct solver used. \n');
elseif (Solvertype == 1)
    fprintf('Dirichlet-Neumann solver (two subdomains version) used, with %d iterations. \n', Nit);
elseif (Solvertype == 2)
    fprintf('Neumann-Neumann solver (two subdomains version) used, with %d iterations. \n', Nit);
elseif (Solvertype == 4)
    fprintf('Dirichlet-Neumann solver (four subdomains version) used, with %d iterations. \n', Nit);
elseif (Solvertype == 5)
    fprintf('Neumann-Neumann solver (four subdomains version) used, with %d iterations. \n', Nit);
end

% Error (if uex known or approximated)
if (true) %(uex ~= false)
    err = abs(u-uex)/norm(uex,'inf');
    fprintf('Infinity norm of the absolute error = %d \n', max(err,[],'all'));
end

% Plot the results
% The solution
figure('Name','Solution');
set(gcf, 'Color', 'w');
mesh(x,y,u);
xlabel('x');
ylabel('y');
zlabel('u');
% Source term
figure('Name','Source term');
set(gcf, 'Color', 'w');
mesh(x,y,f);
xlabel('x');
ylabel('y');
zlabel('f');
% The exact solution (if known) and the error
if (true) %(uex ~= false)
    % Exact solution
    figure('Name','Exact solution');
    set(gcf, 'Color', 'w');
    mesh(x,y,uex);
    xlabel('x');
    ylabel('y');
    zlabel('uex');
    % Exact error
    figure('Name','Absolute error');
    set(gcf, 'Color', 'w');
    mesh(x,y,err);
    xlabel('x');
    ylabel('y');
    zlabel('abs(u-uex)');
% % Else, an approximation of the error
% else
%     figure('Name','Absolute error');
%     mesh(x,y,abs(u-uex));
%     xlabel('x');
%     ylabel('y');
%     zlabel('abs(u-uex)');
end
