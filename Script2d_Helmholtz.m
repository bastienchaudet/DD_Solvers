% Script to solve the Helmholtz equation in two dimensions on the domain
% Omega = (xmin,xmax) x (ymin,ymax) with Jx x Jy grid points.

% --------------- %
% Initialization  %
% --------------- %

% Choice of the test case
Test = 5;
% Initialize all data
[f,k,x,y,xmin,xmax,ymin,ymax,h,Jx,Jy,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex] = InitializeData2d_Helmholtz(Test);

% ---------------------------------------------------------------- %
% Solve the problem using different solvers                        %
% Direct:0, Dirichlet-Neumann:1, Neumann-Neumann:2, Robin-Robin:3  %
% ---------------------------------------------------------------- %

Solvertype = 4;

% Direct solver
% -------------
if (Solvertype == 0)
    tic
    u = Solve2d_Helmholtz(f,k,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    toc
    
% Dirichlet-Neumann solver
% ------------------------
elseif (Solvertype == 1)
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_Helmholtz(f,k,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end
    tic
    % Choice of a in (0,1) that characterizes the position of the vertical interface
    a = 0.2;
    % Corresponding index and coordinate
    idxa = 1+round(a*(Jx-1));
    xa = xmin+h*(idxa-1);
    % Large value to emulate a Dirichlet BC
    lv = 1e12;
    % Subcomponents related to each subdomain
    x1 = x(:,1:idxa);
    x2 = x(:,idxa:end);
    y1 = y(:,1:idxa);
    y2 = y(:,idxa:end);
    f1 = f(:,1:idxa);
    f2 = f(:,idxa:end);
    u2 = zeros(size(f2));
    % Initial guess for u on the interface
    g = zeros(Jy,1);
    % In the case of Dirichlet BC, we can make a better guess
    if (BCtype == 0)
       g = gb(idxa)*ones(Jy,1)+(gt(idxa)-gb(idxa))*(0:1/(Jy-1):1)'; 
    end
    % Boundary conditions for the first (left) subdomain
    pb1 = pb(1:idxa);
    gb1 = gb(1:idxa);
    pt1 = pt(1:idxa);
    gt1 = gt(1:idxa);
    % Boundary conditions for the second (right) subdomain
    pl2 = zeros(Jy,1);
    if (BCtype == 0)
        pb2 = lv*ones(Jx-idxa+1,1);
        gb2 = lv*gb(idxa:end);
        pt2 = pb2;
        gt2 = lv*gt(idxa:end);
        pr = lv*ones(Jy,1);
        gr = lv*gr;
    else
        pb2 = pb(idxa:end);
        gb2 = gb(idxa:end);
        pt2 = pt(idxa:end);
        gt2 = gt(idxa:end);
    end
    % Relaxation parameter and number of iterations
    th = 0.5;
    Nit = 10;
    ErrL2 = zeros(Nit,1);
    ErrH1 = zeros(Nit,1);
    % Normal derivative operator with respect to n2 at the interface
    e = ones(Jy,1);
    NDer = [speye(Jy) -spdiags(0.5*[-e (4-k^2*h^2)*e -e],[-1 0 1],Jy,Jy)]/h;
    % Correct the values at the corners
    if (BCtype == 0)
        NDer(1,Jy+1) = -1/h;
        NDer(1,Jy+2) = 0;
        NDer(end,end) = -1/h;
        NDer(end,end-1) = 0;
    elseif (BCtype == 1)
        NDer(1,Jy+1) = NDer(1,Jy+1)-pb(idxa);
        NDer(1,Jy+2) = 2*NDer(1,Jy+2);
        NDer(end,end) = NDer(end,end)-pt(idxa);
        NDer(end,end-1) = 2*NDer(end,end-1);
    end
    % Fixed-point loop
    for i=1:Nit
        % Dirichlet step
        if (BCtype == 0)
            % Ensure that g is compatible with the other BC
            g(1) = gb1(end);
            g(end) = gt1(end);
            u1 = Solve2d_Helmholtz(f1,k,xmin,xa,ymin,ymax,gl,g,gb1,gt1,0,pl,pr,pb1,pt1);
        else
            u1 = Solve2d_Helmholtz(f1,k,xmin,xa,ymin,ymax,gl,lv*g,gb1,gt1,1,pl,lv*ones(Jy,1),pb1,pt1);
        end
        % Compute the normal derivative du1dn2
        gl2 = NDer*[u1(:,end-1);u1(:,end)]-0.5*h*f1(:,end);
        if (BCtype == 1)
            gl2(1) = gl2(1)+gb(idxa);
            gl2(end) = gl2(end)+gt(idxa);
        end
        % Neumann step
        u2 = Solve2d_Helmholtz(f2,k,xa,xmax,ymin,ymax,gl2,gr,gb2,gt2,1,pl2,pr,pb2,pt2);
        % Update g
        g = th*u2(:,1)+(1-th)*g;
        % Build the solution everywhere u
        u = [u1(:,1:end-1) u2];
        % Compute the L2 norm and the broken H1 norm
        ErrL2(i) = NormL2_FD_2d(h,u-uex)/NormL2_FD_2d(h,uex);
        ErrH1(i) = (NormH1_FD_2d(h,u1-uex(:,1:idxa))+NormH1_FD_2d(h,u2-uex(:,idxa:end)))/NormH1_FD_2d(h,uex);
%         % Plot the recombined solution
%         figure('Name','Real part of u after Neumann step');
%         mesh(x1,y1,real(u1)); hold on; mesh(x2,y2,real(u2)); hold off;
%         xlabel('x'); ylabel('y'); zlabel('Re(u2)');
%         pause
    end
    toc
    % Plot convergence history
    figure('Name','Convergence history');
    semilogy(1:1:Nit,ErrL2); hold on;
    semilogy(1:1:Nit,ErrH1);
    xlabel('iteration'); ylabel('err');
    legend('L^2-norm','Broken H^1-norm');
    
% Neumann-Neumann solver
% ----------------------
elseif (Solvertype == 2)
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_Helmholtz(f,k,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end
    tic
    % Choice of a in (0,1) that characterizes the position of the vertical interface
    a = 0.3;
    % Corresponding index and coordinate
    idxa = 1+round(a*(Jx-1));
    xa = xmin+h*(idxa-1);
    % Large value to emulate a Dirichlet BC
    lv = 1e12;
    % Subcomponents related to each subdomain
    x1 = x(:,1:idxa);
    x2 = x(:,idxa:end);
    y1 = y(:,1:idxa);
    y2 = y(:,idxa:end);
    f1 = f(:,1:idxa);
    f2 = f(:,idxa:end);
    u2 = zeros(size(f2));
    % Definition of useful vectors
    zerox1 = zeros(idxa,1);
    zerox2 = zeros(Jx-idxa+1,1);
    onex1 = ones(idxa,1);
    onex2 = ones(Jx-idxa+1,1);
    zeroy = zeros(Jy,1);
    oney = ones(Jy,1);    
    % Initial guess for u at the interface
    g = zeros(Jy,1);
    % In the case of Dirichlet BC, we can make a better guess
    if (BCtype == 0)
       g = gb(idxa)*oney+(gt(idxa)-gb(idxa))*(0:1/(Jy-1):1)'; 
    end
    % Boundary conditions for the first (left) subdomain
    pb1 = pb(1:idxa);
    gb1 = gb(1:idxa);
    pt1 = pt(1:idxa);
    gt1 = gt(1:idxa);
    % Boundary conditions for the second (right) subdomain
    pb2 = pb(idxa:end);
    gb2 = gb(idxa:end);
    pt2 = pt(idxa:end);
    gt2 = gt(idxa:end);
    % Relaxation parameter and number of iterations
    th = 0.25;
    Nit = 6;
    ErrL2 = zeros(Nit,1);
    ErrH1 = zeros(Nit,1);
    % Normal derivative operator at the interface
    e = ones(Jy,1);
    NDer = [-speye(Jy) spdiags(0.5*[-e (4-k^2*h^2)*e -e],[-1 0 1],Jy,Jy)]/h;
    % Correct the values at the corners
    if (BCtype == 0)
        NDer(1,Jy+1) = 1/h;
        NDer(1,Jy+2) = 0;
        NDer(end,end) = 1/h;
        NDer(end,end-1) = 0;
    elseif (BCtype == 1)
        NDer(1,Jy+1) = NDer(1,Jy+1)+pb(idxa);
        NDer(1,Jy+2) = 2*NDer(1,Jy+2);
        NDer(end,end) = NDer(end,end)+pt(idxa);
        NDer(end,end-1) = 2*NDer(end,end-1);
    end
    % Fixed-point loop
    for i=1:Nit
        % Dirichlet step
        if (BCtype == 0)
            % Ensure that g is compatible with the other BC
            g(1) = gb(idxa);
            g(end) = gt(idxa);
            u1 = Solve2d_Helmholtz(f1,k,xmin,xa,ymin,ymax,gl,g,gb1,gt1,0,pl,pr,pb1,pt1);
            u2 = Solve2d_Helmholtz(f2,k,xa,xmax,ymin,ymax,g,gr,gb2,gt2,0,pl,pr,pb2,pt2);
        else
            u1 = Solve2d_Helmholtz(f1,k,xmin,xa,ymin,ymax,gl,lv*g,gb1,gt1,1,pl,lv*oney,pb1,pt1);
            u2 = Solve2d_Helmholtz(f2,k,xa,xmax,ymin,ymax,lv*g,gr,gb2,gt2,1,lv*oney,pr,pb2,pt2);
        end
%         % Plot u on the domain
%         figure('Name','Real part of u after Dirichlet step');
%         mesh(x1,y1,real(u1)); hold on; mesh(x2,y2,real(u2)); hold off;
        % Compute the normal derivatives of u1 and u2
        du1dn1 = NDer*[u1(:,end-1);u1(:,end)]+0.5*h*f1(:,end);
        du2dn2 = NDer*[u2(:,2);u2(:,1)]+0.5*h*f2(:,1);
        if (BCtype == 1)
            du1dn1(1) = du1dn1(1)-gb(idxa);
            du1dn1(end) = du1dn1(end)-gt(idxa);
            du2dn2(1) = du2dn2(1)-gb(idxa);
            du2dn2(end) = du2dn2(end)-gt(idxa);
        end
        dpsidn = du1dn1+du2dn2;
        % Neumann step - compute the corrections psi1 and psi2
        if (BCtype == 0)
            psi1 = Solve2d_Helmholtz(zeros(size(f1)),k,xmin,xa,ymin,ymax,zeroy,dpsidn,zerox1,zerox1,1,lv*oney,zeroy,lv*onex1,lv*onex1);
            psi2 = Solve2d_Helmholtz(zeros(size(f2)),k,xa,xmax,ymin,ymax,dpsidn,zeroy,zerox2,zerox2,1,zeroy,lv*oney,lv*onex2,lv*onex2);
        else
            psi1 = Solve2d_Helmholtz(zeros(size(f1)),k,xmin,xa,ymin,ymax,zeroy,dpsidn,zerox1,zerox1,1,pl,zeroy,pb1,pt1);
            psi2 = Solve2d_Helmholtz(zeros(size(f2)),k,xa,xmax,ymin,ymax,dpsidn,zeroy,zerox2,zerox2,1,zeroy,pr,pb2,pt2);
        end
%         % Plot psi on the domain
%         figure('Name','Real part of psi after Neumann step');
%         mesh(x1,y1,real(psi1)); hold on; mesh(x2,y2,real(psi2)); hold off;
        % Update g
        g = g-th*(psi1(:,end)+psi2(:,1));
        % Build the solution everywhere u
        u = [u1(:,1:end-1) u2];
        % Compute the L2 norm and the broken H1 norm
        ErrL2(i) = NormL2_FD_2d(h,u-uex)/NormL2_FD_2d(h,uex);
        ErrH1(i) = (NormH1_FD_2d(h,u1-uex(:,1:idxa))+NormH1_FD_2d(h,u2-uex(:,idxa:end)))/NormH1_FD_2d(h,uex);
%        pause
    end
    toc    
    % Plot convergence history
    figure('Name','Convergence history');
    semilogy(1:1:Nit,ErrL2); hold on;
    semilogy(1:1:Nit,ErrH1);
    xlabel('iteration'); ylabel('err');
    legend('L^2-norm','Broken H^1-norm');

% Robin-Robin solver
% ------------------
elseif (Solvertype == 3)
    % Estimate the exact solution (if unknown) using the direct solver
    if (~uex)
        uex = Solve2d_Helmholtz(f,k,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end
    tic
    % Choice of a in (0,1) that characterizes the position of the vertical interface
    a = 0.5;
    % Corresponding index and coordinate
    idxa = 1+round(a*(Jx-1));
    xa = xmin+h*(idxa-1);
    % Large value to emulate a Dirichlet BC
    lv = 1e10;
    % Subcomponents related to each subdomain
    x1 = x(:,1:idxa);
    x2 = x(:,idxa:end);
    y1 = y(:,1:idxa);
    y2 = y(:,idxa:end);
    f1 = f(:,1:idxa);
    f2 = f(:,idxa:end);
    u2 = zeros(size(f2));
    % Definition of useful vectors
    onex1 = ones(idxa,1);
    onex2 = ones(Jx-idxa+1,1);
    oney = ones(Jy,1);  
    % Initial guess for u on the interface
    g1 = zeros(Jy,1);
    p = 1i*k*ones(Jy,1);
    % In the case of Dirichlet BC, we can make a better guess
    if (BCtype == 0)
       g1 = 1i*k*(gb(idxa)*ones(Jy,1)+(gt(idxa)-gb(idxa))*(0:1/(Jy-1):1)'); 
    end
    % Boundary conditions for the first (left) subdomain
    pb1 = pb(1:idxa);
    gb1 = gb(1:idxa);
    pt1 = pt(1:idxa);
    gt1 = gt(1:idxa);
    % Boundary conditions for the second (right) subdomain
    pb2 = pb(idxa:end);
    gb2 = gb(idxa:end);
    pt2 = pt(idxa:end);
    gt2 = gt(idxa:end);
    % Relaxation parameter and number of iterations
    th = 0.5;
    Nit = 30;
    ErrL2 = zeros(Nit,1);
    ErrH1 = zeros(Nit,1);
    % Fixed-point loop
    for i=1:Nit
        % First Robin step
        if (BCtype == 0)
            u1 = Solve2d_Helmholtz(f1,k,xmin,xa,ymin,ymax,lv*gl,g1,lv*gb1,lv*gt1,1,lv*oney,p,lv*onex1,lv*onex1);
        else
            u1 = Solve2d_Helmholtz(f1,k,xmin,xa,ymin,ymax,gl,g1,gb1,gt1,1,pl,p,pb1,pt1);
        end
        % Compute g2
        g2 = -g1+2*1i*k*u1(:,end);
        % Second Robin step
        if (BCtype == 0)
            u2 = Solve2d_Helmholtz(f2,k,xa,xmax,ymin,ymax,g2,lv*gr,lv*gb2,lv*gt2,1,p,lv*oney,lv*onex2,lv*onex2);
        else
            u2 = Solve2d_Helmholtz(f2,k,xa,xmax,ymin,ymax,g2,gr,gb2,gt2,1,p,pr,pb2,pt2);
        end
        % Update g1
        g1 = th*(-g2+2*1i*k*u2(:,1))+(1-th)*g1;
        % Build the solution everywhere u
        u = [u1(:,1:end-1) u2];
        % Compute the L2 norm and the broken H1 norm
        ErrL2(i) = NormL2_FD_2d(h,u-uex)/NormL2_FD_2d(h,uex);
        ErrH1(i) = (NormH1_FD_2d(h,u1-uex(:,1:idxa))+NormH1_FD_2d(h,u2-uex(:,idxa:end)))/NormH1_FD_2d(h,uex);
%         % Plot the recombined solution
%         figure('Name','Real part of u at current iteration');
%         mesh(x1,y1,real(u1)); hold on; mesh(x2,y2,real(u2)); hold off;
%         xlabel('x'); ylabel('y'); zlabel('Re(u2)');
%         pause
    end
    toc
    % Plot convergence history
    figure('Name','Convergence history');
    semilogy(1:1:Nit,ErrL2); hold on;
    semilogy(1:1:Nit,ErrH1);
    xlabel('iteration'); ylabel('err');
    legend('L^2-norm','Broken H^1-norm');
    
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
        uex = Solve2d_Helmholtz(f,k,xmin,xmax,ymin,ymax,gl,gr,gb,gt,BCtype,pl,pr,pb,pt);
    end    
    % Parameters for the DN method
    a = 0.5;                    % a in (0,1) : position of the vertical interface
    b = 0.5;                    % b in (0,1) : position of the horizontal interface
    th = 0.5;                  % relaxation parameter
    Nit = 2;                    % number of iterations
    % Solve
    u = DNSeq2d_FourSubdomains_Helmholtz(f,k,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,b,th,Nit,true,true);

end


% ------------------------ %
% Analysis of the results  %
% ------------------------ %

% Error (if uex known or approximated)
if (uex ~= false)
    err = abs(u-uex);
end

% Display useful information
fprintf('PPW = %d \n', 2*pi/(k*h));
if (Solvertype == 0)
    fprintf('Direct solver used. \n');
elseif (Solvertype == 1)
    fprintf('Dirichlet-Neumann solver used, with %d iterations. \n', Nit);
elseif (Solvertype == 2)
    fprintf('Neumann-Neumann solver used, with %d iterations. \n', Nit);
elseif (Solvertype == 3)
    fprintf('Robin-Robin solver used, with %d iterations. \n', Nit);
end
if (uex ~= false)
    fprintf('L2-norm of the error = %d \n', NormL2_FD_2d(h,u-uex));
end

% Plot the results
% The solution
figure('Name','Real part of the solution');
mesh(x,y,real(u));
xlabel('x');
ylabel('y');
zlabel('Re(u)');
% The exact solution (if known) and the error
if (uex ~= false)
    % Exact solution
    figure('Name','Real part of the exact solution');
    mesh(x,y,real(uex));
    xlabel('x');
    ylabel('y');
    zlabel('Re(uex)');
    % Exact error
    figure('Name','Module of the error');
    mesh(x,y,err);
    xlabel('x');
    ylabel('y');
    zlabel('abs(u-uex)');
end
