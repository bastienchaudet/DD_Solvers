function u = NN2d_TwoSubdomains_EtaMinusDelta(f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,th,Nit,plot_iterates,plot_cvg_curve)

    % Geometric data
    xmin = x(1,1);
    xmax = x(1,end);
    ymin = y(1,1);
    ymax = y(end,1);
    [Jy,Jx] = size(x);
    h = x(1,2)-x(1,1);
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
    ErrL2 = zeros(Nit,1);
    ErrH1 = zeros(Nit,1);
    % Normal derivative operator at the interface
    e = ones(Jy,1);
    NDer = [speye(Jy) -spdiags(0.5*[-e (4+eta*h^2)*e -e],[-1 0 1],Jy,Jy)]/h;
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
        % --------------
        if (BCtype == 0)
            % Ensure that g is compatible with the other BC
            g(1) = gb(idxa);
            g(end) = gt(idxa);
            u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,ymax,gl,g,gb1,gt1,0,pl,pr,pb1,pt1);
            u2 = Solve2d_EtaMinusDelta(f2,eta,xa,xmax,ymin,ymax,g,gr,gb2,gt2,0,pl,pr,pb2,pt2);
        else
            u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,ymax,gl,lv*g,gb1,gt1,1,pl,lv*oney,pb1,pt1);
            u2 = Solve2d_EtaMinusDelta(f2,eta,xa,xmax,ymin,ymax,lv*g,gr,gb2,gt2,1,lv*oney,pr,pb2,pt2);
        end
        
        % Plot u on the domain
        if (plot_iterates)
            figure('Name','u after Dirichlet step');
            mesh(x1,y1,u1); hold on; mesh(x2,y2,u2); hold off;
        end
        
        % Compute the normal derivatives of u1 and u2
        du1dn1 = -NDer*[u1(:,end-1);u1(:,end)]-0.5*h*[0; f1(2:end-1,end); 0];
        du2dn2 = -NDer*[u2(:,2);u2(:,1)]-0.5*h*[0; f2(2:end-1,1); 0];
        if (BCtype == 1)
            du1dn1(1) = du1dn1(1)-gb(idxa)-0.5*h*f1(1,end);
            du1dn1(end) = du1dn1(end)-gt(idxa)-0.5*h*f1(end,end);
            du2dn2(1) = du2dn2(1)-gb(idxa)-0.5*h*f2(1,1);
            du2dn2(end) = du2dn2(end)-gt(idxa)-0.5*h*f2(end,1);
        end
        dpsidn = du1dn1+du2dn2;
        
        % Neumann step - compute the corrections psi1 and psi2
        % ----------------------------------------------------
        if (BCtype == 0)
            psi1 = Solve2d_EtaMinusDelta(zeros(size(f1)),eta,xmin,xa,ymin,ymax,zeroy,dpsidn,zerox1,zerox1,1,lv*oney,zeroy,lv*onex1,lv*onex1);
            psi2 = Solve2d_EtaMinusDelta(zeros(size(f2)),eta,xa,xmax,ymin,ymax,dpsidn,zeroy,zerox2,zerox2,1,zeroy,lv*oney,lv*onex2,lv*onex2);
        else
            psi1 = Solve2d_EtaMinusDelta(zeros(size(f1)),eta,xmin,xa,ymin,ymax,zeroy,dpsidn,zerox1,zerox1,1,pl,zeroy,pb1,pt1);
            psi2 = Solve2d_EtaMinusDelta(zeros(size(f2)),eta,xa,xmax,ymin,ymax,dpsidn,zeroy,zerox2,zerox2,1,zeroy,pr,pb2,pt2);
        end
        
        % Plot psi on the domain
        if (plot_iterates)
            figure('Name','Psi after Neumann step');
            mesh(x1,y1,psi1); hold on; mesh(x2,y2,psi2); hold off;
            pause
        end
        
        % Update g
        % --------
        g = g-th*(psi1(:,end)+psi2(:,1));
        
        % Build the solution everywhere u
        u = [u1(:,1:end-1) u2];
        
        % Compute the L2 norm and the broken H1 norm
        ErrL2(i) = NormL2_FD_2d(h,u-uex)/NormL2_FD_2d(h,uex);
        ErrH1(i) = (NormH1_FD_2d(h,-u1+uex(:,1:idxa))+NormH1_FD_2d(h,-u2+uex(:,idxa:end)))/NormH1_FD_2d(h,uex);

    end
    
    % Plot convergence history
    if (plot_cvg_curve)
        figure('Name','Convergence history');
        semilogy(1:1:Nit,ErrL2); hold on;
        semilogy(1:1:Nit,ErrH1);
        xlabel('iteration'); ylabel('err');
        legend('L^2-norm','Broken H^1-norm');
    end
    
end
