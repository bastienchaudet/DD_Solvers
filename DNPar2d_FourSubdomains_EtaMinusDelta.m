function u = DNPar2d_FourSubdomains_EtaMinusDelta(f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,b,th,Nit,plot_iterates,plot_cvg_curve)

    % Geometric data
    xmin = x(1,1);
    xmax = x(1,end);
    ymin = y(1,1);
    ymax = y(end,1);
    [Jy,Jx] = size(x);
    h = x(1,2)-x(1,1);
    % Corresponding indices and coordinates
    idxa = 1+round(a*(Jx-1));
    idxb = 1+round(b*(Jy-1));
    xa = xmin+h*(idxa-1);
    yb = ymin+h*(idxb-1);
    % Large value to emulate a Dirichlet BC
    lv = 1e10;
    % Subcomponents related to each subdomain (numbering starting from left bottom then anticlockwise)
    x1 = x(1:idxb,1:idxa);
    x2 = x(1:idxb,idxa:end);
    x3 = x(idxb:end,idxa:end);
    x4 = x(idxb:end,1:idxa);
    y1 = y(1:idxb,1:idxa);
    y2 = y(1:idxb,idxa:end);
    y3 = y(idxb:end,idxa:end);
    y4 = y(idxb:end,1:idxa);
    f1 = f(1:idxb,1:idxa);
    f2 = f(1:idxb,idxa:end);
    f3 = f(idxb:end,idxa:end);
    f4 = f(idxb:end,1:idxa);
    % Initial guess for u on the interface
    g12 = zeros(idxb,1);
    g23 = zeros(Jx-idxa+1,1);
    g34 = zeros(Jy-idxb+1,1);
    g41 = zeros(idxa,1);
    % In the case of Dirichlet BC, we can make a better guess
    if (BCtype == 0)
        gcp_vert = (1-b)*gb(idxa)+b*gt(idxa);
        gcp_hor = (1-a)*gl(idxb)+a*gr(idxb);
        gcp = 0.5*(gcp_vert+gcp_hor);
        g12 = gb(idxa)*ones(idxb,1)+(gcp-gb(idxa))*(0:1/(idxb-1):1)'; 
        g23 = gcp*ones(Jx-idxa+1,1)+(gr(idxb)-gcp)*(0:1/(Jx-idxa):1)'; 
        g34 = gcp*ones(Jy-idxb+1,1)+(gt(idxa)-gcp)*(0:1/(Jy-idxb):1)'; 
        g41 = gl(idxb)*ones(idxa,1)+(gcp-gl(idxb))*(0:1/(idxa-1):1)'; 
    end    
    % Boundary conditions for the first (left bottom) subdomain
    pb1 = pb(1:idxa);
    gb1 = gb(1:idxa);
    gl1 = gl(1:idxb);
    pl1 = pl(1:idxb);
    pt1 = zeros(idxa,1);
    pr1 = zeros(idxb,1);
    % Boundary conditions for the second (right bottom) subdomain
    gl2 = zeros(idxb,1);
    gt2 = zeros(Jx-idxa+1,1);
    pl2 = zeros(idxb,1);
    pt2 = zeros(Jx-idxa+1,1);
    if (BCtype == 0)
        pb2 = lv*ones(Jx-idxa+1,1);
        gb2 = lv*gb(idxa:end);
        pr2 = lv*ones(idxb,1);
        gr2 = lv*gr(1:idxb);
    else
        pb2 = pb(idxa:end);
        gb2 = gb(idxa:end);
        pr2 = pr(1:idxb);
        gr2 = gr(1:idxb);
    end
    % Boundary conditions for the third (top right) subdomain
    pt3 = pt(idxa:end);
    gt3 = gt(idxa:end);
    gr3 = gr(idxb:end);
    pr3 = pr(idxb:end);
    pb3 = zeros(Jx-idxa+1,1);
    pl3 = zeros(Jy-idxb+1,1);
    % Boundary conditions for the fourth (top left) subdomain
    gb4 = zeros(idxa,1);
    gr4 = zeros(Jy-idxb+1,1);
    pb4 = zeros(idxa,1);
    pr4 = zeros(Jy-idxb+1,1);
    if (BCtype == 0)
        pt4 = lv*ones(idxa,1);
        gt4 = lv*gt(1:idxa);
        pl4 = lv*ones(Jy-idxb+1,1);
        gl4 = lv*gl(idxb:end);
    else
        pt4 = pt(1:idxa);
        gt4 = gt(1:idxa);
        pl4 = pl(idxb:end);
        gl4 = gl(idxb:end);
    end
    % Initialization for the norms of the error
    ErrL2 = zeros(Nit,1);
    ErrH1 = zeros(Nit,1);
    E22_0 = zeros(Nit,1);
    Slope = zeros(Nit,1);
    % Normal derivative operator at the interface Gamma12
    e12 = ones(idxb,1);
    NDer12 = [speye(idxb) -spdiags(0.5*[-e12 (4+eta*h^2)*e12 -e12],[-1 0 1],idxb,idxb)]/h;
    % Correct the values at the corners
    NDer12(end,end-1) = 2*NDer12(end,end-1);
    if (BCtype == 0)
        NDer12(1,idxb+1) = -1/h;
        NDer12(1,idxb+2) = 0;
    elseif (BCtype == 1)
        NDer12(1,idxb+1) = NDer12(1,idxb+1)-pb(idxa);
        NDer12(1,idxb+2) = 2*NDer12(1,idxb+2);
    end
    % Normal derivative operator at the interface Gamma23
    e23 = ones(Jx-idxa+1,1);
    NDer23 = [speye(Jx-idxa+1) -spdiags(0.5*[-e23 (4+eta*h^2)*e23 -e23],[-1 0 1],Jx-idxa+1,Jx-idxa+1)]/h;
    % Correct the values at the corners
    NDer23(1,Jx-idxa+3) = 2*NDer23(1,Jx-idxa+3);
    if (BCtype == 0)
        NDer23(end,end) = -1/h;
        NDer23(end,end-1) = 0;
    elseif (BCtype == 1)
        NDer23(end,end) = NDer23(end,end)-pr(idxb);
        NDer23(end,end-1) = 2*NDer23(end,end-1);
    end 
    % Normal derivative operator at the interface Gamma34
    e34 = ones(Jy-idxb+1,1);
    NDer34 = [speye(Jy-idxb+1) -spdiags(0.5*[-e34 (4+eta*h^2)*e34 -e34],[-1 0 1],Jy-idxb+1,Jy-idxb+1)]/h;
    % Correct the values at the corners
    NDer34(1,Jy-idxb+3) = 2*NDer34(1,Jy-idxb+3);
    if (BCtype == 0)
        NDer34(end,end) = -1/h;
        NDer34(end,end-1) = 0;
    elseif (BCtype == 1)
        NDer34(end,end) = NDer34(end,end)-pt(idxa);
        NDer34(end,end-1) = 2*NDer34(end,end-1);
    end    
    % Normal derivative operator at the interface Gamma41
    e41 = ones(idxa,1);
    NDer41 = [speye(idxa) -spdiags(0.5*[-e41 (4+eta*h^2)*e41 -e41],[-1 0 1],idxa,idxa)]/h;
    % Correct the values at the corners
    NDer41(end,end-1) = 2*NDer41(end,end-1);
    if (BCtype == 0)
        NDer41(1,idxa+1) = -1/h;
        NDer41(1,idxa+2) = 0;
    elseif (BCtype == 1)
        NDer41(1,idxa+1) = NDer41(1,idxa+1)-pl(idxb);
        NDer41(1,idxa+2) = 2*NDer41(1,idxa+2);
    end
    
    % Fixed-point loop
    for i=1:Nit        
        
        % Dirichlet step (subdomains 1 and 3)
        % -----------------------------------
        if (BCtype == 0)
            % Ensure that the gij are compatible with the other BC, then
            % solve for Omega1
            g12(1) = gb1(end);
            g41(1) = gl1(end);
            %g12(end) = g41(end);
            g12(end) = 0.5*(g12(end)+g41(end));
            g41(end) = g12(end);
            u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,yb,gl1,g12,gb1,g41,0,pl1,pr1,pb1,pt1);
            % Ensure that the gij are compatible with the other BC, then
            % solve for Omega3
            g23(end) = gr3(1);
            g34(end) = gt3(1);
            %g23(1) = g34(1);
            g23(1) = 0.5*(g23(1)+g34(1));
            g34(1) = g23(1);
            u3 = Solve2d_EtaMinusDelta(f3,eta,xa,xmax,yb,ymax,g34,gr3,g23,gt3,0,pl3,pr3,pb3,pt3);
        else
            u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,yb,gl1,lv*g12,gb1,lv*g41,1,pl1,lv*ones(idxb,1),pb1,lv*ones(idxa,1));
            u3 = Solve2d_EtaMinusDelta(f3,eta,xa,xmax,yb,ymax,lv*g34,gr3,lv*g23,gt3,1,lv*ones(Jy-idxb+1,1),pr3,lv*ones(Jx-idxa+1,1),pt3);
        end 
        
        % Neumann step (subdomains 2 and 4)
        % ---------------------------------
        u2 = Solve2d_EtaMinusDelta(f2,eta,xa,xmax,ymin,yb,gl2,gr2,gb2,gt2,1,pl2,pr2,pb2,pt2);
        u4 = Solve2d_EtaMinusDelta(f4,eta,xmin,xa,yb,ymax,gl4,gr4,gb4,gt4,1,pl4,pr4,pb4,pt4); 
        
        % Compute the normal derivative of u1 on Gamma12 and Gamma41 (to be used at the next iteration)
        gl2 = NDer12*[u1(:,end-1);u1(:,end)]+0.5*h*[0; f1(2:end,end)];
        gl2(end) = 0.5*gl2(end);
        if (BCtype == 1)
            gl2(1) = gl2(1)+gb(idxa)+0.5*h*f1(1,end);
        end
        gb4 = NDer41*[u1(end-1,:)';u1(end,:)']+0.5*h*[0; f1(end,2:end)'];
        gb4(end) = 0.5*gb4(end);
        if (BCtype == 1)
            gb4(1) = gb4(1)+gl(idxb)+0.5*h*f1(end,1);
        end
        % Compute the normal derivative of u3 on Gamma23 and Gamma34
        gt2 = NDer23*[u3(2,:)';u3(1,:)']+0.5*h*[f3(1,1:end-1)'; 0];
        gt2(1) = 0.5*gt2(1);
        if (BCtype == 1)
            gt2(end) = gt2(end)+gr(idxb)+0.5*h*f3(1,end);
        end
        gr4 = NDer34*[u3(:,2);u3(:,1)]+0.5*h*[f3(1:end-1,1); 0];
        gr4(1) = 0.5*gr4(1);
        if (BCtype == 1)
            gr4(end) = gr4(end)+gt(idxa)+0.5*h*f3(end,1);
        end
        
        % Update the gij
        % --------------
        g12 = th*u2(:,1)+(1-th)*g12;
        g23 = th*u2(end,:)'+(1-th)*g23;
        g34 = th*u4(:,end)+(1-th)*g34;
        g41 = th*u4(1,:)'+(1-th)*g41;
        
        % Build the solution everywhere u
        u = [[u1(1:end-1,1:end-1) u2(1:end-1,:)]; [u4(1,:) u2(end,2:end)]; [u4(2:end,:) u3(2:end,2:end)]];
        
        % Compute the error
        err1 = -u1+uex(1:idxb,1:idxa);
        err2 = -u2+uex(1:idxb,idxa:end);
        err3 = -u3+uex(idxb:end,idxa:end);
        err4 = -u4+uex(idxb:end,1:idxa);
            
        % Plot the recombined solution (or the error)
        if (plot_iterates) 
            figure('Name','u after Neumann step');
            set(gcf, 'Color', 'w');
            mesh(x1,y1,u1); hold on; mesh(x2,y2,u2); mesh(x3,y3,u3); mesh(x4,y4,u4); hold off;
            xlabel('x'); ylabel('y'); zlabel('u');
%             figure('Name','error after Neumann step');
%             set(gcf, 'Color', 'w');
%             mesh(x1,y1,err1); hold on; mesh(x2,y2,err2); mesh(x3,y3,err3); mesh(x4,y4,err4); hold off;
%             xlabel('x'); ylabel('y'); zlabel('err');
            pause
        end
        
        % Compute the L2 norm and the broken H1 norm
        ErrL2(i) = NormL2_FD_2d(h,u-uex)/NormL2_FD_2d(h,uex);
        ErrH1(i) = (NormH1_FD_2d(h,-u1+uex(1:idxb,1:idxa))+NormH1_FD_2d(h,-u2+uex(1:idxb,idxa:end))+NormH1_FD_2d(h,-u3+uex(idxb:end,idxa:end))+NormH1_FD_2d(h,-u4+uex(idxb:end,1:idxa)))/NormH1_FD_2d(h,uex);
        Slope(i) = abs(1-2*th)^i;       
        E22_0(i) = abs(err2(end,1));
        
%         % Print the value of the jump across the cross-point
%         if (i==1)
%            fprintf('Jump at iteration 1: e_2^1(0,0) = delta = %d \n', err2(end,1)); 
%            fprintf('                     e_4^1(0,0) = %d \n', err4(1,end));                 
%         elseif (i==2)
%            fprintf('Jump at iteration 2: e_2^2(0,0) = %d \n', err2(end,1));                 
%            fprintf('                     e_4^2(0,0) = %d \n', err4(1,end));                 
%         end
        
    end
    
    % Plot convergence history
    if (plot_cvg_curve)
        figure('Name','Convergence history');
        set(gcf, 'Color', 'w');
        semilogy(1:1:Nit,ErrL2,'-'); hold on;
        semilogy(1:1:Nit,ErrH1,'-');
        %semilogy(1:1:Nit,E22_0,'-');
        xlabel('iteration k'); ylabel('err');
        legend(' || err(k) ||_{L^2}',' || err(k) ||_{H^1}');
        %legend(' || err(k) ||_{L^2}',' || err(k) ||_{H^1}',' | e_2^k(0,0) |');
    end
    
    % Compute slope of linear regression for ErrL2 and ErrH1
    pL2 = polyfit(2:Nit,log(ErrL2(2:Nit))',1);
    pH1 = polyfit(2:Nit,log(ErrH1(2:Nit))',1);
    pE22 = polyfit(2:Nit,log(E22_0(2:Nit))',1);
    fprintf('Slope for ErrL2 = %d \n', pL2(1)); 
    fprintf('Slope for ErrH1 = %d \n', pH1(1));
    fprintf('Slope for E22_0 = %d \n', pE22(1));

end