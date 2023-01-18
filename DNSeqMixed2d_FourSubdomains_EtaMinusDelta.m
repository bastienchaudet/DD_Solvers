function u = DNSeqMixed2d_FourSubdomains_EtaMinusDelta(f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex,a,b,th,Nit,plot_iterates,plot_cvg_curve)

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
    % Initial guess for dudn on the interface
    h12 = zeros(idxb,1);
    h23 = zeros(Jx-idxa+1,1);
    h34 = zeros(Jy-idxb+1,1);
    h41 = zeros(idxa,1);
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
    pt1 = lv*ones(idxa,1);
    pl1 = pl(1:idxb);
    gl1 = gl(1:idxb);
    pr1 = zeros(idxb,1);
    % Boundary conditions for the second (right bottom) subdomain
    pb2 = pb(idxa:end);
    gb2 = gb(idxa:end);
    pt2 = zeros(Jx-idxa+1,1);
    pl2 = lv*ones(idxb,1);
    pr2 = pr(1:idxb);
    gr2 = gr(1:idxb);
    % Boundary conditions for the third (top right) subdomain
    pb3 = lv*ones(Jx-idxa+1,1);
    pt3 = pt(idxa:end);
    gt3 = gt(idxa:end);
    pl3 = zeros(Jy-idxb+1,1);
    gr3 = gr(idxb:end);
    pr3 = pr(idxb:end);
    % Boundary conditions for the fourth (top left) subdomain
    pb4 = zeros(idxa,1);
    pt4 = pt(1:idxa);
    gt4 = gt(1:idxa);
    pl4 = pr(idxb:end);
    gl4 = gl(idxb:end);
    pr4 = lv*ones(Jy-idxb+1,1);
    % Initialization for the norms of the error
    ErrL2 = zeros(Nit,1);
    ErrH1 = zeros(Nit,1);
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
    
    % DEBUG
    u2 = zeros(size(f2));
    u4 = zeros(size(f4));
    
    % Fixed-point loop
    for i=1:Nit
        
        % First step (subdomains 1 and 3)
        % -------------------------------
        if (BCtype == 0)
            % Ensure that the gij are compatible with the other BC
            g41(1) = gl(idxb);
            g23(end) = gr(idxb);
            % Then solve in subdomains 1 and 3
            u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,yb,lv*gl1,h12,lv*gb1,lv*g41,1,lv*ones(idxb,1),pr1,lv*ones(idxa,1),pt1);
            u3 = Solve2d_EtaMinusDelta(f3,eta,xa,xmax,yb,ymax,h34,lv*gr3,lv*g23,lv*gt3,1,pl3,lv*ones(Jy-idxb+1,1),pb3,lv*ones(Jx-idxa+1,1));
        else
            u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,yb,gl1,h12,gb1,lv*g41,1,pl1,pr1,pb1,pt1);
            u3 = Solve2d_EtaMinusDelta(f3,eta,xa,xmax,yb,ymax,h34,gr3,lv*g23,gt3,1,pl3,pr3,pb3,pt3);
        end 

        % Compute the normal derivative of u1 on Gamma41
        du1dn14 = -NDer41*[u1(end-1,:)';u1(end,:)']-0.5*h*[0; f1(end,2:end)'];
        du1dn14(end) = 0.5*du1dn14(end);
        if (BCtype == 1)
            du1dn14(1) = du1dn14(1)-gl(idxb)-0.5*h*f1(end,1);
        end
        % Compute the normal derivative of u3 on Gamma23
        du3dn32 = -NDer23*[u3(2,:)';u3(1,:)']-0.5*h*[f3(1,1:end-1)'; 0];
        du3dn32(1) = 0.5*du3dn32(1);
        if (BCtype == 1)
            du3dn32(end) = du3dn32(end)-gr(idxb)-0.5*h*f3(1,end);
        end
        
        % Udpate the gij and hij
        % ----------------------
        g12 = u1(:,end);
        g34 = u3(:,1);
        h41 = -du1dn14;
        h23 = -du3dn32;

        % Second step (subdomains 2 and 4)
        % --------------------------------
        if (BCtype == 0)
            % Ensure that the gij are compatible with the other BC
            g12(1) = gb(idxa);
            g34(end) = gt(idxa);
            % Then solve in subdomains 2 and 4
            u2 = Solve2d_EtaMinusDelta(f2,eta,xa,xmax,ymin,yb,lv*g12,lv*gr2,lv*gb2,h23,1,pl2,lv*ones(idxb,1),lv*ones(Jx-idxa+1,1),pt2);
            u4 = Solve2d_EtaMinusDelta(f4,eta,xmin,xa,yb,ymax,lv*gl4,lv*g34,h41,lv*gt4,1,lv*ones(Jy-idxb+1,1),pr4,pb4,lv*ones(idxa,1));
        else
            u2 = Solve2d_EtaMinusDelta(f2,eta,xa,xmax,ymin,yb,lv*g12,gr2,gb2,h23,1,pl2,pr2,pb2,pt2);
            u4 = Solve2d_EtaMinusDelta(f4,eta,xmin,xa,yb,ymax,gl4,lv*g34,h41,gt4,1,pl4,pr4,pb4,pt4);
        end
        
        % Compute the normal derivatives of u1 and u2 on Gamma12
        du1dn12 = -NDer12*[u1(:,end-1);u1(:,end)]-0.5*h*[0; f1(2:end,end)];
        du1dn12(end) = 0.5*du1dn12(end);
        du2dn21 = -NDer12*[u2(:,2);u2(:,1)]-0.5*h*[0; f2(2:end,1)];
        du2dn21(end) = 0.5*du2dn21(end);
        if (BCtype == 1)
            du1dn12(1) = du1dn12(1)-gb(idxa)-0.5*h*f1(1,end);
            du2dn21(1) = du2dn21(1)-gb(idxa)-0.5*h*f2(1,1);
        end         
        % Compute the normal derivatives of u3 and u4 on Gamma34
        du3dn34 = -NDer34*[u3(:,2);u3(:,1)]-0.5*h*[f3(1:end-1,1); 0];
        du3dn34(1) = 0.5*du3dn34(1);
        du4dn43 = -NDer34*[u4(:,end-1);u4(:,end)]-0.5*h*[f4(1:end-1,end); 0];
        du4dn43(1) = 0.5*du4dn43(1);
        if (BCtype == 1)
            du3dn34(end) = du3dn34(end)-gt(idxa)-0.5*h*f3(end,1);
            du4dn43(end) = du4dn43(end)-gt(idxa)-0.5*h*f4(end,end);
        end 
        
        % Update the gij and hij
        % ----------------------
        g23 = th*u2(end,:)'+(1-th)*g23;
        g41 = th*u4(1,:)'+(1-th)*g41;
        h12 = (1-th)*du1dn12-th*du2dn21;
        h34 = (1-th)*du3dn34-th*du4dn43;        
        
        % Build the solution everywhere u
        u = [[u1(1:end-1,1:end-1) u2(1:end-1,:)]; [u4(1,:) u2(end,2:end)]; [u4(2:end,:) u3(2:end,2:end)]];
        
        % Compute the error
        err1 = -u1+uex(1:idxb,1:idxa);
        err2 = -u2+uex(1:idxb,idxa:end);
        err3 = -u3+uex(idxb:end,idxa:end);
        err4 = -u4+uex(idxb:end,1:idxa);
            
        % Plot the recombined solution (or the error)
        if (plot_iterates) 
            figure('Name','u after second step');
            set(gcf, 'Color', 'w');
            mesh(x1,y1,u1); hold on; mesh(x2,y2,u2); mesh(x3,y3,u3); mesh(x4,y4,u4); hold off;
            xlabel('x'); ylabel('y'); zlabel('u');
            figure('Name','error after second step');
            set(gcf, 'Color', 'w');
            mesh(x1,y1,err1); hold on; mesh(x2,y2,err2); mesh(x3,y3,err3); mesh(x4,y4,err4); hold off;
            xlabel('x'); ylabel('y'); zlabel('err');
            pause
        end
        
        % Compute the L2 norm and the broken H1 norm
        ErrL2(i) = NormL2_FD_2d(h,u-uex)/NormL2_FD_2d(h,uex);
        ErrH1(i) = (NormH1_FD_2d(h,-u1+uex(1:idxb,1:idxa))+NormH1_FD_2d(h,-u2+uex(1:idxb,idxa:end))+NormH1_FD_2d(h,-u3+uex(idxb:end,idxa:end))+NormH1_FD_2d(h,-u4+uex(idxb:end,1:idxa)))/NormH1_FD_2d(h,uex);
        Slope(i) = abs(1-2*th)^i;       
        
    end
    
    % Plot convergence history
    if (plot_cvg_curve)
        figure('Name','Convergence history');
        set(gcf, 'Color', 'w');
        semilogy(1:1:Nit,ErrL2,'-'); hold on;
        semilogy(1:1:Nit,ErrH1,'-');
        xlabel('iteration k'); ylabel('err');
        legend(' || err(k) ||_{L^2}',' || err(k) ||_{H^1}');
    end
    

end