function u = DNOdd2d_FourSubdomains_EtaMinusDelta(nb_level,f,eta,x,y,gl,gb,BCtype,pl,pb)
% Dirichlet-Neumann method (mixed sequential version) optimized for the case of odd symmetric data

    % Geometric data
    a = 0.5;
    b = 0.5;
    th = 0.5;
    Nit = 2;
    xmin = x(1,1);
    xmax = x(1,end);
    ymin = y(1,1);
    [Jy,Jx] = size(x);
    h = x(1,2)-x(1,1);
    % Corresponding indices and coordinates
    idxa = 1+round(a*(Jx-1));
    idxb = 1+round(b*(Jy-1));
    xa = xmin+h*(idxa-1);
    yb = ymin+h*(idxb-1);
    % Large value to emulate a Dirichlet BC
    lv = 1e10;
    % Subcomponents related to subdomains 1 and 2
    x1 = x(1:idxb,1:idxa);
    x2 = x(1:idxb,idxa:end);
    y1 = y(1:idxb,1:idxa);
    y2 = y(1:idxb,idxa:end);
    % Source term for subdomains 1 and 2
    f1 = f(1:idxb,1:idxa);
    f2 = f(1:idxb,idxa:end);
    % Initial guess for u on the interface
    g41 = zeros(idxa,1);
    % Initial guess for dudn on the interface
    h12 = zeros(idxb,1);
    % Build gr, pr, gt by symmetry
    gr = -gl(end:-1:1);
    pr = -pl(end:-1:1);
    gt = -gb(end:-1:1);
    % In the case of Dirichlet BC, we can make a better guess
    if (BCtype == 0)
        gcp_vert = (1-b)*gb(idxa)+b*gt(idxa);
        gcp_hor = (1-a)*gl(idxb)+a*gr(idxb);
        gcp = 0.5*(gcp_vert+gcp_hor);
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
        
        % First step (subdomain 1)
        % ------------------------
        if (BCtype == 0)
            % Ensure that the gij are compatible with the other BC
            g41(1) = gl(idxb);
            % Then solve in subdomain 1
            if (nb_level > 1)
                u1 = NewDN2d_EtaMinusDelta(nb_level-1,f1,eta,x1,y1,lv*gl1,h12,lv*gb1,lv*g41,1,lv*ones(idxb,1),pr1,lv*ones(idxa,1),pt1,true);
            else
                u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,yb,lv*gl1,h12,lv*gb1,lv*g41,1,lv*ones(idxb,1),pr1,lv*ones(idxa,1),pt1);
            end
        else
            if (nb_level > 1)
                u1 = NewDN2d_EtaMinusDelta(nb_level-1,f1,eta,x1,y1,gl1,h12,gb1,lv*g41,1,pl1,pr1,pb1,pt1,true);
            else
                u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,yb,gl1,h12,gb1,lv*g41,1,pl1,pr1,pb1,pt1);
            end
        end 

        % Compute the normal derivative of u1 on Gamma41
        du1dn14 = -NDer41*[u1(end-1,:)';u1(end,:)']-0.5*h*[0; f1(end,2:end)'];
        du1dn14(end) = 0.5*du1dn14(end);
        if (BCtype == 1)
            du1dn14(1) = du1dn14(1)-gl(idxb)-0.5*h*f1(end,1);
        end
        % Deduce the normal derivative of u3 on Gamma23 by symmetry
        du3dn32 = -du1dn14(end:-1:1);
        
        % Udpate the gij and hij
        % ----------------------
        g12 = u1(:,end);
        h23 = -du3dn32;

        % Second step (subdomain 2)
        % -------------------------
        if (BCtype == 0)
            % Ensure that the gij are compatible with the other BC
            g12(1) = gb(idxa);
            % Then solve in subdomain 2
            if (nb_level > 1)
                u2 = NewDN2d_EtaMinusDelta(nb_level-1,f2,eta,x2,y2,lv*g12,lv*gr2,lv*gb2,h23,1,pl2,lv*ones(idxb,1),lv*ones(Jx-idxa+1,1),pt2,true);
            else
                u2 = Solve2d_EtaMinusDelta(f2,eta,xa,xmax,ymin,yb,lv*g12,lv*gr2,lv*gb2,h23,1,pl2,lv*ones(idxb,1),lv*ones(Jx-idxa+1,1),pt2);
            end
        else
            if (nb_level > 1)
                u2 = NewDN2d_EtaMinusDelta(nb_level-1,f2,eta,x2,y2,lv*g12,gr2,gb2,h23,1,pl2,pr2,pb2,pt2,true);
            else
                u2 = Solve2d_EtaMinusDelta(f2,eta,xa,xmax,ymin,yb,lv*g12,gr2,gb2,h23,1,pl2,pr2,pb2,pt2);
            end
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
        
        % Update the gij and hij
        % ----------------------
        g41 = -th*u2(end,end:-1:1)'+(1-th)*g41;
        h12 = (1-th)*du1dn12-th*du2dn21;
        
        % Build the solution everywhere u by symmetry
        u = [[u1(1:end-1,1:end-1) u2(1:end-1,:)]; [u1(end,:) u2(end,2:end)]; [-u2(end-1:-1:1,end:-1:1) -u1(end-1:-1:1,end-1:-1:1)]];
       
    end

end