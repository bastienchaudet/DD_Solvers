function u = DNEven2d_FourSubdomains_EtaMinusDelta(nb_level,f,eta,x,y,gl,gb,BCtype,pl,pb)
% Dirichlet-Neumann method (standard sequential version) optimized for the case of even symmetric data
    
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
    g12 = zeros(idxb,1);
    g41 = zeros(idxa,1);
    % Build gr, pr, gt by symmetry
    gr = gl(end:-1:1);
    pr = pl(end:-1:1);
    gt = gb(end:-1:1);
    % In the case of Dirichlet BC, we can make a better guess
    if (BCtype == 0)
        gcp_vert = (1-b)*gb(idxa)+b*gt(idxa);
        gcp_hor = (1-a)*gl(idxb)+a*gr(idxb);
        gcp = 0.5*(gcp_vert+gcp_hor);
        g12 = gb(idxa)*ones(idxb,1)+(gcp-gb(idxa))*(0:1/(idxb-1):1)'; 
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
        
        % Dirichlet step (subdomain 1)
        % ----------------------------
        if (BCtype == 0)
            % Ensure that the gij are compatible with the other BC, then
            % solve for Omega1
            g12(1) = gb1(end);
            g41(1) = gl1(end);
            %g12(end) = g41(end);
            g12(end) = 0.5*(g12(end)+g41(end));
            g41(end) = g12(end);
            if (nb_level > 1)
                u1 = NewDN2d_EtaMinusDelta(nb_level-1,f1,eta,x1,y1,gl1,g12,gb1,g41,0,pl1,pr1,pb1,pt1,true);
            else
                u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,yb,gl1,g12,gb1,g41,0,pl1,pr1,pb1,pt1);
            end
        else
            if (nb_level > 1)
                u1 = NewDN2d_EtaMinusDelta(nb_level-1,f1,eta,x1,y1,gl1,lv*g12,gb1,lv*g41,1,pl1,lv*ones(idxb,1),pb1,lv*ones(idxa,1),true); 
            else
                u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,yb,gl1,lv*g12,gb1,lv*g41,1,pl1,lv*ones(idxb,1),pb1,lv*ones(idxa,1));
            end 
        end
        
        % Compute the normal derivative of u1 on Gamma12 and Gamma41
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
        % Deduce the normal derivative of u3 on Gamma23 by symmetry
        gt2 = gb4(end:-1:1);
        
        % Neumann step (subdomain 2)
        % --------------------------
        if (nb_level > 1)
            u2 = NewDN2d_EtaMinusDelta(nb_level-1,f2,eta,x2,y2,gl2,gr2,gb2,gt2,1,pl2,pr2,pb2,pt2,false);
        else
            u2 = Solve2d_EtaMinusDelta(f2,eta,xa,xmax,ymin,yb,gl2,gr2,gb2,gt2,1,pl2,pr2,pb2,pt2);
        end
        
        % Update the gij
        % --------------
        g12 = th*u2(:,1)+(1-th)*g12;
        g41 = th*u2(end,end:-1:1)'+(1-th)*g41;
        
        % Build the solution everywhere u by symmetry
        u = [[u1(1:end-1,1:end-1) u2(1:end-1,:)]; [u1(end,:) u2(end,2:end)]; [u2(end-1:-1:1,end:-1:1) u1(end-1:-1:1,end-1:-1:1)]];
        
    end

end