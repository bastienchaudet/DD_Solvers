function [u,nb_solves_out] = NewDNBis2d_EtaMinusDelta(nb_level,f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,nb_solves)
% Domain decomposition into two subdomains which is applied recursively
% nb_level times.
% The first domain is divided vertically into two  subdomains, then each
% of these subdomains is divided horizontally into two subdomains, etc...

    % Keep track of the total number of solves
    nb_solves_out = nb_solves;
    % Geometric data
    a = 0.5;
    th = 0.5;
    Nit = 2;
    if (abs(norm(pl)-norm(pr)) > 1e6)
       Nit = 4;
    end
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
    % Initial guess for u on the interface
    g = zeros(Jy,1);
    % In the case of Dirichlet BC, we can make a better guess
    if (BCtype == 0)
        g = gb(idxa)*ones(Jy,1)+(gt(idxa)-gb(idxa))*(0:1/(Jy-1):1)';
    else
        if(norm(pb) > lv/1e3)
            if(norm(pt) > lv/1e3)
                g = 1/lv*(gb(idxa)*ones(Jy,1)+(gt(idxa)-gb(idxa))*(0:1/(Jy-1):1)');
            else
                g = g(end)+(gb(idxa)/lv-g(end))*(1:-1/(Jy-1):0)';
            end
        elseif(norm(pt) > lv/1e3)
            g = g(1)+(gt(idxa)/lv-g(1))*(0:1/(Jy-1):1)';
        end
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
            g(1) = gb1(end);
            g(end) = gt1(end);
            if (nb_level > 1)
                [u1,nb_solves_out1] = NewDNBis2d_EtaMinusDelta(nb_level-1,f1',eta,y1',x1',gb1,gt1,gl,g,0,pb1,pt1,pl,pr,nb_solves_out);
                u1 = u1';
                nb_solves_out = nb_solves_out1;
            else
                u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,ymax,gl,g,gb1,gt1,0,pl,pr,pb1,pt1);
                nb_solves_out = nb_solves_out+1;
            end
        else
            if (nb_level > 1)
                [u1,nb_solves_out1] = NewDNBis2d_EtaMinusDelta(nb_level-1,f1',eta,y1',x1',gb1,gt1,gl,lv*g,1,pb1,pt1,pl,lv*ones(Jy,1),nb_solves_out);
                u1 = u1';
                nb_solves_out = nb_solves_out1;
            else
                u1 = Solve2d_EtaMinusDelta(f1,eta,xmin,xa,ymin,ymax,gl,lv*g,gb1,gt1,1,pl,lv*ones(Jy,1),pb1,pt1);
                nb_solves_out = nb_solves_out+1;
            end
        end
        % Compute the normal derivative of u1
        gl2 = NDer*[u1(:,end-1);u1(:,end)]+0.5*h*[0; f1(2:end-1,end); 0];
        if (BCtype == 1)
            gl2(1) = gl2(1)+gb(idxa)+0.5*h*f1(1,end);
            gl2(end) = gl2(end)+gt(idxa)+0.5*h*f1(end,end);
        end

        % Neumann step
        % ------------
        % `
        if (BCtype == 0)
            if (nb_level > 1)
                [u2,nb_solves_out2] = NewDNBis2d_EtaMinusDelta(nb_level-1,f2',eta,y2',x2',lv*gb2,lv*gt2,gl2,lv*gr,1,lv*ones(Jx-idxa+1,1),lv*ones(Jx-idxa+1,1),pl2,lv*ones(Jy,1),nb_solves_out);
                u2 = u2';
                nb_solves_out = nb_solves_out2;
            else
                u2 = Solve2d_EtaMinusDelta(f2,eta,xa,xmax,ymin,ymax,gl2,lv*gr,lv*gb2,lv*gt2,1,pl2,lv*ones(Jy,1),lv*ones(Jx-idxa+1,1),lv*ones(Jx-idxa+1,1));
                nb_solves_out = nb_solves_out+1;
            end
        else
            if (nb_level > 1)
                [u2,nb_solves_out2] = NewDNBis2d_EtaMinusDelta(nb_level-1,f2',eta,y2',x2',gb2,gt2,gl2,gr,1,pb2,pt2,pl2,pr,nb_solves_out);
                u2 = u2';
                nb_solves_out = nb_solves_out2;
            else
                u2 = Solve2d_EtaMinusDelta(f2,eta,xa,xmax,ymin,ymax,gl2,gr,gb2,gt2,1,pl2,pr,pb2,pt2);
                nb_solves_out = nb_solves_out+1;
            end
        end
        
        % Update g
        % --------
        g = th*u2(:,1)+(1-th)*g;
        
        % Build the solution everywhere u
        u = [u1(:,1:end-1) u2];       
        
        % Plot the recombined solution
        if (false)%(nb_level == 1)
            figure('Name','Recombined solution u');
            mesh(x1,y1,u1); hold on; mesh(x2,y2,u2); hold off;
            xlabel('x'); ylabel('y'); zlabel('u2');
            pause
        end

    end
    
end
