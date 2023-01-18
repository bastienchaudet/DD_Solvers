function u = DNSeq3d_TwoSubdomains_EtaMinusDelta(f,eta,x,y,z,gl,gri,gb,gt,gre,gf,BCtype,pl,pri,pb,pt,pre,pf,uex,a,th,Nit,plot_iterates,plot_cvg_curve)

    % Geometric data
    xmin = x(1,1,1);
    xmax = x(1,end,1);
    ymin = y(1,1,1);
    ymax = y(end,1,1);
    zmin = z(1,1,1);
    zmax = z(1,1,end);   
    [Jy,Jx,Jz] = size(x);
    h = x(1,2,1)-x(1,1,1);
    
    % Corresponding index and coordinate
    idxa = 1+round(a*(Jx-1));
    xa = xmin+h*(idxa-1);
    
    % Large value to emulate a Dirichlet BC
    lv = 1e10;
    
    % Subcomponents related to each subdomain
    x1 = x(:,1:idxa,:);
    x2 = x(:,idxa:end,:);
    y1 = y(:,1:idxa,:);
    y2 = y(:,idxa:end,:);
    z1 = z(:,1:idxa,:);
    z2 = z(:,idxa:end,:);
    f1 = f(:,1:idxa,:);
    f2 = f(:,idxa:end,:);
    
    % Initial guess for u on the interface
    g = zeros(Jy,Jz);
    
    % In the case of Dirichlet BC, we can make a better guess
    if (BCtype == 0) 
        g = Solve2d_EtaMinusDelta(g,0,zmin,zmax,ymin,ymax,gre(:,idxa),gf(:,idxa),gb(idxa,:)',gt(idxa,:)',0,pre(:,idxa),pf(idxa,:),pb(idxa,:)',pt(idxa,:)');
    end
    
    % Boundary conditions for the first (left) subdomain
    pb1 = pb(1:idxa,:);
    gb1 = gb(1:idxa,:);
    pt1 = pt(1:idxa,:);
    gt1 = gt(1:idxa,:);
    pre1 = pre(:,1:idxa);
    gre1 = gre(:,1:idxa);
    pf1 = pf(:,1:idxa);
    gf1 = gf(:,1:idxa);
    
    % Boundary conditions for the second (right) subdomain
    pl2 = zeros(Jy,Jz);
    if (BCtype == 0)
        pb2 = lv*ones(Jx-idxa+1,Jz);
        gb2 = lv*gb(idxa:end,:);
        pt2 = pb2;
        gt2 = lv*gt(idxa:end,:);
        pri = lv*ones(Jy,Jz);
        gri = lv*gri;
        pre2 = lv*ones(Jy,Jx-idxa+1);
        gre2 = lv*gre(:,idxa:end);
        pf2 = pre2;
        gf2 = lv*gf(:,idxa:end);
    else
        pb2 = pb(idxa:end,:);
        gb2 = gb(idxa:end,:);
        pt2 = pt(idxa:end,:);
        gt2 = gt(idxa:end,:);
        pre2 = pre(:,idxa:end);
        gre2 = gre(:,idxa:end);
        pf2 = pf(:,idxa:end);
        gf2 = gf(:,idxa:end);
    end
    
    % Relaxation parameter and number of iterations
    ErrL2 = zeros(Nit,1);
    ErrH1 = zeros(Nit,1);
    
    % Normal derivative operator at the interface
    ey = ones(Jy,1);
    ez = ones(Jz,1);
    Dyy = spdiags([-ey 0.5*(6+eta*h^2)*ey -ey],[-1 0 1],Jy,Jy);
    Dzz = spdiags([-ez 0.5*(6+eta*h^2)*ez -ez],[-1 0 1],Jz,Jz);
    % Case of Dirichlet BC
    if (BCtype == 0)
        % Build Ayz without taking care of boundary indices
        Ayz = kron(speye(Jz),Dyy)+kron(Dzz,speye(Jy));
    % Case of Robin BC    
    elseif (BCtype == 1)       
        % Vertical parts of the boundary
        Iy = speye(Jy);
        % Boundary condition on the rear side
        Hzre = sparse(1,1,2*h,Jz,Jz);
        Pre = sparse(Jy,Jy);
        Pre = spdiags(pre(:,idxa),0,Pre);
        % Boundary condition on the front side
        Hzf = sparse(Jz,Jz,2*h,Jz,Jz);
        Pf = sparse(Jy,Jy);
        Pf = spdiags(pf(:,idxa),0,Pf);
        % Horizontal parts of the boundary
        Iz = speye(Jz);
        % Boundary condition on the bottom side
        Hyb = sparse(1,1,2*h,Jy,Jy);
        Pb = sparse(Jz,Jz);
        Pb = spdiags(pb(idxa,:)',0,Pb);
        % Boundary condition on the top side
        Hyt = sparse(Jy,Jy,2*h,Jy,Jy);
        Pt = sparse(Jz,Jz);
        Pt = spdiags(pt(idxa,:)',0,Pt);         
        % Build Ayz
        Ayz = kron(Iz,Dyy)+kron(Hzre,Pre)+kron(Hzf,Pf)+kron(Dzz,Iy)+kron(Pb,Hyb)+kron(Pt,Hyt);
        % Correct some boundary indices
        for i=1:Jy
           Ayz(i,Jy+i) = 2*Ayz(i,Jy+i);
           Ayz(Jy*(Jz-1)+i,Jy*(Jz-2)+i) = 2*Ayz(Jy*(Jz-1)+i,Jy*(Jz-2)+i);
        end
        for k=1:Jz
           Ayz((k-1)*Jy+1,(k-1)*Jy+2) = 2*Ayz((k-1)*Jy+1,(k-1)*Jy+2);
           Ayz(k*Jy,k*Jy-1) = 2*Ayz(k*Jy,k*Jy-1);
        end
        % Correct corner indices
        % Bottom-rear corner
        Ayz(1,1) = (6+eta*h^2)+2*(pre(1,idxa)+pb(idxa,1))*h;
        % Bottom-front corner
        Ayz(1+(Jz-1)*Jy,1+(Jz-1)*Jy) = (6+eta*h^2)+2*(pf(1,idxa)+pb(idxa,end))*h;
        % Top-rear corner
        Ayz(Jy,Jy) = (6+eta*h^2)+2*(pre(end,idxa)+pt(idxa,1))*h;
        % Top-front corner
        Ayz(end,end) = (6+eta*h^2)+2*(pf(end,idxa)+pt(idxa,end))*h;
    end      
    % Build NDer     
    NDer = [speye(Jy*Jz) -0.5*Ayz]/h;

    % Fixed-point loop
    for i=1:Nit

        % Dirichlet step
        % --------------
        if (BCtype == 0)
            % Ensure that g is compatible with the other BC
            g(1,:) = gb1(end,:);
            g(end,:) = gt1(end,:);
            g(2:end-1,1) = gre1(2:end-1,end);
            g(2:end-1,end) = gf1(2:end-1,end);
            u1 = Solve3d_EtaMinusDelta(f1,eta,xmin,xa,ymin,ymax,zmin,zmax,gl,g,gb1,gt1,gre1,gf1,0,pl,pri,pb1,pt1,pre1,pf1);
        else
            u1 = Solve3d_EtaMinusDelta(f1,eta,xmin,xa,ymin,ymax,zmin,zmax,gl,lv*g,gb1,gt1,gre1,gf1,1,pl,lv*ones(Jy,Jz),pb1,pt1,pre1,pf1);
        end
        % Compute the normal derivative of u1
        u1 = reshape(u1,Jy,idxa,Jz);
        u1end = u1(:,end,:);
        u1endm1 = u1(:,end-1,:);
        gl2 = NDer*[u1endm1(:); u1end(:)];
        gl2 = reshape(gl2,Jy,Jz);
        f1p = permute(f1,[1 3 2]);
        gl2(2:end-1,2:end-1) = gl2(2:end-1,2:end-1)+0.5*h*f1p(2:end-1,2:end-1,end);
        % Correct boundary values in case of Robin BC
        if (BCtype == 1)
            % Rear and front edges
            gl2(2:end-1,1) = gl2(2:end-1,1)+gre(2:end-1,idxa)+0.5*h*f1p(2:end-1,1,end);
            gl2(2:end-1,end) = gl2(2:end-1,end)+gf(2:end-1,idxa)+0.5*h*f1p(2:end-1,end,end);            
            % Bottom and top edges
            gl2(1,2:end-1) = gl2(1,2:end-1)+gb(idxa,2:end-1)+0.5*h*f1p(1,2:end-1,end);
            gl2(end,2:end-1) = gl2(end,2:end-1)+gt(idxa,2:end-1)+0.5*h*f1p(end,2:end-1,end);
            % Corners
            gl2(1,1) = gl2(1,1)+(gre(1,idxa)+gb(idxa,1))+0.5*h*f1p(1,1,end);
            gl2(end,1) = gl2(end,1)+(gre(end,idxa)+gt(idxa,1))+0.5*h*f1p(end,1,end);
            gl2(1,end) = gl2(1,end)+(gf(1,idxa)+gb(idxa,end))+0.5*h*f1p(1,end,end);
            gl2(end,end) = gl2(end,end)+(gf(end,idxa)+gt(idxa,end))+0.5*h*f1p(end,end,end);
        end

        % Neumann step
        % ------------
        u2 = Solve3d_EtaMinusDelta(f2,eta,xa,xmax,ymin,ymax,zmin,zmax,gl2,gri,gb2,gt2,gre2,gf2,1,pl2,pri,pb2,pt2,pre2,pf2);
        
        % Update g
        % --------
        u2 = reshape(u2,Jy,Jx-idxa+1,Jz);
        u2p = permute(u2,[1 3 2]);
        g = th*u2p(:,:,1)+(1-th)*g;
        
        % Build the solution everywhere u
        u1pint = permute(u1(:,1:end-1,:),[1 3 2]);
        up = [reshape(u1pint,Jy,(idxa-1)*Jz) reshape(u2p,Jy,(Jx-idxa+1)*Jz)];
        u = ipermute(reshape(up,Jy,Jx,Jz),[1 3 2]);

        % Compute the L2 norm and the broken H1 norm
        ErrL2(i) = NormL2_FD_3d(h,u-uex)/NormL2_FD_3d(h,uex);
        ErrH1(i) = (NormH1_FD_3d(h,-u1+uex(:,1:idxa,:))+NormH1_FD_3d(h,-u2+uex(:,idxa:end,:)))/NormH1_FD_3d(h,uex);

        % Plot the recombined solution
        if (plot_iterates)
            % Define the planes to display the solution
            xs1 = xmin+(xmax-xmin)/3;
            xs2 = xmin+2*(xmax-xmin)/3;
            ys1 = ymin+(ymax-ymin)/3;
            ys2 = ymin+2*(ymax-ymin)/3;
            zs0 = zmin+(zmax-zmin)/2;
            xs11 = xmin+(xa-xmin)/3;
            xs12 = xmin+2*(xa-xmin)/3;
            xs21 = xa+(xmax-xa)/3;
            xs22 = xa+2*(xmax-xa)/3;
            % Define the slices
%             xslice = [xs1,xs2];
%             yslice = [ys1,ys2];
%             zslice = zs0;
            x1slice = [xmin,xa];
            x2slice = [xa,xmax];
            xslice = [xmin,xmax];
            yslice = [ymin,ymax];
            zslice = [zmin,zmax];
            % The solution
%             u1 = reshape(u1,Jy,idxa,Jz);
%             figure('Name','u1 after Neumann step');
%             set(gcf, 'Color', 'w');
%             slice(x1,y1,z1,u1,x1slice,yslice,zslice);
%             %
%             u2 = reshape(u2,Jy,Jx-idxa+1,Jz);
%             figure('Name','u2 after Neumann step');
%             set(gcf, 'Color', 'w');
%             slice(x2,y2,z2,u2,x2slice,yslice,zslice);
            %
            figure('Name','u after Neumann step');
            set(gcf, 'Color', 'w');
            slice(x,y,z,u,xslice,yslice,zslice);
            pause
        end

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