function u = DNOdd3d_FourSubdomains_EtaMinusDelta(f,eta,x,y,z,gl,gb,gre,gf,BCtype,pl,pb,pre,pf)
% Dirichlet-Neumann method (mixed version) optimized for the case of odd symmetric data
    
    % Geometric data
    a = 0.5;
    b = 0.5;
    th = 0.5;
    Nit = 2;
    xmin = x(1,1,1);
    xmax = x(1,end,1);
    ymin = y(1,1,1);
    zmin = z(1,1,1);
    zmax = z(1,1,end);
    [Jy,Jx,Jz] = size(x);
    h = x(1,2,1)-x(1,1,1);
    
    % Corresponding indices and coordinates
    idxa = 1+round(a*(Jx-1));
    idyb = 1+round(b*(Jy-1));
    xa = xmin+h*(idxa-1);
    yb = ymin+h*(idyb-1);
    
    % Large value to emulate a Dirichlet BC
    lv = 1e10;
    
    % Subcomponents related to subdomains 1 and 2
    x1 = x(1:idyb,1:idxa,:);
    x2 = x(1:idyb,idxa:end,:);
    y1 = y(1:idyb,1:idxa,:);
    y2 = y(1:idyb,idxa:end,:);
    z1 = z(1:idyb,1:idxa,:);
    z2 = z(1:idyb,idxa:end,:);
    f1 = f(1:idyb,1:idxa,:);
    f2 = f(1:idyb,idxa:end,:);
    
    % Initial guess for u on the interface
    h12 = zeros(idyb,Jz);
    g41 = zeros(idxa,Jz);
    
    % Build gr, pr by symmetry
    gri = gl(end:-1:1,:);
    pri = pl(end:-1:1,:);
    
    % In the case of Dirichlet BC, we can make a better guess
    if (BCtype == 0)
        gce = gre(idyb,idxa)*ones(Jz,1)+(gf(idyb,idxa)-gre(idyb,idxa))*(0:1/(Jz-1):1)';
        g41 = Solve2d_EtaMinusDelta(g41,0,zmin,zmax,xmin,xa,gre(idyb,1:idxa)',gf(idyb,1:idxa)',gl(idyb,:)',gce,0,pre(idyb,1:idxa)',pf(idyb,1:idxa)',pl(idyb,:)',pl(idyb,:)');
    end  
    
    % Boundary conditions for the first (left bottom) subdomain
    gb1 = gb(1:idxa,:);
    gl1 = gl(1:idyb,:);
    gre1 = gre(1:idyb,1:idxa);
    gf1 = gf(1:idyb,1:idxa); 
    if (BCtype == 0)
        pb1 = lv*ones(idxa,Jz);
        pl1 = lv*ones(idyb,Jz);
        pre1 = lv*ones(idyb,idxa);
        pf1 = pre1;
    elseif (BCtype == 1)
        pb1 = pb(1:idxa,:);
        pl1 = pl(1:idyb,:);
        pre1 = pre(1:idyb,1:idxa);
        pf1 = pf(1:idyb,1:idxa);
    end
    pt1 = lv*ones(idxa,Jz);
    pr1 = zeros(idyb,Jz);   
    
    % Boundary conditions for the second (right bottom) subdomain
    gb2 = gb(idxa:end,:);
    gri2 = gri(1:idyb,:);
    gre2 = gre(1:idyb,idxa:end);
    gf2 = gf(1:idyb,idxa:end);
    if (BCtype == 0)
        pb2 = lv*ones(Jx-idxa+1,Jz);
        pri2 = lv*ones(Jy-idyb+1,Jz);
        pre2 = lv*ones(Jy-idyb+1,Jx-idxa+1);
        pf2 = pre2;       
    elseif (BCtype == 1)
        pb2 = pb(idxa:end,:);
        pri2 = pri(1:idyb,:);
        pre2 = pre(1:idyb,idxa:end);
        pf2 = pf(1:idyb,idxa:end);
    end
    pl2 = lv*ones(idyb,Jz);
    pt2 = zeros(Jx-idxa+1,Jz);
    
    % Normal derivative operator at the interface Gamma12 (from subdomain 2)
    e12 = ones(idyb,1);
    ez = ones(Jz,1);
    Dyy = spdiags([-e12 0.5*(6+eta*h^2)*e12 -e12],[-1 0 1],idyb,idyb);
    Dzz = spdiags([-ez 0.5*(6+eta*h^2)*ez -ez],[-1 0 1],Jz,Jz);
    Iy = speye(idyb);
    Iz = speye(Jz);
    % Case of Dirichlet BC
    if (BCtype == 0)
        % Build Ayz without taking care of boundary indices
        Ayz = kron(Iz,Dyy)+kron(Dzz,Iy);
    % Case of Robin BC    
    elseif (BCtype == 1)       
        % Boundary condition on the rear side
        Hzre = sparse(1,1,2*h,Jz,Jz);
        Prey = sparse(idyb,idyb);
        Prey = spdiags(pre(1:idyb,idxa),0,Prey);
        % Boundary condition on the front side
        Hzf = sparse(Jz,Jz,2*h,Jz,Jz);
        Pfy = sparse(idyb,idyb);
        Pfy = spdiags(pf(1:idyb,idxa),0,Pfy);
        % Boundary condition on the bottom side
        Hyb = sparse(1,1,2*h,idyb,idyb);
        Pbz = sparse(Jz,Jz);
        Pbz = spdiags(pb(idxa,:)',0,Pbz);       
        % Build Ayz
        Ayz = kron(Iz,Dyy)+kron(Hzre,Prey)+kron(Hzf,Pfy)+kron(Dzz,Iy)+kron(Pbz,Hyb);
        % Correct some boundary indices
        % Rear and front edges 
        for i=1:idyb
           Ayz(i,idyb+i) = 2*Ayz(i,idyb+i);
           Ayz(idyb*(Jz-1)+i,idyb*(Jz-2)+i) = 2*Ayz(idyb*(Jz-1)+i,idyb*(Jz-2)+i);
        end
        % Bottom edge
        for k=1:Jz
           Ayz((k-1)*idyb+1,(k-1)*idyb+2) = 2*Ayz((k-1)*idyb+1,(k-1)*idyb+2);
        end
        % Correct corner indices
        % Bottom-rear corner
        Ayz(1,1) = (6+eta*h^2)+2*(pre(1,idxa)+pb(idxa,1))*h;
        % Bottom-front corner
        Ayz(1+(Jz-1)*idyb,1+(Jz-1)*idyb) = (6+eta*h^2)+2*(pf(1,idxa)+pb(idxa,end))*h;
        % Top-rear corner
        Ayz(idyb,idyb) = (6+eta*h^2)+2*pre(idyb,idxa)*h;
        % Top-front corner
        Ayz(Jz*idyb,Jz*idyb) = (6+eta*h^2)+2*pf(idyb,idxa)*h;
    end   
%     % Take into account Neumann BC on the top side of subdomain 2
%     % Top edge
%     for k=1:Jz
%         Ayz(k*idyb,k*idyb-1) = 2*Ayz(k*idyb,k*idyb-1);
%     end
        
    % Normal derivative operator at the interface Gamma41 (from subdomain 1)
    e41 = ones(idxa,1);
    Dxx = spdiags([-e41 0.5*(6+eta*h^2)*e41 -e41],[-1 0 1],idxa,idxa);
    Ix = speye(idxa);
    % Case of Dirichlet BC
    if (BCtype == 0)
        % Build Axz without taking care of boundary indices
        Axz = kron(Iz,Dxx)+kron(Dzz,Ix);
    % Case of Robin BC    
    elseif (BCtype == 1)       
        % Boundary condition on the rear side
        Hzre = sparse(1,1,2*h,Jz,Jz);
        Prex = sparse(idxa,idxa);
        Prex = spdiags(pre(idyb,1:idxa)',0,Prex);
        % Boundary condition on the front side
        Hzf = sparse(Jz,Jz,2*h,Jz,Jz);
        Pfx = sparse(idxa,idxa);
        Pfx = spdiags(pf(idyb,1:idxa)',0,Pfx);
        % Boundary condition on the left side
        Hxl = sparse(1,1,2*h,idxa,idxa);
        Plz = sparse(Jz,Jz);
        Plz = spdiags(pl(idyb,:)',0,Plz);       
        % Build Axz
        Axz = kron(Iz,Dxx)+kron(Hzre,Prex)+kron(Hzf,Pfx)+kron(Dzz,Ix)+kron(Plz,Hxl);
        % Correct some boundary indices
        % Rear and front edges
        for j=1:idxa
           Axz(j,idxa+j) = 2*Axz(j,idxa+j);
           Axz(idxa*(Jz-1)+j,idxa*(Jz-2)+j) = 2*Axz(idxa*(Jz-1)+j,idxa*(Jz-2)+j);
        end
        % Left edge
        for k=1:Jz
           Axz((k-1)*idxa+1,(k-1)*idxa+2) = 2*Axz((k-1)*idxa+1,(k-1)*idxa+2);
        end
        % Correct corner indices
        % Left-rear corner
        Axz(1,1) = (6+eta*h^2)+2*(pre(idyb,1)+pl(idyb,1))*h;
        % Left-front corner
        Axz(1+(Jz-1)*idxa,1+(Jz-1)*idxa) = (6+eta*h^2)+2*(pf(idyb,1)+pl(idyb,end))*h;
        % Right-rear corner
        Axz(idxa,idxa) = (6+eta*h^2)+2*pre(idyb,idxa)*h;
        % Right-front corner
        Axz(Jz*idxa,Jz*idxa) = (6+eta*h^2)+2*pf(idyb,idxa)*h;
    end 
%     % Take into account Neumann BC on the right side of subdomain 1
%     % Right edge
%     for k=1:Jz
%         Axz(k*idxa,k*idxa-1) = 2*Axz(k*idxa,k*idxa-1);
%     end
    
    % Correct the values of Ayz and Axz at the cross-edge
    for k=1:Jz
        Ayz(k*idyb,k*idyb-1) = 2*Ayz(k*idyb,k*idyb-1);
        Axz(k*idxa,k*idxa-1) = 2*Axz(k*idxa,k*idxa-1);
    end
    
    % Build NDer12 and NDer41     
    NDer12 = [speye(idyb*Jz) -0.5*Ayz]/h;
    NDer41 = [speye(idxa*Jz) -0.5*Axz]/h;    
    
    % Fixed-point loop
    for i=1:Nit
        
        % First step (subdomain 1)
        % ------------------------
        if (BCtype == 0)
            % Ensure that the gij are compatible with the other BC, then
            % solve for Omega1
            g41(1,2:end-1) = gl1(end,2:end-1);
            g41(:,1) = gre1(end,:)';
            g41(:,end) = gf1(end,:)';
            % Solve
            u1 = Solve3d_EtaMinusDelta(f1,eta,xmin,xa,ymin,yb,zmin,zmax,lv*gl1,h12,lv*gb1,lv*g41,lv*gre1,lv*gf1,1,pl1,pr1,pb1,pt1,pre1,pf1);
        else
            u1 = Solve3d_EtaMinusDelta(f1,eta,xmin,xa,ymin,yb,zmin,zmax,gl1,h12,gb1,lv*g41,gre1,gf1,1,pl1,pr1,pb1,pt1,pre1,pf1);
        end        
                   
        % Compute the normal derivative of u1 on Gamma41       
        u1endy = u1(end,:,:);
        u1endym1 = u1(end-1,:,:);
        du1dn14 = NDer41*[u1endym1(:); u1endy(:)];
        du1dn14 = reshape(du1dn14,idxa,Jz);
        f1py = permute(f1,[2 3 1]);
        du1dn14(2:end,:) = du1dn14(2:end,:)+0.5*h*f1py(2:end,:,end);
        % Correct boundary values
%        du1dn14(end,:) = 0.5*du1dn14(end,:);
        % Right edge
        du1dn14(end,:) = du1dn14(end,:)+h12(end,:);
        if (BCtype == 1)
            % Rear and front edges          
            du1dn14(2:end-1,1) = du1dn14(2:end-1,1)+gre1(end,2:end-1)';
            du1dn14(2:end-1,end) = du1dn14(2:end-1,end)+gf1(end,2:end-1)';  
            % Left edge
            du1dn14(1,2:end-1) = du1dn14(1,2:end-1)+gl1(end,2:end-1)+0.5*h*f1py(1,2:end-1,end);
            % Corners
            du1dn14(1,1) = du1dn14(1,1)+(gre1(end,1)+gl1(end,1))+0.5*h*f1py(1,1,end);
            du1dn14(1,end) = du1dn14(1,end)+(gf1(end,1)+gl1(end,end))+0.5*h*f1py(1,end,end);
            %
            du1dn14(end,1) = du1dn14(end,1)+gre1(end,end);
            du1dn14(end,end) = du1dn14(end,end)+gf1(end,end);
        end
 
        % Deduce the normal derivative of u3 on Gamma23 by (odd) symmetry
        du3dn32 = -du1dn14(end:-1:1,:);
        
        % Udpate the gij and hij
        % ----------------------
        u1 = reshape(u1,idyb,idxa,Jz);
        u1px = permute(u1,[1 3 2]);
        g12 = u1px(:,:,end);
        h23 = -du3dn32;
        
        % Second step (subdomain 2)
        % -------------------------
        if (BCtype == 0)
            % Ensure that the gij are compatible with the other BC, then
            % solve for Omega1
            g12(1,2:end-1) = gb1(end,2:end-1);
            g12(:,1) = gre1(:,end);
            g12(:,end) = gf1(:,end);
            % Solve
            u2 = Solve3d_EtaMinusDelta(f2,eta,xa,xmax,ymin,yb,zmin,zmax,lv*g12,lv*gri2,lv*gb2,h23,lv*gre2,lv*gf2,1,pl2,pri2,pb2,pt2,pre2,pf2);
        elseif (BCtype == 1)
            u2 = Solve3d_EtaMinusDelta(f2,eta,xa,xmax,ymin,yb,zmin,zmax,g12,gri2,gb2,h23,gre2,gf2,1,pl2,pri2,pb2,pt2,pre2,pf2);
        end
        
        % Compute the normal derivative of u2 on Gamma12
        u2 = reshape(u2,idyb,Jx-idxa+1,Jz);
        u21x = u2(:,1,:);
        u22x = u2(:,2,:);
        du2dn21 = NDer12*[u22x(:); u21x(:)];
        du2dn21 = reshape(du2dn21,idyb,Jz);
        f2px = permute(f2,[1 3 2]);
        du2dn21(2:end,:) = du2dn21(2:end,:)+0.5*h*f2px(2:end,:,1);
        % Correct boundary values
%         du2dn21(end,:) = 0.5*du2dn21(end,:);
        % Top edge
        du2dn21(end,:) = du2dn21(end,:)+h23(end,:);
        if (BCtype == 1)
            % Rear and front edges                      
            du2dn21(2:end-1,1) = du2dn21(2:end-1,1)+gre2(2:end-1,1);
            du2dn21(2:end-1,end) = du2dn21(2:end-1,end)+gf2(2:end-1,1);
            % Bottom edge
            du2dn21(1,2:end-1) = du2dn21(1,2:end-1)+gb2(1,2:end-1)+0.5*h*f2px(1,2:end-1,1);
            % Corners
            du2dn21(1,1) = du2dn21(1,1)+(gre2(1,1)+gb2(1,1))+0.5*h*f2px(1,1,1);
            du2dn21(1,end) = du2dn21(1,end)+(gf2(1,1)+gb2(1,end))+0.5*h*f2px(1,end,1);
            %
            du2dn21(end,1) = du2dn21(end,1)+gre2(end,1);
            du2dn21(end,end) = du2dn21(end,end)+gf2(end,1);
        end   
        
        du2dn21
        
        % Plot the results
        % ----------------
        if (true)
            % Define the planes to display the solution
            zs0 = zmin+(zmax-zmin)/2;
            xs11 = xmin+(xa-xmin)/3;
            xs12 = xmin+2*(xa-xmin)/3;
            xs21 = xa+(xmax-xa)/3;
            xs22 = xa+2*(xmax-xa)/3;
            ys1 = ymin+(yb-ymin)/3;
            ys2 = ymin+2*(yb-ymin)/3;
            % Define the slices
            x1slice = [xmin,xa];
            x2slice = [xa,xmax];
            yslice = [ymin,yb];
            zslice = [zmin,zmax];
%             x1slice = [xs11,xs12];
%             x2slice = [xs21,xs22];
%             yslice = [ys1,ys2];
%             zslice = [zs0];
            % The solution
            figure('Name','u1 after Neumann step');
            set(gcf, 'Color', 'w');
            slice(x1,y1,z1,u1,x1slice,yslice,zslice);
            %
            figure('Name','u2 after Neumann step');
            set(gcf, 'Color', 'w');
            slice(x2,y2,z2,u2,x2slice,yslice,zslice);
            pause                         
        end
        
        % Update the gij and hij
        % ----------------------
        u2py = permute(u2,[2 3 1]); 
        g41 = (1-th)*g41-th*u2py(end:-1:1,:,end);
        h12 = (1-th)*h12-th*du2dn21;
        
        % Build the solution everywhere u by symmetry
        % In subdomains 1 and 2
        u2px = permute(u2,[1 3 2]);
        u1pxint = permute(u1(:,1:end-1,:),[1 3 2]);
        u12px = [reshape(u1pxint,idyb,(idxa-1)*Jz) reshape(u2px,idyb,(Jx-idxa+1)*Jz)];
        u12 = ipermute(reshape(u12px,idyb,Jx,Jz),[1 3 2]);
        % In the whole domain
        u12pyint = permute(u12(1:end-1,:,:),[2 3 1]);
        u34py = -permute(u12(end:-1:1,end:-1:1,:),[2 3 1]);
        upy = [reshape(u12pyint,Jx,(idyb-1)*Jz) reshape(u34py,Jx,(Jy-idyb+1)*Jz)];        
        u = ipermute(reshape(upy,Jx,Jz,Jy),[2 3 1]);
                
    end

end