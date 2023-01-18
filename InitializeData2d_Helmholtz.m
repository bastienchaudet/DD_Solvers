function [f,k,x,y,xmin,xmax,ymin,ymax,h,Jx,Jy,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex] = InitializeData2d_Helmholtz(Test)
% Initializes the data required to run Script2d_Helmholtz
% Test 0: uex = exp(i(kx*x+ky*y))
% Test 1: uex = cos(kx*x)sin(ky*y)
% Test 2: Laplacian with Dirichlet BC 0,0,0,1
% Test 3: uex = sin(pi*x)sin(k*pi*y)
% Test 4: waveguide (0,2)x(0,1) with punctual source
% Test 5: free space (0,1)x(0,1) with punctual source
% Test 6: uex = (x(x-1)y(y-1))^2

    % Test 0: uex = exp(i(kx*x+ky*y)) on the unit square
    % --------------------------------------------------
    if (Test == 0)
        % Geometrical data
        xmin = 0;
        xmax = 1;
        ymin = 0;
        ymax = 1;
        Jx = 200;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Physical data
        n = 1;
        k = n*pi*2/sqrt(2);
        theta = pi/4;
        kx = k*cos(theta);
        ky = k*sin(theta);
        % Exact solution
        uex = exp(1i*(kx*x+ky*y));
        % Source term
        f = zeros(Jy,Jx);
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 0;
        % Dirichlet case
        if (BCtype == 0)
            % Left side
            pl = zeros(Jy,1);
            gl = exp(1i*ky*y(:,1));
            % Right side
            pr = zeros(Jy,1);
            gr = exp(1i*(kx+ky*y(:,1)));
            % Bottom side
            pb = zeros(Jx,1);
            gb = exp(1i*kx*x(1,:).');
            % Top side
            pt = zeros(Jx,1);
            gt = exp(1i*(kx*x(1,:).'+ky));
        % Robin case    
        elseif (BCtype == 1)
            BCtype = 1;
            % Left side
            pl = 1i*k*ones(Jy,1);
            gl = 1i*(k-kx)*exp(1i*ky*y(:,1));
            % Right side
            pr = 1i*k*ones(Jy,1);
            gr = 1i*(k+kx)*exp(1i*(kx+ky*y(:,1)));
            % Bottom side
            pb = 1i*k*ones(Jx,1);
            gb = 1i*(k-ky)*exp(1i*kx*x(1,:).');
            % Top side
            pt = 1i*k*ones(Jx,1);
            gt = 1i*(k+ky)*exp(1i*(kx*x(1,:).'+ky));
        end
        
    % Test 1: uex = cos(kx*x)sin(ky*y) on the unit square
    % ---------------------------------------------------
    elseif (Test == 1)
        % Geometrical data
        xmin = 0;
        xmax = 1;
        ymin = 0;
        ymax = 1;
        Jx = 200;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Physical data
        n = 1;
        k = n*pi*2/sqrt(2);
        theta = pi/4;
        kx = k*cos(theta);
        ky = k*sin(theta);
        % Exact solution
        uex = cos(kx*x).*sin(ky*y);
        % Source term
        f = zeros(Jy,Jx);
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 1;
        % Dirichlet case
        if (BCtype == 0)
            % Left side
            pl = zeros(Jy,1);
            gl = sin(ky*y(:,1));
            % Right side
            pr = zeros(Jy,1);
            gr = (-1)^n*sin(ky*y(:,1));
            % Bottom side
            pb = zeros(Jx,1);
            gb = pb;
            % Top side
            pt = zeros(Jx,1);
            gt = pt;
        % Robin case    
        elseif (BCtype == 1)
        end
        
    % Test 2: Laplacian with Dirichlet BC 0,0,0,1 on the unit square
    % --------------------------------------------------------------
    elseif (Test == 2)
        % Geometrical data
        xmin = 0;
        xmax = 1;
        ymin = 0;
        ymax = 1;
        Jx = 200;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Physical data
        k = 0;
        % Exact solution (unknown)
        uex = false;
        % Source term
        f = zeros(Jy,Jx);
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 0;
        % Dirichlet case
        if (BCtype == 0)
            % Left side
            pl = zeros(Jy,1);
            gl = pl;
            gl(end) = 1;
            % Right side
            pr = zeros(Jy,1);
            gr = pr;
            gr(end) = 1;
            % Bottom side
            pb = zeros(Jx,1);
            gb = pb;
            % Top side
            pt = zeros(Jx,1);
            gt = ones(Jx,1);
        end
        
    % Test 3: uex = sin(pi*x)sin(k*pi*y) on the unit square
    % -----------------------------------------------------
    elseif (Test == 3)
        % Geometrical data
        xmin = 0;
        xmax = 1;
        ymin = 0;
        ymax = 1;
        Jx = 100;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Physical data
        k = 10;
        % Exact solution
        %uex = sin(pi*x).*sin(k*pi*y);
        uex = false;
        % Source term
        f = (k^2-pi^2-k^2*pi^2)*sin(pi*x).*sin(k*pi*y);
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 0;
        % Dirichlet case
        if (BCtype == 0)
            % Left side
            pl = zeros(Jy,1);
            gl = pl;
            % Right side
            pr = zeros(Jy,1);
            gr = pr;
            % Bottom side
            pb = zeros(Jx,1);
            gb = pb;
            % Top side
            pt = zeros(Jx,1);
            gt = pt;
        end
        
    % Test 4: waveguide (0,2)x(0,1) with punctual source term  
    % -------------------------------------------------------
    elseif (Test == 4)
        % Geometrical data
        xmin = 0;
        xmax = 2;
        ymin = 0;
        ymax = 1;
        Jx = 301;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Physical data
        k = 20;
        % Exact solution (unknown)
        uex = false;
        % Source term
        a = 0.01;
        b = 0.1;
        r = sqrt((x-1).^2+(y-0.5).^2);
        fab = (2*r.^3-3*(a+b)*r.^2+6*a*b*r+b^2*(b-3*a))/(b-a)^3;
        f = zeros(Jy,Jx);
        f(r<=a) = 10;
        f((a<r)&(r<b)) = 10*fab((a<r)&(r<b));
%         mesh(x,y,f);
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 1;
        % Robin case    
        if (BCtype == 1)
            % Left side
            pl = 1i*k*ones(Jy,1);
            gl = zeros(Jy,1);
            % Right side
            pr = 1i*k*ones(Jy,1);
            gr = zeros(Jy,1);
            % Bottom side
            pb = 1e12*ones(Jx,1);
            gb = zeros(Jx,1);
            % Top side
            pt = 1e12*ones(Jx,1);
            gt = zeros(Jx,1);
        end
    
    % Test 5: free space (0,1)x(0,1) with punctual source term  
    % --------------------------------------------------------
    elseif (Test == 5)
        % Geometrical data
        xmin = 0;
        xmax = 1;
        ymin = 0;
        ymax = 1;
        % Physical data
        k = 30;
        Nppw = 20;
        % Build the grid
        Jx = 1+2*floor((1+Nppw*k/(2*pi)*(xmax-xmin))/2)
        %Jx = 201;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Exact solution (unknown)
        uex = false;
        % Source term
        a = 0.01;
        b = 0.1;
        r = sqrt((x-0.5).^2+(y-0.5).^2);
        fab = (2*r.^3-3*(a+b)*r.^2+6*a*b*r+b^2*(b-3*a))/(b-a)^3;
        f = zeros(Jy,Jx);
        f(r<=a) = 10;
        f((a<r)&(r<b)) = 10*fab((a<r)&(r<b));
%        mesh(x,y,f);
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 1;
        % Robin case    
        if (BCtype == 1)
            % Left side
            pl = 1i*k*ones(Jy,1);
            gl = zeros(Jy,1);
            % Right side
            pr = 1i*k*ones(Jy,1);
            gr = zeros(Jy,1);
            % Bottom side
            pb = 1i*k*ones(Jx,1);
            gb = zeros(Jx,1);
            % Top side
            pt = 1i*k*ones(Jx,1);
            gt = zeros(Jx,1);
        end
   
    % Test 6: uex = (x(x-1)y(y-1))^2  
    % -------------------------------
    elseif (Test == 6)
        % Geometrical data
        xmin = 0;
        xmax = 1;
        ymin = 0;
        ymax = 1;
        % Physical data
        k = 80;
        Nppw = 10;
        % Build the grid
        Jx = 1+floor(Nppw*k/(2*pi)*(xmax-xmin));
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Exact solution
        uex = (x.*(x-1).*y.*(y-1)).^2;
        % Source term
        f = 2*(3*x.^2-3*x+1).*(y.*(y-1)).^2+2*(3*y.^2-3*y+1).*(x.*(x-1)).^2+k^2*uex;
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 1;
        % Robin case    
        if (BCtype == 1)
            % Left side
            pl = 1i*k*ones(Jy,1);
            gl = zeros(Jy,1);
            % Right side
            pr = 1i*k*ones(Jy,1);
            gr = zeros(Jy,1);
            % Bottom side
            pb = 1i*k*ones(Jx,1);
            gb = zeros(Jx,1);
            % Top side
            pt = 1i*k*ones(Jx,1);
            gt = zeros(Jx,1);
        end
    end


end