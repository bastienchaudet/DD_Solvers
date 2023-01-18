function [f,eta,x,y,z,xmin,xmax,ymin,ymax,zmin,zmax,h,Jx,Jy,Jz,gl,gri,gb,gt,gre,gf,BCtype,pl,pri,pb,pt,pre,pf,uex] = InitializeData3d_EtaMinusDelta(Test)
% Initializes the data required to run Script3d_EtaMinusDelta
% Test 0: f=1 and u=0 on the boudary, Omega=(-1,1)x(-1,1)x(-1,1)

    % Test 0: f=1 and homogeneous Dirichlet BC on the unit cube
    % -----------------------------------------------------------
    if (Test == 0)
        % Geometrical data
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        zmin = -1;
        zmax = 1;
        Jx = 31;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        Jz = 1+(zmax-zmin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y,z] = meshgrid(xmin:h:xmax,ymin:h:ymax,zmin:h:zmax);
        % Choice of eta
        eta = 0;
        % Source term
        f = ones(Jy,Jx,Jz);
        % Exact solution (unknown)
        uex = false;
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 0;
        % Dirichlet case
        if (BCtype == 0)
            % Left side
            pl = zeros(Jy,Jz);
            gl = zeros(Jy,Jz);
            % Right side
            pri = zeros(Jy,Jz);
            gri = zeros(Jy,Jz);
            % Bottom side
            pb = zeros(Jx,Jz);
            gb = zeros(Jx,Jz);
            % Top side
            pt = zeros(Jx,Jz);
            gt = zeros(Jx,Jz);
            % Rear side
            pre = zeros(Jy,Jx);
            gre = zeros(Jy,Jx);
            % Front side
            pf = zeros(Jy,Jx);
            gf = zeros(Jy,Jx);
        % Robin case
        elseif (BCtype == 1)
            % Left side
            pl = ones(Jy,Jz);
            gl = zeros(Jy,Jz);
            % Right side
            pri = ones(Jy,Jz);
            gri = zeros(Jy,Jz);
            % Bottom side
            pb = ones(Jx,Jz);
            gb = zeros(Jx,Jz);
            % Top side
            pt = ones(Jx,Jz);
            gt = zeros(Jx,Jz);
            % Rear side
            pre = ones(Jy,Jx);
            gre = zeros(Jy,Jx);
            % Front side
            pf = ones(Jy,Jx);
            gf = zeros(Jy,Jx);
        end
        
    % Test 1: Laplacian with Dirichlet BC 0,0,0,1 on the unit square
    % --------------------------------------------------------------
    elseif (Test == 1)
        % Geometrical data
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        zmin = -1;
        zmax = 1;
        Jx = 31;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        Jz = 1+(zmax-zmin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y,z] = meshgrid(xmin:h:xmax,ymin:h:ymax,zmin:h:zmax);
        % Physical data
        eta = 0;
        % Exact solution
        uex = false;
        % Source term
        f = zeros(Jy,Jx,Jz);
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 0;
        % Dirichlet case
        if (BCtype == 0)
            % Left side
            pl = zeros(Jy,Jz);
            gl = zeros(Jy,Jz);
            gl(2:end-1,2:end-1) = ones(Jy-2,Jz-2);
            % Right side
            pri = zeros(Jy,Jz);
            gri = zeros(Jy,Jz);
            % Bottom side
            pb = zeros(Jx,Jz);
            gb = zeros(Jx,Jz);
            % Top side
            pt = zeros(Jx,Jz);
            gt = zeros(Jx,Jz);
            % Rear side
            pre = zeros(Jy,Jx);
            gre = zeros(Jy,Jx);
            % Front side
            pf = zeros(Jy,Jx);
            gf = zeros(Jy,Jx);
        end        
        
    % Test 2: uex = sin(pi*x)sin(pi*y)sin(pi*z) on the unit cube
    % ----------------------------------------------------------
    elseif (Test == 2)
        % Geometrical data
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        zmin = -1;
        zmax = 1;
        Jx = 31;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        Jz = 1+(zmax-zmin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y,z] = meshgrid(xmin:h:xmax,ymin:h:ymax,zmin:h:zmax);
        % Physical data
        eta = 0;
        % Exact solution
        %uex = sin(pi*x).*sin(pi*y).*sin(pi*z);
        %uex = sin(pi*x).*cos(0.5*pi*y).*sin(pi*z);
        uex = false;
        % Source term
        %f = 3*pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z);
        f = 5*pi^2/4*sin(pi*x).*cos(0.5*pi*y).*sin(pi*z);
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 0;
        % Dirichlet case
        if (BCtype == 0)
            % Left side
            pl = zeros(Jy,Jz);
            gl = zeros(Jy,Jz);
            % Right side
            pri = zeros(Jy,Jz);
            gri = zeros(Jy,Jz);
            % Bottom side
            pb = zeros(Jx,Jz);
            gb = zeros(Jx,Jz);
            % Top side
            pt = zeros(Jx,Jz);
            gt = zeros(Jx,Jz);
            % Rear side
            pre = zeros(Jy,Jx);
            gre = zeros(Jy,Jx);
            % Front side
            pf = zeros(Jy,Jx);
            gf = zeros(Jy,Jx);
        end
        
    % Test 3: uex = xy-z^2/2 on the unit cube
    % -----------------------------------------------------------
    elseif (Test == 3)
        % Geometrical data
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        zmin = -1;
        zmax = 1;
        Jx = 31;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        Jz = 1+(zmax-zmin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y,z] = meshgrid(xmin:h:xmax,ymin:h:ymax,zmin:h:zmax);
        % Choice of eta
        eta = 0;
        % Source term
        f = ones(Jy,Jx,Jz);
        % Exact solution (unknown)
        uex = x.*y-0.5*z.^2;
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 0;
        % Dirichlet case
        if (BCtype == 0)
            % Left side
            pl = zeros(Jy,Jz);
            yp = permute(y,[1 3 2]);
            zp = permute(z,[1 3 2]);
            gl = -yp(:,:,1)-0.5*zp(:,:,1).^2;
            % Right side
            pri = zeros(Jy,Jz);
            gri = yp(:,:,1)-0.5*zp(:,:,1).^2;
            % Bottom side
            pb = zeros(Jx,Jz);
            xp = permute(x,[2 3 1]);
            zp = permute(z,[2 3 1]); 
            gb = -xp(:,:,1)-0.5*zp(:,:,1).^2;
            % Top side
            pt = zeros(Jx,Jz);
            gt = xp(:,:,1)-0.5*zp(:,:,1).^2;
            % Rear side
            pre = zeros(Jy,Jx);
            gre = x(:,:,1).*y(:,:,1)-0.5*ones(Jy,Jx);
            % Front side
            pf = zeros(Jy,Jx);
            gf = x(:,:,1).*y(:,:,1)-0.5*ones(Jy,Jx);
        % Robin case
        elseif (BCtype == 1)
            % Left side
            pl = ones(Jy,Jz);
            yp = permute(y,[1 3 2]);
            zp = permute(z,[1 3 2]);
            gl = -2*yp(:,:,1)-0.5*zp(:,:,1).^2;
            % Right side
            pri = ones(Jy,Jz);
            gri = 2*yp(:,:,1)-0.5*zp(:,:,1).^2;
            % Bottom side
            pb = ones(Jx,Jz);
            xp = permute(x,[2 3 1]);
            zp = permute(z,[2 3 1]); 
            gb = -2*xp(:,:,1)-0.5*zp(:,:,1).^2;
            % Top side
            pt = ones(Jx,Jz);
            gt = 2*xp(:,:,1)-0.5*zp(:,:,1).^2;
            % Rear side
            pre = ones(Jy,Jx);
            gre = x(:,:,1).*y(:,:,1)-1.5*ones(Jy,Jx);
            % Front side
            pf = ones(Jy,Jx);
            gf = x(:,:,1).*y(:,:,1)-1.5*ones(Jy,Jx);
        end
        
 % Test 4: uex = x+y-z on the unit cube
 % -----------------------------------------------------------
    elseif (Test == 4)
        % Geometrical data
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        zmin = -1;
        zmax = 1;
        Jx = 11;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        Jz = 1+(zmax-zmin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y,z] = meshgrid(xmin:h:xmax,ymin:h:ymax,zmin:h:zmax);
        % Choice of eta
        eta = 0;
        % Source term
        f = zeros(Jy,Jx,Jz);
        % Exact solution (unknown)
        %uex = x+y-z;
        %uex = x+y;
        uex = false;
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 0;
        % Robin case
        if (BCtype == 0)
            % Left side
            pl = zeros(Jy,Jz);
            yp = permute(y,[1 3 2]);
            zp = permute(z,[1 3 2]);
            %gl = -ones(Jy,Jz)+yp(:,:,1)-zp(:,:,1);
            gl = -ones(Jy,Jz)+yp(:,:,1);
            % Right side
            pri = zeros(Jy,Jz);
            %gri = ones(Jy,Jz)+yp(:,:,1)-zp(:,:,1);
            gri = ones(Jy,Jz)+yp(:,:,1);
            % Bottom side
            pb = zeros(Jx,Jz);
            xp = permute(x,[2 3 1]);
            zp = permute(z,[2 3 1]); 
            %gb = -ones(Jx,Jz)+xp(:,:,1)-zp(:,:,1);
            gb = -ones(Jx,Jz)+xp(:,:,1);
            % Top side
            pt = zeros(Jx,Jz);
            %gt = ones(Jx,Jz)+xp(:,:,1)-zp(:,:,1);
            gt = ones(Jx,Jz)+xp(:,:,1);
            % Rear side
            pre = zeros(Jy,Jx);
            %gre = x(:,:,1)+y(:,:,1)+ones(Jy,Jx);
            gre = x(:,:,1)+y(:,:,1);
            % Front side
            pf = zeros(Jy,Jx);
            %gf = x(:,:,1)+y(:,:,1)-ones(Jy,Jx);
            gf = x(:,:,1)+y(:,:,1);            
        elseif (BCtype == 1)
            % Left side
            pl = ones(Jy,Jz);
            yp = permute(y,[1 3 2]);
            zp = permute(z,[1 3 2]);
            %gl = -2*ones(Jy,Jz)+yp(:,:,1)-zp(:,:,1);
            gl = -2*ones(Jy,Jz)+yp(:,:,1);
            % Right side
            pri = ones(Jy,Jz);
            %gri = 2*ones(Jy,Jz)+yp(:,:,1)-zp(:,:,1);
            gri = 2*ones(Jy,Jz)+yp(:,:,1);
            % Bottom side
            pb = ones(Jx,Jz);
            xp = permute(x,[2 3 1]);
            zp = permute(z,[2 3 1]); 
            %gb = -2*ones(Jx,Jz)+xp(:,:,1)-zp(:,:,1);
            gb = -2*ones(Jx,Jz)+xp(:,:,1);
            % Top side
            pt = ones(Jx,Jz);
            %gt = 2*ones(Jx,Jz)+xp(:,:,1)-zp(:,:,1);
            gt = 2*ones(Jx,Jz)+xp(:,:,1);
            % Rear side
            pre = ones(Jy,Jx);
            %gre = x(:,:,1)+y(:,:,1)+2*ones(Jy,Jx);
            gre = x(:,:,1)+y(:,:,1);
            % Front side
            pf = ones(Jy,Jx);
            %gf = x(:,:,1)+y(:,:,1)-2*ones(Jy,Jx);
            gf = x(:,:,1)+y(:,:,1);
        end        
        
    end
    
end