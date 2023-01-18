function [f,eta,x,y,xmin,xmax,ymin,ymax,h,Jx,Jy,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,uex] = InitializeData2d_EtaMinusDelta(Test)
% Initializes the data required to run Script2d_EtaMinusDelta
% Test 0: f=1 and u=0 on the boudary, Omega=(-1,1)x(-1,1)
% Test 1: f=0, eta=0 with Dirichlet BC 0,0,0,1, Omega=(0,1)x(0,1)
% Test 2: uex = sin(pi*x)sin(pi*y) on the unit square

    % Test 0: f=1 and homogeneous Dirichlet BC on the unit square
    % -----------------------------------------------------------
    if (Test == 0)
        % Geometrical data
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        Jx = 101;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Choice of eta
        eta = 0;
        % Source term
        %f = ones(Jy,Jx);
        f = x+y;
        idxa = 1+round(0.5*(Jx-1));
        idxb = 1+round(0.5*(Jy-1));
        f(1:idxb-1,1:idxa-1) = f(1:idxb-1,1:idxa-1)+sin(2*atan(y(1:idxb-1,1:idxa-1)./x(1:idxb-1,1:idxa-1)));
        f(idxb+1:end,idxa+1:end) = f(idxb+1:end,idxa+1:end)-sin(2*atan(y(idxb+1:end,idxa+1:end)./x(idxb+1:end,idxa+1:end)));
        % Exact solution (unknown)
        uex = false;
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 0;
        % Dirichlet case
        if (BCtype == 0)
            % Left side
            pl = zeros(Jy,1);
            gl = zeros(Jy,1);
            % Right side
            pr = zeros(Jy,1);
            gr = zeros(Jy,1);
            % Bottom side
            pb = zeros(Jx,1);
            gb = zeros(Jx,1);
            % Top side
            pt = zeros(Jx,1);
            gt = zeros(Jx,1);
        % Robin case
        elseif (BCtype == 1)
            % Left side
            pl = 1e10*ones(Jy,1);
            gl = zeros(Jy,1);
            % Right side
            pr = 1e10*ones(Jy,1);
            gr = zeros(Jy,1);
            % Bottom side
            pb = zeros(Jx,1);
            gb = zeros(Jx,1);
            % Top side
            pt = zeros(Jx,1);
            gt = zeros(Jx,1);
        end
    
    % Test 1: Laplacian with Dirichlet BC 0,0,0,1 on the unit square
    % --------------------------------------------------------------
    elseif (Test == 1)
        % Geometrical data
        xmin = 0;
        xmax = 1;
        ymin = 0;
        ymax = 1;
        Jx = 201;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Physical data
        eta = 0;
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
       
    % Test 2: uex = sin(pi*x)sin(pi*y) on the unit square
    % ----------------------------------------------------
    elseif (Test == 2)
        % Geometrical data
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        Jx = 101;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Physical data
        eta = 0;
        % Exact solution
        %uex = sin(pi*x).*sin(pi*y);
        uex = false;
        % Source term
        f = 2*pi^2*sin(pi*x).*sin(pi*y);
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 1;
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
        elseif (BCtype == 1)
            % Left side
            pl = ones(Jy,1);
            gl = -pi*sin(pi*y(:,1));
            % Right side
            pr = ones(Jy,1);
            gr = -pi*sin(pi*y(:,1));
            % Bottom side
            pb = ones(Jx,1);
            gb = -pi*sin(pi*x(1,:)');
            % Top side
            pt = ones(Jx,1);
            gt = -pi*sin(pi*x(1,:)');           
        end
        
    % Test 3: uex = sin(pi*x).*cos(pi/2*y) on the unit square
    % -----------------------------------------------------------
    elseif (Test == 3)
        % Geometrical data
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        Jx = 101;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Choice of eta
        eta = 1;
        % Source term
        f = 2*ones(Jy,Jx)+5*pi^2/4*sin(pi*x).*cos(pi/2*y);
        % Exact solution (unknown)
        %uex = sin(pi*x).*cos(pi/2*y);
        uex = false;
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 0;
        % Dirichlet case
        if (BCtype == 0)
            % Left side
            pl = zeros(Jy,1);
            gl = zeros(Jy,1);
            % Right side
            pr = zeros(Jy,1);
            gr = zeros(Jy,1);
            % Bottom side
            pb = zeros(Jx,1);
            gb = zeros(Jx,1);
            % Top side
            pt = zeros(Jx,1);
            gt = zeros(Jx,1);
        end

    % Test 4: other test on the unit square
    % -----------------------------------------------------------
    elseif (Test == 4)
        % Geometrical data
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        Jx = 101;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Choice of eta
        eta = 0;
        % Source term
        f = 10*ones(Jy,Jx)-6*x+pi^2*sin(pi*y);
        %f = -6*x.*(y.^2+1)-2*x.^3;
        % Exact solution (unknown)
        uex = false;
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 0;
        % Dirichlet case
        if (BCtype == 0)
%             % Left side
%             pl = zeros(Jy,1);
%             gl = zeros(Jy,1);
%             % Right side
%             pr = zeros(Jy,1);
%             gr = zeros(Jy,1);
%             % Bottom side
%             pb = zeros(Jx,1);
%             gb = zeros(Jx,1);
%             % Top side
%             pt = zeros(Jx,1);
%             gt = zeros(Jx,1);
            % Left side
            pl = zeros(Jy,1);
            gl = -ones(Jy,1)+sin(pi*y(:,1));
            % Right side
            pr = zeros(Jy,1);
            gr = ones(Jy,1)+sin(pi*y(:,1));
            % Bottom side
            pb = zeros(Jx,1);
            gb = x(1,:).^3';
            % Top side
            pt = zeros(Jx,1);
            gt = x(1,:).^3';
%             % Left side
%             pl = zeros(Jy,1);
%             gl = -(y(:,1).^2+ones(Jy,1));
%             % Right side
%             pr = zeros(Jy,1);
%             gr = y(:,1).^2+ones(Jy,1);
%             % Bottom side
%             pb = zeros(Jx,1);
%             gb = 2*x(1,:).^3';
%             % Top side
%             pt = zeros(Jx,1);
%             gt = 2*x(1,:).^3';
        elseif (BCtype == 1)
            % Left side
            pl = ones(Jy,1);
            gl = y(:,1).^2;
            % Right side
            pr = ones(Jy,1);
            gr = y(:,1).^2;
            % Bottom side
            pb = ones(Jx,1);
            gb = pb;
            % Top side
            pt = ones(Jx,1);
            gt = pt;
%             % Left side
%             pl = 1e14*ones(Jy,1);
%             gl = zeros(Jy,1);
%             % Right side
%             pr = 1e14*ones(Jy,1);
%             gr = zeros(Jy,1);
%             % Bottom side
%             pb = 1e14*ones(Jx,1);
%             gb = zeros(Jx,1);
%             % Top side
%             pt = 1e14*ones(Jx,1);
%             gt = zeros(Jx,1);
        end
        
    % Test 5: f=1 unit square, zero Dirichlet BC on 3 sides, zero Neumann
    % on right side
    % -------------------------------------------------------------------
    elseif (Test == 5)
        % Geometrical data
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        Jx = 101;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Choice of eta
        eta = 0;
        % Source term
        %f = zeros(Jy,Jx);
        f = -(20*x.^3-12*x).*(y.^2-1)-2*x.*(x.^2-1).^2;
        % Exact solution (unknown)
        uex = false;
        %uex = x.*(x.^2-1).^2.*(y.^2-1);
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 1;
        % Robin case
        if (BCtype == 1)
            % Left side
            pl = 1e10*ones(Jy,1);
            gl = zeros(Jy,1);
            % Right side
            pr = zeros(Jy,1);
            gr = zeros(Jx,1);
            % Bottom side
            pb = 1e10*ones(Jx,1);
            gb = zeros(Jx,1);
            % Top side
            pt = 1e10*ones(Jx,1);
            gt = zeros(Jx,1);
        end

    % Test 6: uex = x^2+xy-y on the unit square
    % -------------------------------------------------------------------
    elseif (Test == 6)
        % Geometrical data
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        Jx = 101;
        Jy = 1+(ymax-ymin)/(xmax-xmin)*(Jx-1);
        % Build the grid
        h = (xmax-xmin)/(Jx-1);
        [x,y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
        % Choice of eta
        eta = 0;
        % Source term
        f = -2*ones(Jy,Jx);
        % Exact solution (unknown)
        uex = x.^2+x.*y-y;
        % Boundary conditions (Dirichlet:0 or Robin:1)
        BCtype = 1;
        % Robin case
        if (BCtype == 1)
            % Left side
            pl = ones(Jy,1);
            gl = 3*ones(Jy,1)-3*y(:,1);
            % Right side
            pr = ones(Jy,1);
            gr = 3*ones(Jy,1)+y(:,1);
            % Bottom side
            pb = ones(Jx,1);
            gb = x(1,:).^2'-2*x(1,:)'+2*ones(Jx,1);
            % Top side
            pt = ones(Jx,1);
            gt = x(1,:).^2'+2*x(1,:)'-2*ones(Jx,1);
        end
        
    end
    
end