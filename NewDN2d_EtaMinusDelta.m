function [u,u_e,u_o] = NewDN2d_EtaMinusDelta(nb_level,f,eta,x,y,gl,gr,gb,gt,BCtype,pl,pr,pb,pt,plot_evenodd)
% Domain decomposition into four subdomains (using an even/odd decomposition 
% of the data) which is applied recursively nb_level times.
% The first domain is divided into four subdomains, then each
% of these subdomains is divided into four subdomains, etc...

    % Decompose the data into the even/odd symmetric parts
    % Source term
    f_e = 0.5*(f+f(end:-1:1,end:-1:1));
    f_o = 0.5*(f-f(end:-1:1,end:-1:1));
    % Boundary conditions
    gl_e = 0.5*(gl+gr(end:-1:1));
    gl_o = 0.5*(gl-gr(end:-1:1));
    gb_e = 0.5*(gb+gt(end:-1:1));
    gb_o = 0.5*(gb-gt(end:-1:1));       
    pl_e = 0.5*(pl+pr(end:-1:1));
    pl_o = 0.5*(pl-pr(end:-1:1));
    pb_e = 0.5*(pb+pt(end:-1:1));
    pb_o = 0.5*(pb-pt(end:-1:1));    
    
    % Solve the even symmetric part using the DNEven method
    u_e = DNEven2d_FourSubdomains_EtaMinusDelta(nb_level,f_e,eta,x,y,gl_e,gb_e,BCtype,pl_e,pb_e);
%    u_e = DNEven_FourSubdomains_EtaMinusDelta(nb_level,f_e,eta,x,y,gl_e,gb_e,BCtype,pl,pb);

    % Solve the odd symmetric part using the DNOdd method
    u_o = DNOdd2d_FourSubdomains_EtaMinusDelta(nb_level,f_o,eta,x,y,gl_o,gb_o,BCtype,pl_o,pb_o);
%    u_o = DNOdd_FourSubdomains_EtaMinusDelta(nb_level,f_o,eta,x,y,gl_o,gb_o,BCtype,pl,pb);
    
    % Recompose the solution
    u = u_e + u_o;
    
    % Plot the even/odd symmetric parts of the solution
    if (plot_evenodd)
        % Even symmetric part    
        figure('Name','Even symmetric part of the solution');
        set(gcf, 'Color', 'w');
        mesh(x,y,u_e);
        xlabel('x');
        ylabel('y');
        zlabel('u_e');
        % Odd symmetric part
        figure('Name','Odd symmetric part of the solution');
        set(gcf, 'Color', 'w');
        mesh(x,y,u_o);
        xlabel('x');
        ylabel('y');
        zlabel('u_o');
    end
    
end