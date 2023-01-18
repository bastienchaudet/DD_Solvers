function res = NormH1_FD_2d(h,u)
% Computes the L2 norm (in the finite difference sense) of u written in
% matrix form
        
    % Initialization    
    res = 0;
    % Loop on the entries of u
    [Jy,Jx] = size(u);
    for i=2:Jy-1
        res = res+0.25*norm(u(i+1,:)-u(i-1,:))^2;
    end
    for j=2:Jx-1
        res = res+0.25*norm(u(:,j+1)-u(:,j-1))^2;
    end
    res = sqrt(res)+NormL2_FD_2d(h,u);

end