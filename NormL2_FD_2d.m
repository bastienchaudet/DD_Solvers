function res = NormL2_FD_2d(h,u)
% Computes the L2 norm (in the finite difference sense) of u written in
% matrix form
        
    % Initialization    
    res = 0;
    % Loop on the entries of u
    [Jy,Jx] = size(u);
    for i=1:Jy
        for j=1:Jx
            res = res+h^2*abs(u(i,j))^2;
        end
    end
    res = sqrt(res);

end