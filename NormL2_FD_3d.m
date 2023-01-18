function res = NormL2_FD_3d(h,u)
% Computes the L2 norm (in the finite difference sense) of u written in
% matrix form
        
    % Initialization    
    res = 0;
    % Loop on the entries of u
    [Jy,Jx,Jz] = size(u);
    for i=1:Jy
        for j=1:Jx
            for k=1:Jz
                res = res+h^3*abs(u(i,j,k))^2;
            end
        end
    end
    res = sqrt(res);

end