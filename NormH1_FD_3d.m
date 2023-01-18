function res = NormH1_FD_3d(h,u)
% Computes the L2 norm (in the finite difference sense) of u written in
% matrix form
        
    % Initialization    
    res = 0;
    % Loop on the entries of u
    [Jy,Jx,Jz] = size(u);
    for i=2:Jy-1
        duy = u(i+1,:,:)-u(i-1,:,:);
        res = res+0.25*norm(duy(:))^2;
    end
    for j=2:Jx-1
        dux = u(:,j+1,:)-u(:,j-1,:);
        res = res+0.25*norm(dux(:))^2;
    end
    for k=2:Jz-1
        duz = u(:,:,k+1)-u(:,:,k-1);
        res = res+0.25*norm(duz(:))^2;
    end
    res = sqrt(res)+NormL2_FD_3d(h,u);

end