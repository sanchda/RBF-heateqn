function U = heateqnT(U0, x, h, M, alpha, eps)
% In this method, the differential operators
% are expressed in R3, but restricted to acting on the surface of the
% sphere.  The centers x are assumed to be on the surface of the sphere
% (in a "narrow band" about this surface), thus preserving the action of
% the operators.  These operators (in this case, only the scalar Laplacian)
% have been expressed in terms of differential operators on the surface of
% the sphere (TODO: link to Mathematica notebook or derived .pdf).  The
% differential operators are implemented as RBF differentiation matrices
% following the method described in Fasshauer "RBF Collocation Methods and
% Pseudospectral Methods" (http://www.math.iit.edu/~fass/RBFPS.pdf as of
% 9/22/10).  The kernel used is the Gaussian.

% TODO: lit citation

    % Sense the size of the problem from the initial conditions.  The
    % scalars in U0 are assumed to identify with the R3 vectors in X.  If
    % these are not the same length, function should return error.
    % TODO:  return error
    N = size(U0,1);
    n = (1/2)*(sqrt(1 + 4*N)-1);
    n = round(n);
    
    % Parameter to all the RBF calls.  This parameter will have to be
    % changed if a different RBF kernel is used.
    twoeps2 = 2*eps*eps;
   

    % Build the distance matrix.
    r2 = (repmat(x(:,1),[1 size(x(:,1),1)]) - repmat(x(:,1)',[size(x(:,1),1) 1])).^2;
    r2 = r2 + (repmat(x(:,2),[1 size(x(:,2),1)]) - repmat(x(:,2)',[size(x(:,2),1) 1])).^2;
    r2 = r2 + (repmat(x(:,3),[1 size(x(:,3),1)]) - repmat(x(:,3)',[size(x(:,3),1) 1])).^2;
    r2 = eps*eps*r2;
    A = exp(-r2);
    
    % Initialize the differentiation matrices
    % ASSERT:  X contains only points on the surface of the unit sphere
    
    % In order to write things like (x1-y1)*A (see accompanying Mathematica
    % notebook), we need to form some distance matrices.  The following
    % method was inspired by Joseph Kirk, author of distmat.
    % (http://www.mathworks.com/matlabcentral/fileexchange/15145-distance-m
    % atrix/content/distmat.m case 2)

    % TODO: time whether a similar method would yield more efficient
    % initialization for A.
    xdiff1 = reshape(x(:,1),1,N,1);
    xdiff2 = reshape(x(:,1),N,1,1);
    
    ydiff1 = reshape(x(:,2),1,N,1);
    ydiff2 = reshape(x(:,2),N,1,1);
    
    zdiff1 = reshape(x(:,3),1,N,1);
    zdiff2 = reshape(x(:,3),N,1,1);
    
    
    xdiff = -xdiff1(ones(N,1),:,:) + xdiff2(:,ones(N,1),:);
    ydiff = -ydiff1(ones(N,1),:,:) + ydiff2(:,ones(N,1),:);
    zdiff = -zdiff1(ones(N,1),:,:) + zdiff2(:,ones(N,1),:);
    
    % To form Lx, a matrix differential operator, the same process as
    % before is used, but where the RBF kernel is the derivative (in x)
    % of the one before.  Since this function is the exponential, we can
    % reuse some of the previous data in the derivative (this intermediate
    % matrix should really be called something like Bx, but we reuse Lx
    % here).  Finally, Lx = Bx*inv(A), but inv() isn't stable so we do a
    % right-division.
   
    Lx = -twoeps2*A.*xdiff;
    Lx = Lx/A;

    Ly = -twoeps2*A.*ydiff;  	
    Ly = Ly/A;
 	
    Lz = -twoeps2*A.*zdiff;  	
    Lz = Lz/A;
    
    
    Lxy = twoeps2*twoeps2*A.*xdiff.*ydiff;
    Lxy = Lxy/A;
  	
    Lxz = twoeps2*twoeps2*A.*xdiff.*zdiff;
    Lxz = Lxz/A;
    
    Lyz = twoeps2*twoeps2*A.*ydiff.*zdiff;
    Lyz = Lyz/A;

    
    Lxx = twoeps2*A.*(-1 + twoeps2*(xdiff.*xdiff)); 	
    Lxx = Lxx/A;	  	
    
    Lyy = twoeps2*A.*(-1 + twoeps2*(ydiff.*ydiff));  	
    Lyy = Lyy/A;

    Lzz = twoeps2*A.*(-1 + twoeps2*(zdiff.*zdiff));
    Lzz = Lzz/A;

    
    % Not useful any longer
    clear('xdiff1','xdiff2','xdiff');
    clear('ydiff1','ydiff2','ydiff');
    clear('zdiff1','zdiff2','zdiff');
    

    % Initialize the scalar Laplacian.
    % ASSERT:  X contains only points on the surface of the unit sphere
    %
    % On the surface of the sphere, the
    % scalar Laplacian has the following form:
    %
    % lapu = - 2x	 (dU/dx)   - 2y    (dU/dy)   - 2z    (dU/dz)
    %        - 2xy   (ddU/dxy) - 2xz   (ddU/dxz) - 2yz   (ddU/dyz)
    %        + (xx-1)(ddU/dxx) + (yy-1)(ddU/dyy) + (zz-1)(ddU/dzz)
    %
    % We can write this as a matrix since we're using matrix operators and
    % matrix multiplication is associative.  However, multiplication by the
    % coordinate should scale the proper element of the resulting vector
    % (i.e., -x*(dU/dx) means the ith element of dU/dx is multiplied by the
    % x coordinate of the ith RBF center.
    
    X = diag(x(:,1));
    Y = diag(x(:,2));
    Z = diag(x(:,3));
  
    lapu = - 2*X*Lx       - 2*Y*Ly       - 2*Z*Lz ...
           - 2*X*Y*Lxy    - 2*X*Z*Lxz    - 2*Y*Z*Lyz ...
           + (eye(N)-X*X)*Lxx  + (eye(N)-Y*Y)*Lyy  + (eye(N)-Z*Z)*Lzz;

    % Clear up the unused matrices
    clear('X','Y','Z');
    clear('Lx','Ly','Lz');
    clear('Lxy','Lxz','Lyz');
    clear('Lxx','Lyy','Lzz');
 
	% Get the first timestep
    U = U0;
    gmin = min(min(U));
    gmax = max(max(U));
   
 
    [xp, yp, zp] = sphere((n+1)^2);
	% Now commence timestepping.  This is the standard RK4 method, which
	% may or may not be well-suited for this problem.
    for i=1:M
        disp(i);
%         rk1 = h*rkfun(U, lapu);
%         rk2 = h*rkfun(U + (1/2)*rk1, lapu);
%         rk3 = h*rkfun(U + (1/2)*rk2, lapu);
%         rk4 = h*rkfun(U + rk3, lapu);
%         
%         U = U + (1/6)*(rk1 + 2*rk2 + 2*rk3 + rk4);
        U = U + alpha*h*lapu*U;
        % Visualize the output
        W = scatteredScalarMake(x(:,1), x(:,2), x(:,3), U, xp, yp, zp, 0.75);

        figure('Color',[1 1 1]);
        htop = surf(xp,yp,zp,W);
        daspect([1 1 1]);
        %caxis([gmin gmax]);
        set(htop, 'edgecolor','none')
        set(gca,'Visible','off')
        set(gcf,'PaperPositionMode','auto')
        
        name = sprintf('../output2/imageUNSCALED_%03d.png',i-1);
        print('-dpng',name);
        close all force;
        close all hidden;
    end


    % TODO: visualization, analysis, etc.  Currently on the final timestep
    % is visualized from the driver.

end