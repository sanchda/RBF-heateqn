function [A,V,D] = heateqn(U0, x, h, M, alpha, eps, surfeps)
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
   

    % Build the distance matrix.  Simplification is due to the centers
    % being on the surface of the sphere.
    
    r2 = 2*(1 - x(:,1)*x(:,1).' - x(:,2)*x(:,2).' - x(:,3)*x(:,3).');
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

    Ly = -twoeps2*A.*ydiff;  	
 	
    Lz = -twoeps2*A.*zdiff;  	
    
    
    Lxy = twoeps2*twoeps2*A.*xdiff.*ydiff;
  	
    Lxz = twoeps2*twoeps2*A.*xdiff.*zdiff;
    
    Lyz = twoeps2*twoeps2*A.*ydiff.*zdiff;

    
    Lxx = twoeps2*A.*(-1 + twoeps2*(xdiff.*xdiff)); 		  	
    
    Lyy = twoeps2*A.*(-1 + twoeps2*(ydiff.*ydiff));  	

    Lzz = twoeps2*A.*(-1 + twoeps2*(zdiff.*zdiff));

    
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
    
    [V, D] = eig(lapu);
    
%     for i = 1:size(V,1)
%        eigU = A*V(:,i);
%        W = scatteredScalarMake(x(:,1), x(:,2), x(:,3), eigU, xp, yp, zp, surfeps);
% 
%        figure('Color',[1 1 1]);
%        htop = surf(xp,yp,zp,W);
%        daspect([1 1 1]);
%        %caxis([gmin gmax]);
%        set(htop, 'edgecolor','none')
%        set(gca,'Visible','off')
%        set(gcf,'PaperPositionMode','auto')
%         
%        name = sprintf('../lapout2/image_%03d_%33d.png',i-1,D(i,i));
%        print('-dpng',name);
%        close all force;
%        close all hidden;
%        
%     end


end