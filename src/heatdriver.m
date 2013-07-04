cd 'C:\Users\david\Documents\GitHub\RBF-heat\src\';
alpha = 0.1;  % Parameter for the heat equation
eps = 1.5 ;  % Shape paramater for the RBF kernel
surfeps=1.75;% Shape parameter for the visualization interpolation
N = 20;     % Somehow related to the number of centers.  For the ME points,
            % the number of centers is (N+1)^2.
M = 100;      % how many iterations to run the simulation for
h = 1/(1.5*N);  % timestep

X = getMEPoints(N);
X = sortrows(X,3);

% Generate some initial data.
U0 = X(:,1).^2 + X(:,3).^2;
for i = 1:size(U0,1)
   if X(i,3) > 0
       U0(i) = 0.8;
   else
       U0(i) = 0;
   end
end

% Try some other initial data
nn=10;
lambda = -nn*(nn+1);
U1 = dsph(nn,X(:,1),X(:,2),X(:,3));
U2 = dsph(7,X(:,1),X(:,2),X(:,3));
U3 = dsph(13,X(:,1),X(:,2),X(:,3));
U0 = 1*U1(:,8) + 0*U2(:,1) + 0*U3(:,4);
U=U0;
% Run the sim
[V,D] = heateqn(U0, X, h, M, alpha, lambda, eps, surfeps);

plot(err)


[A,V,D] = eigenLap(U0, X, h, M, alpha, eps, surfeps);
eigU = A*V(:,2);

lambda = D(2,2);
U0 = eigU;

% Visualize the output
[Xp, Yp, Zp] = sphere((N+1)^2);
W = scatteredScalarMake(X(:,1), X(:,2), X(:,3), eigU, Xp, Yp, Zp, surfeps);

htop = surf(Xp,Yp,Zp,W);
set(htop, 'edgecolor','none')
set(gca,'Visible','off')
daspect([1 1 1]);
drawnow
