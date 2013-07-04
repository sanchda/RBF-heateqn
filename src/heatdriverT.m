cd 'C:\Users\david\Documents\GitHub\RBF-heat\src\';
alpha = 0.1;  % Parameter for the heat equation
eps = 0.7;  % Shape paramater for the RBF kernel
surfeps=0.75;% Shape parameter for the visualization interpolation
M = 100;      % how many iterations to run the simulation for
N = 14;
h = 1/(2*N);  % timestep

X = getTpoints();
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
U1 = dsph(3,X(:,1),X(:,2),X(:,3));
U0 = U1(:,2);
U=U0;
% Run the sim
U = heateqnT(U0, X, h, M, alpha, eps);


% Visualize the output
[Xp, Yp, Zp] = torus(1,235,2);
if size(U,2) == 3
    W1 = scatteredScalarMake(X(:,1), X(:,2), X(:,3), U(:,1), Xp, Yp, Zp, surfeps);
    W2 = scatteredScalarMake(X(:,1), X(:,2), X(:,3), U(:,2), Xp, Yp, Zp, surfeps);
    W3 = scatteredScalarMake(X(:,1), X(:,2), X(:,3), U(:,3), Xp, Yp, Zp, surfeps);

    W = zeros(size(W1,1),size(W1,1),3);
    W(:,:,1) = (W1-min(min(W1)))/max(max(W1));
    W(:,:,2) = (W2-min(min(W2)))/max(max(W2));
    W(:,:,3) = (W3-min(min(W2)))/max(max(W3));
    W = W/max(max(max(W)));
elseif size(U,2) == 1
    W = scatteredScalarMake(X(:,1), X(:,2), X(:,3), U, Xp, Yp, Zp, surfeps);
end

htop = surf(Xp,Yp,Zp,W);
set(htop, 'edgecolor','none')
set(gca,'Visible','off')
daspect([1 1 1]);
drawnow
