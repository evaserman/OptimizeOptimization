function [F,U,freedofs,xPhys_on,x_passive,E0] = MBB(nelx, nely, rmin)

E0 = 100;
NumForces   = 11;


% Multiple (NumForces) Forces in Bottom Right
AddForce = [];
for j = 1:NumForces
    AddForce = [AddForce; (2*(nelx+1)*(nely+1)-(floor(nelx/40)-5)*2*(nely+1)*0-2*(rmin-1)) - 2*(nely+1)*(j-1), 1, 1];
end
F = sparse(AddForce(:,1), AddForce(:,2), AddForce(:,3), 2*(nely+1)*(nelx+1),1);

U = zeros(2*(nely+1)*(nelx+1),1);

fixeddofs = union([(2-1):2:(2*(nely+1)-(2-1))],[2*rmin]);
% fixeddofs = union([(2*rmin-1):2:(2*(nely+1)-(2*rmin-1))],[2*rmin]);
% fixeddofs = union([(2*rmin-1):2:(2*(nely+1)-(2*rmin-1))],[2*rmin:2*(nely+1):2*(nely+1)*1+2*rmin]);
save('fixeddofs.mat','fixeddofs')

alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

%% DESIGNATE PASSIVE DESIGN VARIABLES
x_passive = zeros(nely,nelx); 
for i = 1:rmin-1
     x_passive(i,:) = 1;
end

for i = nely-rmin+2:nely 
      x_passive(i,:) = 1; 
end

for i = nelx-rmin+2:nelx 
      x_passive(:,i) = 1; 
end

%% DESIGNATE ACTIVE ELEMENTS
%Create base supports, add elements parallel to x axis

xPhys_on = zeros(nely,nelx);

%turn on feet at suports
xPhys_on(nely-rmin*2:nely,nelx-round(nelx/20):nelx) = 1;

% Turn on entire base
% xPhys_on(nely-rmin*2:nely,1:nelx) = 1;


%% Visualize Fixed Dofs locations and force
visualize(nelx, nely, F, fixeddofs, xPhys_on)



