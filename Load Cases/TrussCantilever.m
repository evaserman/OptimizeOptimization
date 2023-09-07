
function [F,freedofs,num_nod,num_ele,L,C,S,ele_nod,nod_coor,...
    ntoud,ntolr,szt,ndof_c,ele_dof_t,mNodes] =...
    TrussCantilever(nbglr,nbgtb,nelx,nely,cScale,cover)

%Generate szt, ntolr, ntoud based on number of baes and nelx, nely
ntoud = nbglr+1;
ntolr = nbgtb+1;
% szt = min(floor(nely/nbgtb),floor(nelx/nbglr));
szt = cScale;
ndof_c = (nelx+1)*(nely+1)*2;

% element nodes
ele_nod=[];
for i = 1:ntoud-1
    for j = 1:ntolr-1   
        % connect current node to node below
        ele_nod=[ele_nod; j+(i-1)*ntolr,  j+(i-1)*ntolr+1];
        % connect current node to node to the right
        ele_nod=[ele_nod; j+(i-1)*ntolr,  j+(i)*ntolr];
        % connect current node to node down and right
        ele_nod=[ele_nod; j+(i-1)*ntolr,  j+(i)*ntolr+1];
        % connect node below to node to the right
        ele_nod=[ele_nod; j+(i-1)*ntolr+1,  j+(i)*ntolr];
    end
end

% add bottom horz trusses
for i = 1:ntoud-1
    % connect current node to node to the right
    ele_nod=[ele_nod; (i)*ntolr,  (i+1)*ntolr];
end

% add far-right vert trusses (spatially, not politically)
for i = 1:ntolr-1
    % connect current node to node below
    ele_nod=[ele_nod;  (ntoud-1)*ntolr+i,  (ntoud-1)*ntolr+i+1];
end

%number of elements
num_ele=size(ele_nod,1);

%number of nodes
num_nod=length(unique(ele_nod));

% establish node parameters
% nodes coordinates
nod_coor=[];
ysize = szt*(ntolr-1);

for i = 1:ntoud
    for j = 1:ntolr  
        % construct node locations, considering node numbers are built from
        % top to bottom, but using up=positive convention in y direction
        nod_coor=[nod_coor; szt*(i-1),  ysize-szt*(j-1)+cover];
    end
end

ele_dof_1 = ele_nod(:,1)*2-1;
ele_dof_2 = ele_nod(:,1)*2;
ele_dof_3 = ele_nod(:,2)*2-1;
ele_dof_4 = ele_nod(:,2)*2;

ele_dof_t = [ele_dof_1 ele_dof_2 ele_dof_3 ele_dof_4];

%applied loads at DOFs
F = zeros(ndof_c,1);
mag = -10;
% assign point load, mod is used to always assign to even value to ensure y direction load
% F( 2*(nely+1)*(nelx-0)+nely+2+ mod(nely, 2)-0 ) = -10;

% fixed nodes
fixedNodes1 = 1:ntolr;
fixedNodes2 = ntolr*ntoud-floor(ntolr/2);

if mod(ntolr,2) == 0  % even number of bars
   % F( 2*(nely+1)*(nelx-0) + (cover+1)*2 + szt*(ntolr/2-1)*2 ) = mag/2;
   % F( 2*(nely+1)*(nelx-0) + (cover+1)*2 + szt*(ntolr/2+0)*2 ) = mag/2;
   F( 2*szt*(nbglr)*(nely+1) + (cover+1)*2 + szt*(ntolr/2-1)*2 ) = mag/2;
   F( 2*szt*(nbglr)*(nely+1) + (cover+1)*2 + szt*(ntolr/2+0)*2 ) = mag/2;
   fixedNodes2 = [fixedNodes2 ntolr*ntoud-ntolr/2+1];
else                  % odd number of bars
    F( 2*(nely+1)*(nelx-0)+nely+2 ) = mag;
end

%Boundary conditions
fixeddofs = (1:(nely+1)*2)';
alldofs = [1:ndof_c]';
freedofs = setdiff(alldofs,fixeddofs);

% fixedNodes2 = [];

mNodes = setdiff(1:num_nod,unique([fixedNodes1 fixedNodes2]));

L = zeros(num_ele,1);
C = zeros(num_ele,1);
S = zeros(num_ele,1);
% FEA Parameters
for e=1:num_ele
 L(e)=sqrt((nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))^2+...
      (nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))^2);
 C(e)=(nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))/L(e);
 S(e)=(nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))/L(e);
end

end










