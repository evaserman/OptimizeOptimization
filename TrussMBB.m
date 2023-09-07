
function [F,ele_dof_c,freedofs,num_nod,num_ele,L,C,S,ele_nod,nod_coor,...
    ntoud,ntolr,szt,ndof,ele_dof_t,mNodes] =...
    TrussMBB(nbglr,nbgtb,nelx,nely)

%Generate szt, ntolr, ntoud based on number of baes and nelx, nely
ntoud = nbglr+1;
ntolr = nbgtb+1;
szt = min(floor(nely/nbgtb),floor(nelx/nbglr));
ndof = (nelx+1)*(nely+1)*2;

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
mapped_dof = [];
for i = 1:ntoud
    for j = 1:ntolr  
        % construct node locations, considering node numbers are built from
        % top to bottom, but using up=positive convention in y direction
        nod_coor=[nod_coor; szt*(i-1),  ysize-szt*(j-1)];

        % Also make map of how DoFs correspond to nodes
        mapped_dof = [mapped_dof; (j-1)*szt*2 + (i-1)*szt*(szt*2*(ntolr-1)+2)+1,  (j-1)*szt*2 + (i-1)*szt*(szt*2*(ntolr-1)+2)+2 ];
    end
end


% elements degree of freedom (DOF) 
% dofs of each node of each truss. convert node map to dof map, then concatenate
% change this based on szt, nelx, nely
ele_dof_c= [mapped_dof(ele_nod(:,1),:)  mapped_dof(ele_nod(:,2),:)];


ele_dof_1 = ele_nod(:,1)*2-1;
ele_dof_2 = ele_nod(:,1)*2;
ele_dof_3 = ele_nod(:,2)*2-1;
ele_dof_4 = ele_nod(:,2)*2;

ele_dof_t = [ele_dof_1 ele_dof_2 ele_dof_3 ele_dof_4];



ele_dof_conn = zeros(ndof);
% vertical dof maps
for i2 = 1:num_nod
    
    % add node to right
    if nod_coor(i2,1) < szt*(ntoud-1)
        neighbor1 = mapped_dof(i2,1)+(nely+1)*2 ;
        ele_dof_conn(mapped_dof(i2,1), neighbor1) = 1;
    end

    % add node to left
    if nod_coor(i2,1) > 0
        neighbor2 = mapped_dof(i2,1)-(nely+1)*2;
        ele_dof_conn(mapped_dof(i2,1), neighbor2) = 1;
    end

    % add node above
    if nod_coor(i2,2) < szt*(ntolr-1)
        neighbor3 = mapped_dof(i2,2)-2;
        ele_dof_conn(mapped_dof(i2,2), neighbor3) = 1;
    end

    %add node below
    if nod_coor(i2,2) > 0
        neighbor4 = mapped_dof(i2,2)+2;
        ele_dof_conn(mapped_dof(i2,2), neighbor4) = 1;
    end

end

%applied loads at DOFs
F = zeros(ndof,1);
% assign point load, mod is used to always assign to even value to ensure y direction load
F(2) = -10;

%Boundary conditions
fixeddofs1 = (1:2:(nely+1)*2-1)';
fixeddofs2 = ((nely+1)*(nelx+1)*2)';
fixeddofs = unique([fixeddofs1;fixeddofs2]);
alldofs = [1:ndof]';
freedofs = setdiff(alldofs,fixeddofs);

% fixed nodes
fixedNodes1 = 1:ntolr;
fixedNodes2 = ntolr*ntoud;

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










