
function [L,C,S,nod_coor] = moveNodes(xp,num_ele,nod_coor,nod_coor0,szmv,mNodes,ele_nod)


for i = 1:length(mNodes)

    % construct node locations, considering node numbers are built from
    % top to bottom, but using up=positive convention in y direction
    % add displacements from placement variables xp

    nod_coor(mNodes(i),1) = nod_coor0(mNodes(i),1)+szmv*(xp( (i)*2-1 )-0.5); %szmv*(xp( (i)*2-1 )-0.5);
    nod_coor(mNodes(i),2) = nod_coor0(mNodes(i),2)+szmv*(xp( (i)*2   )-0.5); %szmv*(xp( (i)*2   )-0.5)

end

% FEA Parameters
L = zeros(num_ele,1);
C = zeros(num_ele,1);
S = zeros(num_ele,1);

for e=1:num_ele
 L(e)=sqrt((nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))^2+...
      (nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))^2);
 C(e)=(nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))/L(e);
 S(e)=(nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))/L(e);
end

end









