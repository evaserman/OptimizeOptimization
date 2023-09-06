function [K_t,dKdx,Ke,N,dNdy,dNdx] = trussK(nele_t,C,S,A0,x_t,E_t,L,stiffness,...
    num_nod,nelx,nely,nod_coor,ele_nod,Rc)

dKdx    = zeros(4,4,nele_t);
Ke      = zeros(nele_t,1);
N0      = zeros(4,(nelx+1)*(nely+1)*2,nele_t);
N       = zeros(4,(nelx+1)*(nely+1)*2,nele_t);
dNdx_0  = zeros(4,(nelx+1)*(nely+1)*2);
dNdx    = cell(nele_t,num_nod);         %accessible bins for sens matrices
dNdy_0  = zeros(4,(nelx+1)*(nely+1)*2);
dNdy    = cell(nele_t,num_nod);         %accessible bins for sens matrices

% create coordinates of continuum elements
[x,y] = meshgrid(1:nelx+1,1:nely+1);
s = [x(:)-1  nely+1-y(:)];

for e=1:nele_t

 Ce = C(e);
 Se = S(e);

Ke(e)=(A0*x_t(e)*E_t/L(e));
k=Ke(e)*[Ce*Ce Ce*Se -Ce*Ce -Ce*Se;Ce*Se Se*Se -Ce*Se -Se*Se;...
        -Ce*Ce -Ce*Se Ce*Ce Ce*Se; -Ce*Se -Se*Se Ce*Se Se*Se];

 % save dK_dRho for sens calc, = (Te' * dKedRhoe * Te)    
dKdx(:,:,e) = (A0*E_t/L(e)*[Ce*Ce Ce*Se -Ce*Ce -Ce*Se;Ce*Se Se*Se -Ce*Se -Se*Se;...
   -Ce*Ce -Ce*Se Ce*Ce Ce*Se; -Ce*Se -Se*Se Ce*Se Se*Se]);
   
% Build N
for i=1:2     % i is node 1, then node 2 of truss element
    for j=1:(nelx+1)*(nely+1)    % j is continuum nodes
        
        % establish coords of truss and cont nodes
        x1 = nod_coor(ele_nod(e,i),1); x2 = s(j,1);
        y1 = nod_coor(ele_nod(e,i),2); y2 = s(j,2); 

        dsr = sqrt( (x2-x1)^2 + (y2-y1)^2 ) ; %distance

        if dsr <= Rc 

        N0(i*2-1,j*2-1,e) = 1-sin(dsr*pi/(2*Rc));
        N0(i*2  ,j*2  ,e) = 1-sin(dsr*pi/(2*Rc));

            if dsr > 0 % avoid 0/0=nan if dis = 0
                % sens calc
                % https://www.wolframalpha.com/input?i=derivative+of+-sin%28sqrt%28%28x2-t%29%5E2%2B%28y2-y1%29%5E2%29*pi%2F%282*R%29%29+with+respect+to+t    
                dNdx_0(i*2-1,j*2-1) = ( (pi/(2*Rc)) * cos(dsr*pi/(2*Rc)) )*(x2-x1)/dsr;
                dNdx_0(i*2  ,j*2  ) = ( (pi/(2*Rc)) * cos(dsr*pi/(2*Rc)) )*(x2-x1)/dsr;
                % https://www.wolframalpha.com/input?i=derivative+of+-sin%28sqrt%28%28x2-x1%29%5E2%2B%28y2-t%29%5E2%29*pi%2F%282*R%29%29+with+respect+to+t
                dNdy_0(i*2-1,j*2-1) = ( (pi/(2*Rc)) * cos(dsr*pi/(2*Rc)) )*(y2-y1)/dsr;
                dNdy_0(i*2  ,j*2  ) = ( (pi/(2*Rc)) * cos(dsr*pi/(2*Rc)) )*(y2-y1)/dsr;
            end

        end
 
    end

dNdx{e,ele_nod(e,i)} = dNdx_0; %stash sens in proper bin
dNdy{e,ele_nod(e,i)} = dNdy_0; 
dNdx_0  = zeros(4,(nelx+1)*(nely+1)*2); % reset matrix
dNdy_0  = zeros(4,(nelx+1)*(nely+1)*2); 

end

% normalize magnitudes of convolution values
mags = sum(N0(:,:,e),2);
% N(:,:,e) = N0(:,:,e);
N(:,:,e) = N0(:,:,e)./repmat(mags,1,(nelx+1)*(nely+1)*2);

% add element stiffness to global stiffness using convolution
stiffness = stiffness + N(:,:,e)'*k*N(:,:,e);

% modify sensitivities using normalization mags as well
for i2 = 1:2
    dNdx{e,ele_nod(e,i2)} = dNdx{e,ele_nod(e,i2)}./repmat(mags,1,(nelx+1)*(nely+1)*2);
    dNdy{e,ele_nod(e,i2)} = dNdy{e,ele_nod(e,i2)}./repmat(mags,1,(nelx+1)*(nely+1)*2);
end

end

K_t = stiffness;











