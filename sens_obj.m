function df = sens_obj(num_ele,dKdRho,U,N)

df = zeros(num_ele,1);
for e = 1:num_ele
%   from , dCdRho_e = -Ue' * (Te' * dKedRhoe * Te) * Ue, 
%   (Te' * dKedRhoe * Te) is inhereted as dKdRho
% FEA from https://mecheng.iisc.ac.in/suresh/me237/feaNotes/Chapter6.pdf
    df(e,1) = -1*U'*N(:,:,e)'*dKdRho(:,:,e)* N(:,:,e)*U; %???

end

