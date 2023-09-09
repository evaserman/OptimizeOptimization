function [df, dv] = sens_xy(U, ele_nod, nod_coor, L, Ke, x_t, A0, E, ...
    mNodes, szmv, N, dNdx, dNdy)

n_mNodes = length(mNodes);
temp_df = cell(n_mNodes, 1);
temp_dv = cell(n_mNodes, 1);
Mm = [1 0 -1 0;
      0 0 0 0;
      -1 0 1 0;
      0 0 0 0];

parfor ii = 1:n_mNodes
    df_local = zeros(n_mNodes*2, 1);
    dv_local = zeros(n_mNodes*2, 1);
    
    % All elements attached to node n
    i = mNodes(ii);
    Ev1 = find(ele_nod(:,1)==i);
    Ev2 = find(ele_nod(:,2)==i);
    Ev = unique([Ev1; Ev2]);

    % preallocate sum for obj sens
    EvSum_x = zeros(length(Ev),1);
    EvSum_y = zeros(length(Ev),1);

    % preallocate sum for constr sens
    EvVolSum_x = zeros(length(Ev),1);
    EvVolSum_y = zeros(length(Ev),1);

        for e = 1:length(Ev)

            % coordinates of first node of element j
            coor_1 = nod_coor(ele_nod(Ev(e),1),:);
            % coordinates of second node of element j
            coor_2 = nod_coor(ele_nod(Ev(e),2),:);

            % simplify variable names
            x1 = coor_1(1);  y1 = coor_1(2);  
            x2 = coor_2(1);  y2 = coor_2(2);
            % Make Te
            C = (x2-x1)/L(Ev(e));
            S = (y2-y1)/L(Ev(e));
            Te = [ C S  0 0;
                  -S C  0 0
                   0 0  C S
                   0 0 -S C];

            if ele_nod(Ev(e),1) == i %(i-1)*2+1
            % Calculate dTdx and dTdy with respect to x1 and y1
            % these are from wolfram alpha
% https://www.wolframalpha.com/input?i=derivative+of+++%28x2-x1%29%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+x1
            dCdx = szmv*(-1)*( (y2-y1)^2 )/ ...
              ( ( (x2-x1)^2 + (y2-y1)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28x2-x1%29%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+y1
            dCdy = szmv*(-1)*( (x2-x1)*(y1-y2) )/ ...
              ( ( (x1-x2)^2 + (y1-y2)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28y2-y1%29%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+x1
            dSdx = szmv*(-1)*( (x2-x1)*(y1-y2) )/ ...
              ( ( (x1-x2)^2 + (y1-y2)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28y2-y1%29%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+y1
            dSdy = szmv*(-1)*( (x2-x1)^2 )/ ...
              ( ( (x2-x1)^2 + (y2-y1)^2 )^(3/2) );

            % calculate dKdx and dKdy
            % Ke is local stiffness, so dKe/dx is just affecting L in denom. of EAx/L
            dKdx = szmv*(A0*E)* x_t(Ev(e)) * 1*(x2-x1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(3/2)  );
            dKdy = szmv*(A0*E)* x_t(Ev(e)) * 1*(y2-y1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(3/2)  ); %%

            % Calculate volume constraint gradients, capturing change in L
            EvVolSum_x(e) = szmv*(A0)* x_t(Ev(e))* (-1) *(x2-x1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(1/2)  );
            EvVolSum_y(e) = szmv*(A0)* x_t(Ev(e))* (-1) *(y2-y1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(1/2)  );

            else
            % Calculate dTdx and dTdy with respect to x2 and y2
            % these are from wolfram alpha
% https://www.wolframalpha.com/input?i=derivative+of+++%28x2-x1%29+%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+x2
            dCdx = szmv*( 1)*( (y2-y1)^2 )/ ...
              ( ( (x2-x1)^2 + (y2-y1)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28x2-x1%29+%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+y2
            dCdy = szmv*( 1)*( (x2-x1)*(y1-y2) )/ ...
              ( ( (x1-x2)^2 + (y1-y2)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28y2-y1%29+%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+x2
            dSdx = szmv*( 1)*( (x2-x1)*(y1-y2) )/ ...
              ( ( (x1-x2)^2 + (y1-y2)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28y2-y1%29+%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+y2
            dSdy = szmv*( 1)*( (x2-x1)^2 )/ ...
              ( ( (x2-x1)^2 + (y2-y1)^2 )^(3/2) );

            % calculate dKdx and dKdy
            % Ke is local stiffness, so dKe/dx is just affecting L in denom. of EA/L
            dKdx = szmv*(A0*E)* x_t(Ev(e)) * -1*(x2-x1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(3/2)  );
            dKdy = szmv*(A0*E)* x_t(Ev(e)) * -1*(y2-y1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(3/2)  ); %%

            % Calculate volume constraint gradients, capturing change in L
            EvVolSum_x(e) = szmv*(A0)* x_t(Ev(e)) * (x2-x1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(1/2)  );
            EvVolSum_y(e) = szmv*(A0)* x_t(Ev(e)) * (y2-y1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(1/2)  );

            end

            dTdx = [ dCdx dSdx    0    0; 
                    -dSdx dCdx    0    0;
                     0     0    dCdx dSdx;
                     0     0   -dSdx dCdx];

            dTdy = [ dCdy dSdy    0     0;
                    -dSdy dCdy    0     0;
                      0     0    dCdy  dSdy
                      0     0   -dSdy  dCdy];

            % Calculate inside of sensitivity eq, before convolution
            local_x = zeros(4,4,4);
            local_x(:,:,1) = Te'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*Te;
            local_x(:,:,2) = dTdx'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*Te;
            local_x(:,:,3) = Te' *(dKdx*Mm)*Te; %already has x_t baked in
            local_x(:,:,4) = Te'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*dTdx;

            % repeat for y
            local_y = zeros(4,4,4);
            local_y(:,:,1) = Te'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*Te;
            local_y(:,:,2) = dTdy'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*Te;
            local_y(:,:,3) = Te' * (dKdy*Mm)*Te;
            local_y(:,:,4) = Te'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*dTdy;

            % Incorprorate covolution into sensitivity calc
            N1 = N(:,:,Ev(e));
            dNdx1 = dNdx{Ev(e),i};
            part1_x = dNdx1'*local_x(:,:,1)*N1;
            part2_x = N1'   *local_x(:,:,2)*N1;
            part3_x = N1'   *local_x(:,:,3)*N1;
            part4_x = N1'   *local_x(:,:,4)*N1;
            part5_x = N1'   *local_x(:,:,1)*dNdx1;
           
            % repeat for y
            dNdy1 = dNdy{Ev(e),i};
            part1_y = dNdy1'*local_y(:,:,1)*N1;
            part2_y = N1'   *local_y(:,:,2)*N1;
            part3_y = N1'   *local_y(:,:,3)*N1;
            part4_y = N1'   *local_y(:,:,4)*N1;
            part5_y = N1'   *local_y(:,:,1)*dNdy1;

            % sum 5 components above
            dKedxn = part1_x+part2_x+part3_x+part4_x+part5_x;
            dKedyn = part1_y+part2_y+part3_y+part4_y+part5_y;

            % mult by Ue
            EvSum_x(e) = -1*U'*dKedxn * U;
            EvSum_y(e) = -1*U'*dKedyn * U;

        end
     % Aggregate your results from the inner loop
    df_local((ii-1)*2+1) = sum(EvSum_x);
    df_local(ii*2) = sum(EvSum_y);

    dv_local((ii-1)*2+1) = sum(EvVolSum_x);
    dv_local(ii*2) = sum(EvVolSum_y)
    
    temp_df{ii} = df_local;
    temp_dv{ii} = dv_local;
end

% Aggregating the results from the parfor loop
df = zeros(n_mNodes*2, 1);
dv = zeros(n_mNodes*2, 1);
for ii = 1:n_mNodes
    df = df + temp_df{ii};
    dv = dv + temp_dv{ii};
end
end
%}
%{
df = zeros(length(mNodes)*2,1);
dv = zeros(length(mNodes)*2,1);
Mm = [ 1 0 -1 0; 
       0 0  0 0; 
      -1 0  1 0; 
       0 0  0 0];

% derivative with convolution is:
% dK+/dx = dN'/dx*K*N + N'*dK/dx*N + N'*K*dN/dx
% First, calculate dKdX, then use results to calculate dK+/dx
% computation following Xia et al 2013 
% Then dK+/dx computation following Zegard and Paulino 2013

% dK/dx
for ii = 1:length(mNodes)                %PARALLELIZE?!?!?!?
    % All elements attached to node n
    i = mNodes(ii);
    Ev1 = find(ele_nod(:,1)==i);
    Ev2 = find(ele_nod(:,2)==i);
    Ev = unique([Ev1; Ev2]);

    % preallocate sum for obj sens
    EvSum_x = zeros(length(Ev),1);
    EvSum_y = zeros(length(Ev),1);

    % preallocate sum for constr sens
    EvVolSum_x = zeros(length(Ev),1);
    EvVolSum_y = zeros(length(Ev),1);

        for e = 1:length(Ev)

            % coordinates of first node of element j
            coor_1 = nod_coor(ele_nod(Ev(e),1),:);
            % coordinates of second node of element j
            coor_2 = nod_coor(ele_nod(Ev(e),2),:);

            % simplify variable names
            x1 = coor_1(1);  y1 = coor_1(2);  
            x2 = coor_2(1);  y2 = coor_2(2);
            % Make Te
            C = (x2-x1)/L(Ev(e));
            S = (y2-y1)/L(Ev(e));
            Te = [ C S  0 0;
                  -S C  0 0
                   0 0  C S
                   0 0 -S C];

            if ele_nod(Ev(e),1) == i %(i-1)*2+1
            % Calculate dTdx and dTdy with respect to x1 and y1
            % these are from wolfram alpha
% https://www.wolframalpha.com/input?i=derivative+of+++%28x2-x1%29%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+x1
            dCdx = szmv*(-1)*( (y2-y1)^2 )/ ...
              ( ( (x2-x1)^2 + (y2-y1)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28x2-x1%29%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+y1
            dCdy = szmv*(-1)*( (x2-x1)*(y1-y2) )/ ...
              ( ( (x1-x2)^2 + (y1-y2)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28y2-y1%29%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+x1
            dSdx = szmv*(-1)*( (x2-x1)*(y1-y2) )/ ...
              ( ( (x1-x2)^2 + (y1-y2)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28y2-y1%29%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+y1
            dSdy = szmv*(-1)*( (x2-x1)^2 )/ ...
              ( ( (x2-x1)^2 + (y2-y1)^2 )^(3/2) );

            % calculate dKdx and dKdy
            % Ke is local stiffness, so dKe/dx is just affecting L in denom. of EAx/L
            dKdx = szmv*(A0*E)* x_t(Ev(e)) * 1*(x2-x1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(3/2)  );
            dKdy = szmv*(A0*E)* x_t(Ev(e)) * 1*(y2-y1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(3/2)  ); %%

            % Calculate volume constraint gradients, capturing change in L
            EvVolSum_x(e) = szmv*(A0)* x_t(Ev(e))* (-1) *(x2-x1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(1/2)  );
            EvVolSum_y(e) = szmv*(A0)* x_t(Ev(e))* (-1) *(y2-y1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(1/2)  );

            else
            % Calculate dTdx and dTdy with respect to x2 and y2
            % these are from wolfram alpha
% https://www.wolframalpha.com/input?i=derivative+of+++%28x2-x1%29+%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+x2
            dCdx = szmv*( 1)*( (y2-y1)^2 )/ ...
              ( ( (x2-x1)^2 + (y2-y1)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28x2-x1%29+%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+y2
            dCdy = szmv*( 1)*( (x2-x1)*(y1-y2) )/ ...
              ( ( (x1-x2)^2 + (y1-y2)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28y2-y1%29+%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+x2
            dSdx = szmv*( 1)*( (x2-x1)*(y1-y2) )/ ...
              ( ( (x1-x2)^2 + (y1-y2)^2 )^(3/2) );
% https://www.wolframalpha.com/input?i=derivative+of+++%28y2-y1%29+%2F+%28+%28+%28x2-x1%29%5E2+%2B+%28y2-y1%29%5E2+%29%5E%281%2F2%29+%29+with+respect+to+y2
            dSdy = szmv*( 1)*( (x2-x1)^2 )/ ...
              ( ( (x2-x1)^2 + (y2-y1)^2 )^(3/2) );

            % calculate dKdx and dKdy
            % Ke is local stiffness, so dKe/dx is just affecting L in denom. of EA/L
            dKdx = szmv*(A0*E)* x_t(Ev(e)) * -1*(x2-x1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(3/2)  );
            dKdy = szmv*(A0*E)* x_t(Ev(e)) * -1*(y2-y1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(3/2)  ); %%

            % Calculate volume constraint gradients, capturing change in L
            EvVolSum_x(e) = szmv*(A0)* x_t(Ev(e)) * (x2-x1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(1/2)  );
            EvVolSum_y(e) = szmv*(A0)* x_t(Ev(e)) * (y2-y1) / (  ( (x2-x1)^2 + (y2-y1)^2 )^(1/2)  );

            end

            dTdx = [ dCdx dSdx    0    0; 
                    -dSdx dCdx    0    0;
                     0     0    dCdx dSdx;
                     0     0   -dSdx dCdx];

            dTdy = [ dCdy dSdy    0     0;
                    -dSdy dCdy    0     0;
                      0     0    dCdy  dSdy
                      0     0   -dSdy  dCdy];

            % Calculate inside of sensitivity eq, before convolution
            local_x = zeros(4,4,4);
            local_x(:,:,1) = Te'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*Te;
            local_x(:,:,2) = dTdx'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*Te;
            local_x(:,:,3) = Te' *(dKdx*Mm)*Te; %already has x_t baked in
            local_x(:,:,4) = Te'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*dTdx;

            % repeat for y
            local_y = zeros(4,4,4);
            local_y(:,:,1) = Te'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*Te;
            local_y(:,:,2) = dTdy'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*Te;
            local_y(:,:,3) = Te' * (dKdy*Mm)*Te;
            local_y(:,:,4) = Te'*(Ke(Ev(e))*x_t(Ev(e))*Mm)*dTdy;

            % Incorprorate covolution into sensitivity calc
            N1 = N(:,:,Ev(e));
            dNdx1 = dNdx{Ev(e),i};
            part1_x = dNdx1'*local_x(:,:,1)*N1;
            part2_x = N1'   *local_x(:,:,2)*N1;
            part3_x = N1'   *local_x(:,:,3)*N1;
            part4_x = N1'   *local_x(:,:,4)*N1;
            part5_x = N1'   *local_x(:,:,1)*dNdx1;
           
            % repeat for y
            dNdy1 = dNdy{Ev(e),i};
            part1_y = dNdy1'*local_y(:,:,1)*N1;
            part2_y = N1'   *local_y(:,:,2)*N1;
            part3_y = N1'   *local_y(:,:,3)*N1;
            part4_y = N1'   *local_y(:,:,4)*N1;
            part5_y = N1'   *local_y(:,:,1)*dNdy1;

            % sum 5 components above
            dKedxn = part1_x+part2_x+part3_x+part4_x+part5_x;
            dKedyn = part1_y+part2_y+part3_y+part4_y+part5_y;

            % mult by Ue
            EvSum_x(e) = -1*U'*dKedxn * U;
            EvSum_y(e) = -1*U'*dKedyn * U;

        end

       % then sum all of these 3-part-sums  with U mult to get dLdxi and dLdyi
       df((ii-1)*2+1) = sum(EvSum_x); 
       df(ii*2)       = sum(EvSum_y); 

       % then sum all vol gradient calcs to get dVdxi and dVdyi
       dv((ii-1)*2+1) = sum(EvVolSum_x); 
       dv(ii*2)       = sum(EvVolSum_y); 

end

% % scale dv
% dv = dv/(vmax_t*sum(L0*A0));
dv = dv;

end
%}
