function [ MatStruct, r_Cmodeuvb,k_iter_v,k_iter_u ] = MultiVMatMHKS_fun( n1, n2,matri_sample, comb_mode, view_num, InputPar, L )


total_iter = 500;  
M_row=zeros(view_num,1);
M_col=zeros(view_num,1);
Y_v = cell(view_num,1);
Y_u = cell(view_num,1);
t_v = zeros(1,total_iter);  
t_u = zeros(1,total_iter);  

t_Lv =zeros(1,total_iter);
t_Lu =zeros(1,total_iter);
L_v = zeros(1,total_iter);
L_u = zeros(1,total_iter);
for p_view=1:view_num  
        M_row(p_view) = comb_mode(p_view,1);
        M_col(p_view) = comb_mode(p_view,2);
        
        u{p_view} = zeros(M_row(p_view) + 1,total_iter);     
        v{p_view} = zeros(M_col(p_view) + 1,total_iter);
        b_v{p_view} = zeros(n1+n2,total_iter ); 
        b_u{p_view} = zeros(n1+n2,total_iter ); 
        
        e_v{p_view} = zeros(n1+n2,total_iter);
        e_u{p_view} = zeros(n1+n2,total_iter);

        
        S_1{p_view} = eye(M_row(p_view));
        S_2{p_view} = eye(M_col(p_view));
        S1{p_view} = [S_1{p_view} zeros(size(S_1{p_view},1),1);zeros(1,size(S_1{p_view},2)) 0];
        S2{p_view} = [S_2{p_view} zeros(size(S_2{p_view},1),1);zeros(1,size(S_2{p_view},2)) 0];
end



I = ones(n1+n2,1);
phi = [eye(n1) zeros(n1,n2) ; zeros(n2,n1) (-1)*eye(n2)]; 
D = diag( [ones(1,n1),(n1/n2)*ones(1,n2)] );



for i_mat = 1 : view_num 
r_Cmodeuvb{i_mat}=zeros(3 +M_col(i_mat)+1  +M_row(i_mat)+1 +n1+n2,  +total_iter*2);
end


%train for v
for p_view=1:view_num  
    u{p_view}(:,1) = [InputPar.u_d * ones(M_row(p_view),1);1];
    v{p_view}(:,1) = [InputPar.v_d * ones(M_col(p_view)+1,1)];
    b_v{p_view}(:,1) = InputPar.b_d * ones(n1+n2,1);

end


for c_v=1:view_num
    for p_v = 1:(n1+n2) 

        A_view=matri_sample{c_v}{p_v}; 
        y_view=u{c_v}(:,1)'*A_view; 
        Y_v{c_v}=[Y_v{c_v};y_view];
    end
end

for c_v=1:view_num
    v_pinv_view{c_v} = pinv( Y_v{c_v}'*phi*D*phi*Y_v{c_v} + InputPar.C*S2{c_v} + Y_v{c_v}'*L{c_v, c_v}*Y_v{c_v} );
end


k_iter_v = 2;
while (k_iter_v < total_iter)
    
    for p_view=1:view_num  

    v{p_view}(:,k_iter_v) = get_v(Y_v, p_view,v_pinv_view, phi, I, view_num, v,b_v, k_iter_v, L, D );
    
    e_v{p_view}(:,k_iter_v) = Y_v{p_view}*v{p_view}(:,k_iter_v) - I -b_v{p_view}(:,k_iter_v-1);
    
    b_v{p_view}(:,k_iter_v) = (b_v{p_view}(:,k_iter_v-1) + InputPar.ratio*(e_v{p_view}(:,k_iter_v) + abs(e_v{p_view}(:,k_iter_v))));
  
    r_Cmodeuvb{p_view}(:,k_iter_v)=[InputPar.C;M_row(p_view);M_col(p_view);u{p_view}(:,1);v{p_view}(:,k_iter_v);b_v{p_view}(:,k_iter_v)];


    
    end
    
%Judge termination conditions
    if k_iter_v == total_iter
        break;
    end          
        t_v(k_iter_v) = termination_jude(b_v ,k_iter_v, view_num) ;  

        [t_Lv(k_iter_v),L_v(k_iter_v)] = v_termination_jude( view_num, k_iter_v, phi, Y_v, D ,S1, S2, L, InputPar, u, v, b_v );

    if t_v(k_iter_v) < InputPar.eta
        break;
    end
    
    k_iter_v = k_iter_v + 1;
    
end

%train for u
for final_v=1:view_num
    if k_iter_v == total_iter
    MatStruct.v{final_v} = v{final_v}(:,k_iter_v-1);
    b_u{final_v}(:,1) = b_v{final_v}(:,k_iter_v-1);
    else
    MatStruct.v{final_v} = v{final_v}(:,k_iter_v); 
    b_u{final_v}(:,1) = b_v{final_v}(:,k_iter_v);
    
    end
end



for c_v=1:view_num
    for p_v = 1:(n1+n2) 

        A_view = matri_sample{c_v}{p_v};
        y_view = (A_view *  MatStruct.v{c_v} )'; 
        Y_u{c_v}=[Y_u{c_v};y_view];
    end
end

for c_v=1:view_num
    u_pinv_view{c_v} = pinv( Y_u{c_v}'*phi*D*phi*Y_u{c_v} + InputPar.C*S1{c_v} + Y_u{c_v}'*L{c_v, c_v}*Y_u{c_v} );
end


k_iter_u = 2;
while (k_iter_u < total_iter)
    
    for p_view=1:view_num  

    u{p_view}(:,k_iter_u) = get_u(Y_u, p_view,u_pinv_view, phi, I, view_num, u,b_u, k_iter_u, L, D );
    
    e_u{p_view}(:,k_iter_u) = Y_u{p_view}*u{p_view}(:,k_iter_u) - I -b_u{p_view}(:,k_iter_u-1);
    b_u{p_view}(:,k_iter_u) = (b_u{p_view}(:,k_iter_u-1) + InputPar.ratio*(e_u{p_view}(:,k_iter_u) + abs(e_u{p_view}(:,k_iter_u))));
    r_Cmodeuvb{p_view}(:,total_iter+k_iter_u)=[InputPar.C;M_row(p_view);M_col(p_view);u{p_view}(:,k_iter_u);MatStruct.v{p_view};b_u{p_view}(:,k_iter_u)];
  



    
    end
    
%Judge termination conditions
    if k_iter_u == total_iter
        break;
    end
        
        t_u(k_iter_u) = termination_jude(b_u ,k_iter_u, view_num) ; 
        [t_Lu(k_iter_u),L_u(k_iter_u)] = u_termination_jude( view_num, k_iter_u, phi, Y_u, D ,S1, S2, L, InputPar, u, MatStruct.v, b_u );

    

    if t_u(k_iter_u) < InputPar.eta
        break;
    end
    
    k_iter_u = k_iter_u + 1;
    
end

for final_v=1:view_num
    if k_iter_u == total_iter
    MatStruct.u{final_v} = u{final_v}(:,k_iter_u-1);
   
    else
    MatStruct.u{final_v} = u{final_v}(:,k_iter_u);
  
    end
end


end

