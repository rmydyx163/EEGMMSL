function [ u1 ] = get_u( Y, p_view,u_pinv_view, phi, I, view_num, u,b_u, k_iter_u, L, D )


term_sum=zeros( size( Y{1},1), 1 );

for c_v=1:view_num
    
    if c_v==p_view
        continue
    end
    
    term_sum = term_sum + L{p_view,c_v}*Y{c_v}*u{c_v}(:,k_iter_u-1) ;  
    
end

u1 = u_pinv_view{p_view}*Y{p_view}'*( phi*D*(I + b_u{p_view}(:,k_iter_u-1)) - term_sum );
u1=u1./norm(u1);  



end

