function [ v1 ] = get_v( Y, p_view,v_pinv_view, phi, I, view_num, v,b_v, k_iter_v,L,D )


term_sum=zeros( size( Y{1},1), 1 );

for c_v=1:view_num
    
    if c_v==p_view
        continue
    end
    
    term_sum = term_sum + L{p_view,c_v}*Y{c_v}*v{c_v}(:,k_iter_v-1);  
    
end


v1 = v_pinv_view{p_view}*Y{p_view}'*( phi*D*(I + b_v{p_view}(:,k_iter_v-1)) - term_sum );
v1=v1./norm(v1);  


end

