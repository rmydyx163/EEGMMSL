function [ t, L2 ] = v_termination_jude( view_num, k_iter_v, phi, Y, D ,S1, S2, L, InputPar, u, v, b )

L1 = 0;
for c_v=1:view_num
    temp = phi*Y{c_v}*v{c_v}(:,k_iter_v-1) - 1 -b{c_v}(:,k_iter_v-1);
    L1 = L1 + temp'*D*temp + InputPar.C*(u{c_v}(:,1)'*S1{c_v}*u{c_v}(:,1) + v{c_v}(:,k_iter_v-1)'*S2{c_v}*v{c_v}(:,k_iter_v-1) ) ;     
end
for i = 1:view_num
    for j = 1:view_num
        temp = ( Y{i}*v{i}(:,k_iter_v-1) )' * L{i, j} * ( Y{j}*v{j}(:,k_iter_v-1) );
    end
end
L1 = L1 + temp;

L2 = 0;
for c_v=1:view_num
    temp = phi*Y{c_v}*v{c_v}(:,k_iter_v) - 1 -b{c_v}(:,k_iter_v);
    L2 = L2 + temp'*D*temp + InputPar.C*(u{c_v}(:,1)'*S1{c_v}*u{c_v}(:,1) + v{c_v}(:,k_iter_v)'*S2{c_v}*v{c_v}(:,k_iter_v) ) ;     
end
for i = 1:view_num
    for j = 1:view_num
        temp = ( Y{i}*v{i}(:,k_iter_v) )' * L{i, j} * ( Y{j}*v{j}(:,k_iter_v) );
    end
end
L2 = L2 + temp;

t = abs( L2 - L1 )/abs(L1);

end

