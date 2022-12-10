function [ t, L2 ] = u_termination_jude( view_num, k_iter_u, phi, Y, D ,S1, S2, L, InputPar, u, v, b )


L1 = 0;
for c_v=1:view_num
    temp = phi*Y{c_v}*u{c_v}(:,k_iter_u-1) - 1 -b{c_v}(:,k_iter_u-1);
    L1 = L1 + temp'*D*temp + InputPar.C*(u{c_v}(:,k_iter_u-1)'*S1{c_v}*u{c_v}(:,k_iter_u-1) + v{c_v}'*S2{c_v}*v{c_v} ) ;     
end
for i = 1:view_num
    for j = 1:view_num
        temp = ( Y{i}*u{i}(:,k_iter_u-1) )' * L{i, j} * ( Y{j}*u{j}(:,k_iter_u-1) );
    end
end
L1 = L1 + temp;

L2 = 0;
for c_v=1:view_num
    temp = phi*Y{c_v}*u{c_v}(:,k_iter_u) - 1 -b{c_v}(:,k_iter_u);
    L2 = L2 + temp'*D*temp + InputPar.C*(u{c_v}(:,k_iter_u)'*S1{c_v}*u{c_v}(:,k_iter_u) + v{c_v}'*S2{c_v}*v{c_v} ) ;     
end
for i = 1:view_num
    for j = 1:view_num
        temp = ( Y{i}*u{i}(:,k_iter_u) )' * L{i, j} * ( Y{j}*u{j}(:,k_iter_u) );
    end
end
L2 = L2 + temp;

t = abs( L2 - L1 )/abs(L1);

end

