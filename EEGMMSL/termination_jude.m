function [ b ] = termination_jude( b_v, k_iter_v, view_num )


b1 = 0;
b2 = 0;
b = 0;
for p_v = 1:view_num
    b1 = b1 + norm(b_v{p_v}(:, k_iter_v));
    b2 = b2 + norm(b_v{p_v}(:, k_iter_v-1));
end
b = abs(b2 -b1)/b2;

end

