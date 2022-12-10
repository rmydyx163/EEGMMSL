function [ matricized_sample ] = matricize_andbroad_data( train_binary_data, comb_mod, view_num )


for c_view_num=1:view_num
    for c_sample=1:size(train_binary_data,1)
        M_row_view = comb_mod(c_view_num,1);
        M_col_view = comb_mod(c_view_num,2);
        A_view=reshape(train_binary_data(c_sample,:),M_row_view,M_col_view);
        matricized_sample{c_view_num}{c_sample}=[A_view zeros(size(A_view,1),1);zeros(1,size(A_view,2)) 1];
    end
end

end

