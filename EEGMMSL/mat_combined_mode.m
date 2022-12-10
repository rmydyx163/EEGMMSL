function [ valid_comb_mode,valid_comb_num ] = mat_combined_mode( vec )


while(1)
    
    [r_vec,c_vec]=size(vec);
    c_vec = c_vec-1; 

    ele_num = 0;   
    if(c_vec == 1)
        ele_num = r_vec;
    elseif(r_vec == 1)
        ele_num = c_vec;
    else
        disp('The input value is not a vector! It is a matrix! ');
        break;
    end
    
    
    valid_comb_mode = []; 
    valid_comb_num = 0;
    for n = 1:ele_num 
         tempt_value = mod(ele_num,n);
        if(tempt_value == 0)
            valid_comb_num = valid_comb_num+1;
            r_combine_num = ele_num/n;
            c_combine_num = n;
            valid_comb_mode = [valid_comb_mode; [r_combine_num, c_combine_num] ];
        else
            continue;
        end 
    end 


    
    if(valid_comb_num == 0) 
        disp('This vector can not be matrixlized!');
        break
    end
    
    break;
   
end
