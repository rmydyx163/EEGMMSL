function [ L ] = getW1andW2andL( view_num, n1, n2, InputPar ) 



m_W = cell(view_num);
W_r = cell(view_num,1);
D = cell(view_num,1);



for p = 1:view_num  
    for q = 1:view_num
        
        if p == q


        m_W{p,q} = 0;
        elseif p ~= q

            m_W{p,q} = 0;
           
         end   
       
    end
end 



L = m_W;


end

