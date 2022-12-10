function [ L ] = v_getW1andW2andL( view_num, n1, n2, InputPar,matri_sample, MatStruct) 



m_W = cell(view_num);
W_r = cell(view_num,1);
D = cell(view_num,1);
Y = cell(view_num,1);


for c_v=1:view_num
    for p_v = 1:(n1+n2) 

        A_view=matri_sample{c_v}{p_v}; 
        y_view=MatStruct.u{c_v}'*A_view; 
        Y{c_v}=[Y{c_v};y_view];
    end
end


for p = 1:view_num  
    output_1 = Y{p}(1:n1,:)*MatStruct.v{p};
    output_2 = Y{p}(n1+1:end,:)*MatStruct.v{p} ;  
    
    for q = 1:view_num
        if p == q
            
            m_W{p,q} =  zeros(n1+n2,n1+n2);
            
            for i = 1:n1+n2 
                if i<=n1
                    
                    index_1 = find( sign(output_1 )~=sign(Y{q}(i,:)*MatStruct.v{q}) ); 
                    index_2 = find( sign(output_2)==sign(Y{q}(i,:)*MatStruct.v{q}) );
                m_W{p,q}(  index_1  , i  ) = abs(output_1(index_1)) / norm( output_1(index_1) ) * InputPar.alpha;%Y{q}(i,:)*MatStruct.v{q} 
                m_W{p,q}(  n1+index_2   , i  ) = -abs(output_2(index_2)) / norm( output_2(index_2) ) * InputPar.beta;
                else
                    index_1 = find( sign(output_1 )==sign(Y{q}(i,:)*MatStruct.v{q}) );
                    index_2 = find( sign(output_2)~=sign(Y{q}(i,:)*MatStruct.v{q}) );
                m_W{p,q}(  index_1  , i  ) = -abs(output_1(index_1)) / norm( output_1(index_1) ) * InputPar.beta;
                m_W{p,q}(  n1+index_2   , i  ) = abs(output_2(index_2)) / norm( output_2(index_2) ) * InputPar.alpha;
                 end
            end
            m_W{p,q} =  - m_W{p,q};
            m_W{p,q} = sparse( m_W{p,q});
        
        elseif p ~= q
            
            m_W{p,q} =zeros(n1+n2,n1+n2);
            
            for i = 1:n1+n2 
                if i<=n1
                    
                    index_1 = find( sign(output_1(i))~=sign(Y{q}(i,:)*MatStruct.v{q}) );
                    if index_1 ==1
                        m_W{p,q}(  i  , i  ) = 1;
                    end
                
                
                else
                    
                    index_2 = find( sign(output_2(i-n1))~=sign(Y{q}(i,:)*MatStruct.v{q}) );
                    if index_2 ==1
                        m_W{p,q}(  i  , i  ) = 1;
                    end
                    
                end
            end
            m_W{p,q} =  - InputPar.gamma * m_W{p,q};
            m_W{p,q} = sparse( m_W{p,q});
            
        end
    end
end
        


for p = 1:view_num
    W_r{p} = -m_W{p,1};
    for q = 2:view_num
        W_r{p} = W_r{p} - m_W{p,q};
    end
    D{p} = sum(W_r{p},2);
end


for p =1: view_num
    for i = 1:n1+n2
    m_W{p,p}(i,i) = D{p}(i) + m_W{p,p}(i,i);
    end     
    
end

L = m_W;


end

