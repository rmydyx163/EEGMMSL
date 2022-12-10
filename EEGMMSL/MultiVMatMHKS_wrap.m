warning off;

clear

foldsTimes = 5;
ratio = 0.80;
eta = 10^(-6);
C = [0.01;0.1;1;10;100];

alpha = [0.01;0.1;1;10];
beta =  [0.01;0.1;1;10];
gamma = [0.01;0.1;1;10]; 


u_default = 0.5;
v_default = 0.5;
b_default = 10^(-8);

for i_datanum =8 
    [ data, dataname ] = import_data( i_datanum );
    

    
    [comb_mode,view_num] = mat_combined_mode(data{1}(1,:));

    for i_C = 1:length(C) 
        for i_alpha = 1:length(alpha)
            for i_beta = 1 :length(beta)       
         for i_gamma=1:length(gamma) 
            
            InputPar.C = C(i_C);
            InputPar.ratio = ratio;
            InputPar.u_d = u_default;
            InputPar.v_d = v_default;
            InputPar.b_d = b_default;
            
            InputPar.alpha=alpha(i_alpha);  
            InputPar.beta=beta(i_beta);
            InputPar.gamma=gamma(i_gamma);
            
            InputPar.eta = eta;
            
            disp(['Setting: foldsTimes-',num2str(foldsTimes),', Training_Ratio(ratio)-',num2str(ratio),' ,C-',num2str(InputPar.C),' ,view num-',num2str(view_num),' ,alpha-',num2str(InputPar.alpha),' ,beta-',num2str(InputPar.beta),' ,Dataset-',dataname]);
            disp(['--------------------------------------']); 
            

            file_name_result=['..\result\',dataname,'_Mv_Mat_all','.txt'];
            fp_result=fopen(file_name_result,'at+');
            
            
            time = zeros(foldsTimes,1);

            
            accuracy = zeros(foldsTimes,1);
            
            AUC=zeros(foldsTimes,1);
            GMM=ones(foldsTimes,1);
            
            cur_acc = zeros(foldsTimes,size(data,1) );
            AAcc=zeros(foldsTimes,1);

            
%---------------------------------Iteration loop------------------------------------------
            for k_fold = 1:foldsTimes   
                %Pairwise class training
                tic;

                testSet = data(:, k_fold); 
                trainSet = data; 
                trainSet(:, k_fold) = [];

                [train test totalClass] = getTrainAndTest(trainSet, testSet);

                for i_classone = 1:(totalClass-1)
                    for i_classtwo = (i_classone+1):totalClass
                        
                        classone = train{i_classone};
                        classtwo = train{i_classtwo};
                        
                        train_classone_num = size(classone,1);
                        train_classtwo_num = size(classtwo,1);
                        train_binary_data = [classone; classtwo];

                        matri_sample = matricize_andbroad_data( train_binary_data, comb_mode, view_num );
                        
                        [ L ] = getW1andW2andL( view_num,train_classone_num,  train_classtwo_num,  InputPar);
                        
                        [ MatStruct(i_classone,i_classtwo), r_Cmodeuvb,k_iter_v,k_iter_u ] = MultiVMatMHKS_fun(train_classone_num, train_classtwo_num, matri_sample, comb_mode, view_num, InputPar,L);
                        

%---------------------------------Second process---------------------------------
                        [ L ] = v_getW1andW2andL( view_num,train_classone_num,  train_classtwo_num,  InputPar, matri_sample, MatStruct(i_classone,i_classtwo));
                        
                        [ MatStruct(i_classone,i_classtwo), r_Cmodeuvb,k_iter_v,k_iter_u ] = MultiVMatMHKS_fun(train_classone_num, train_classtwo_num, matri_sample, comb_mode, view_num, InputPar,L);
                       

                    
                    end
                end  
                
                time(foldsTimes) = toc;
    
        %test
      
        test_data = test(:,1:end-1);
        test_label = test(:,end);
        classifier_num = totalClass * (totalClass-1)/2;
        matrix_vote = zeros(1,classifier_num+1);
        
        for p_test=1:size(test_data,1)
            p_matrix_vote = 0;
            
        for i_testone = 1:(totalClass -1)                
            for i_testtwo = (i_testone +1):totalClass    
                sum_all_view = 0;
                sigle_view = zeros(view_num,1);
                for p_view = 1:view_num 
                    A = reshape(test_data(p_test,:),comb_mode(p_view,1),comb_mode(p_view,2));
                     B=[A zeros(size(A,1),1);zeros(1,size(A,2)) 1];
                     sigle_view(p_view) = MatStruct(i_testone,i_testtwo).u{p_view}'*B*MatStruct(i_testone,i_testtwo).v{p_view};
                end 
                sum_all_view = sum(sigle_view);
                
                
                if sum_all_view >= 0
                    p_matrix_vote = cat(2,p_matrix_vote,i_testone);
                else
                    p_matrix_vote = cat(2,p_matrix_vote,i_testtwo);
                end
                
            end
        end
        matrix_vote = cat(1,matrix_vote,p_matrix_vote);

        end
        matrix_vote(1,:) = [];
        
        for i_poll = 1:size(test_data,1)
            vector_vote = matrix_vote(i_poll,2:end);
            matrix_vote(i_poll,1) = mode(vector_vote);
        end%for_i_poll
    
        accuracy(k_fold,1) = 100*(1- (length(find((test_label - matrix_vote(:,1))~=0)) / length(test_label)) );
        [~,tempt_location] = unique(test_label);
           tempt_location1=[0;tempt_location];
          for i=1:totalClass
          index_AUC = find((test_label(1+tempt_location1(i):tempt_location1(i+1)) - matrix_vote(1+tempt_location1(i):tempt_location1(i+1),1))~=0);
          AUC(k_fold)=AUC(k_fold)+100*(1-(length(index_AUC)/length(test_label(1+tempt_location1(i):tempt_location1(i+1)))));
          end
          AUC(k_fold)=AUC(k_fold)/totalClass;
          GMM(k_fold) = 0;
          
          results = test_label - matrix_vote(:,1);
          for i =1:totalClass
              ind = find(test_label == i);
              len_ind = length(ind);
              cur_acc(k_fold,i) = length(find(results(ind) == 0));
              AAcc(k_fold) = AAcc(k_fold) + cur_acc(k_fold,i)/len_ind;
          end
          AAcc(k_fold) = 100 * AAcc(k_fold)/totalClass ;
        
       
            end
    disp(['The average accuracy is: ',num2str(mean(accuracy))]);

    disp(['The average AAcc is: ',num2str(mean(AAcc))]);


    fprintf(fp_result,' %6.3f  %6.3f  %6.3f  %6.3f  %d   %3.4f   %3.4f   %3.4f  %3.4f  %3.4f  %3.4f  %3.4f  %3.4f  %3.4f\r\n',       InputPar.C,InputPar.alpha,InputPar.beta,InputPar.gamma,view_num,mean(accuracy),mean(AAcc),mean(time),mean(AUC),mean(GMM),std(accuracy),std(AUC),std(GMM),std(AAcc));

    
    clear accuracy;clear time;


    fclose(fp_result);
    delete file_id_result;
       
        
         end 
            end
        
   
    
        end
     
    end
    

end
