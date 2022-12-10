function [ data, dataname ] = import_data( i_datanum )



switch i_datanum 

    case 8
        load ..\datasets_UCI\iris_v.mat; 
        data = iris_v;
        dataname='iris_v';

    case 41
        load ..\datasets_UCI\ionosphere_v.mat;
        data = ionosphere_v;
        dataname='ionosphere_v';

    otherwise
        disp('You have input wrong dataset number!');
end



end

