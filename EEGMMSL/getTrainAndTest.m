function [train test, totalClass] = getTrainAndTest(trainSet, testSet)
    [totalClass, section] = size(trainSet);
    test=[];
    for i = 1:totalClass
        test = [test; testSet{i, 1}];
    end
    for i = 1:totalClass
        col_train = [];
        for j = 1:section
            col_train = [col_train; trainSet{i, j}];
        end
        train{1, i} = col_train(:, 1:end-1);
    end
end