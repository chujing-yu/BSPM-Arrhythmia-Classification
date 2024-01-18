function [y_predict,models] = MultiSvm(X_train, y_train, X_test)
% multi svm
% one vs all 模型
% Input:
% X_train: n*m矩阵 n为训练集样本数 m为特征数
% y_train: n*1向量 为训练集label，支持任意多种类
% X_test: n*m矩阵 n为测试集样本数 m为特征数
% Output:
% y_predict: n*1向量 测试集的预测结果
% 
% Copyright(c) lihaoyang 2020
%

    y_labels = unique(y_train);
    n_class = size(y_labels, 1);
    models = cell(n_class, 1);
    % 训练n个模型
    for i = 1:n_class
        %{
        class_i_place = find(y_train == y_labels(i));
        svm_train_x = X_train(class_i_place,:);
        sample_num = numel(class_i_place);
        class_others = find(y_train ~= y_labels(i));
        randp = randperm(numel(class_others));
        svm_train_minus = randp(1:sample_num)';
        svm_train_x = [svm_train_x; X_train(svm_train_minus,:)];
        svm_train_y = [ones(sample_num, 1); -1*ones(sample_num, 1)];
        disp(['生成模型：', num2str(i)])
        models{i} = fitcsvm(svm_train_x, svm_train_y);
        %}
        class_i_place=find(y_train==y_labels(i));
        positive_train_x=X_train(class_i_place,:);
        positive_num=numel(class_i_place);
        class_others_place=find(y_train~=y_labels(i));
        negative_train_x=X_train(class_others_place,:);
        negative_num=numel(class_others_place);
        svm_train_x=[positive_train_x;negative_train_x];
        svm_train_y=[ones(positive_num,1);-1.*ones(negative_num,1)];
        disp(['生成模型：', num2str(i)]);
        models{i}=fitcsvm(svm_train_x,svm_train_y);
    end
    test_num = size(X_test, 1);
    y_predict = zeros(test_num, 1);
    % 对每条数据，n个模型分别进行预测，选择label为1且概率最大的一个作为预测类别
    for i = 1:test_num
        if mod(i, 100) == 0
            disp(['预测个数：', num2str(i)])
        end
        bagging = zeros(n_class, 1);
        for j = 1:n_class
            model = models{j};
            [label, rat] = predict(model, X_test(i,:));
            bagging(j) = bagging(j) + rat(2);
        end
        [maxn, maxp] = max(bagging);
        y_predict(i) = y_labels(maxp);
    end
end