clear 
clc
fs=500;          %采样频率
% 设置文件夹路径
folderPath = 'E:\新华医院202107心电信号\dataset';
%信号长度改成10S，MATLAB自带的SVM等分类器分类
% 获取文件夹中的所有文件
fileList = dir(fullfile(folderPath, '*.dat'));
%特征矩阵
dataset=[];

%每位患者取的信号段数
segment_num=10;
% 遍历文件列表
for i = 1:length(fileList)
    % 获取当前文件的完整路径
    filePath = fullfile(folderPath, fileList(i).name);
    
    % 读取ECG信号数据
    fid=fopen(filePath);
    A = fread(fid,inf,'float');
    lenperchan = floor(length(A)/133);
    fprintf(1,'%g\n',lenperchan);
    x_date = reshape(A,133,lenperchan);
    channel_number_x=24; %选取通道编号25.50
    x_signal=x_date(channel_number_x,1:end);
    x_signal=x_signal+10000;
    %每位患者取10段10s的信号
    preprocessed_signal=[];
    segment_start=1;
    segment_len=5000;
    segment_interval=1000;
        
    %每位患者的小波包SVD,MAX,STD特征
    decompose_level=4;
    %WPD_Charact_type=3;
    WPD_Charact_num=power(2,decompose_level);
    
    total_charact_num=57;
    one_paient=zeros(segment_num,total_charact_num);

    for n=1:segment_num
        segment_start_present=segment_start+(n-1)*(segment_len+segment_interval)
        segment_signal=x_signal(1,segment_start_present:(segment_start_present+segment_len-1));
        preprocessed_segment_signal=preprocess(segment_signal,fs);
        preprocessed_signal=[preprocessed_signal;preprocessed_segment_signal]; %每一行是一段10s的信号

        flag=0;
        [segment_ecg_filter,segment_Signal_filter_index,segment_Signal_filter]=Find_R_Peaks(preprocessed_segment_signal,fs,flag);
        [seg_Amp,seg_SD,seg_CF,seg_SF,seg_IF]=Seg_Time_Charact(segment_ecg_filter,segment_Signal_filter,segment_len);
        
        one_paient(n,1)=seg_Amp;
        one_paient(n,2)=seg_SD;
        one_paient(n,3)=seg_CF;
        one_paient(n,4)=seg_SF;
        one_paient(n,5)=seg_IF;

        [seg_fft_SM,seg_fft_SSD,seg_fft_SD,seg_fft_SK]=Seg_Freq_Charact(segment_ecg_filter,segment_len);
        
        one_paient(n,6)=seg_fft_SM;
        one_paient(n,7)=seg_fft_SSD;
        one_paient(n,8)=seg_fft_SD;
        one_paient(n,9)=seg_fft_SK;

        wavename='db6';
        entropy_name='shannon';
        [seg_SVD,seg_MAX,seg_STD]=Seg_WPD_Charact(segment_ecg_filter,wavename,decompose_level,entropy_name);

        one_paient(n,10:25)=seg_SVD;
        one_paient(n,26:41)=seg_MAX;
        one_paient(n,42:57)=seg_STD;

    end

    %% 时频域、小波域特征
    dataset=[dataset;one_paient];

end
%% 归一化
sample_size=length(fileList)*segment_num;      %总的样本数量
charact_type=57;                   %总的特征数量 9+16*3
%分类标签
labels=importdata('Labels.mat');
labels=repmat(labels,segment_num,1);
labels=reshape(labels,[],1);
%dataset
dataset_row_num=sample_size;   %行数为患者数
dataset_col_num=charact_type;  %列数为特征数
normalized_dataset=zeros(dataset_row_num,dataset_col_num); %归一化的数据集
for k=1:dataset_col_num
    col_mean=mean(dataset(:,k)); %列均值
    col_max=max(dataset(:,k));   %列最大值
    col_min=min(dataset(:,k));   %列最小值
    normalized_dataset(:,k)=(dataset(:,k)-col_mean)./(col_max-col_min);
end

%% SVM
%从中抽取0.2作为测试集
test_set_num=floor(dataset_row_num*0.2);
x_test=zeros(test_set_num,dataset_col_num);
[x_test,test_index]=datasample(normalized_dataset,test_set_num,1,'Replace',false);
y_test=labels(test_index,:);
x_train=normalized_dataset;
x_train(test_index,:)=[];
y_train=labels;
y_train(test_index,:)=[];
x_train=real(x_train);
x_test=real(x_test);
% 特征归一化
% label采用二进制编码
[y_predict,models]=MultiSvm(x_train,y_train,x_test);

[results,~]=confusionmat(y_test,y_predict);
% % 计算-1类的评价值
% c1_precise = results(1,1)/(results(1,1) + results(2,1));
% c1_recall = results(1,1)/(results(1,1) + results(1,2));
% c1_F1 = 2 * c1_precise * c1_recall/(c1_precise + c1_recall);
%  
% % 计算1类的评价值
% c2_precise = results(2,2)/(results(1,2) + results(2,2));
% c2_recall = results(2,2)/(results(2,1) + results(2,2));
% c2_F1 = 2 * c2_precise * c2_recall/(c2_precise + c2_recall);
%% 评价
figure;
cm=confusionchart(y_test,y_predict);

result=confusionmat(y_test,y_predict);
%% 随机森林
%拆分为训练集和测试集
trainRatio=0.8;
[trainInd,~,testInd]=dividerand(size(normalized_dataset,1),trainRatio,0,1);
XTrain=normalized_dataset(trainInd,:);
YTrain=labels(trainInd);
XTest=normalized_dataset(testInd,:);
YTest=labels(testInd);

%设置随机森林模型的参数
ntrees=100; %决策树数量
rng(1); %设定种子，保证结果可重复性
XTrain=real(XTrain);
XTest=real(XTest);
model=TreeBagger(ntrees,XTrain,YTrain,'Method','classification');

YTestPredict=predict(model,XTest);
YTestPredict=cell2mat(YTestPredict);
YTestPredict=str2num(YTestPredict);

figure;
cm=confusionchart(YTest,YTestPredict);

result=confusionmat(YTest,YTestPredict);