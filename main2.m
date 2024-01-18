clear 
clc
fs=500;          %采样频率
% 设置文件夹路径
folderPath = 'E:\新华医院202107心电信号\dataset';
%信号长度改成10S，MATLAB自带的SVM等分类器分类
% 获取文件夹中的所有文件
fileList = dir(fullfile(folderPath, '*.dat'));
%特征矩阵
%时域特征
Amp=[];    %R峰振幅
SD=[];     %标准差
CF=[];     %峰值因子
SF=[];     %形状因子
IF=[];     %脉冲因子
%频域特征
fft_SM=[]; %频谱均值
fft_SSD=[]; %频谱标准差
fft_SD=[];  %频谱偏度
fft_SK=[];  %频谱峭度

%小波域特征
wSVD=[];  %小波系数奇异值
wMAX=[];  %小波系数最大值
wSTD=[];  %小波系数标准差

%每位患者取的信号段数
segment_num=12;
%% 

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
    %每位患者取1段10s的信号
    preprocessed_signal=[];
    segment_start=1;
    segment_len=5000;
    segment_interval=2000;

    %每位患者的Amp,SD,CF,SF,IF特征
    Amp_segs=zeros(1,segment_num);
    SD_segs=zeros(1,segment_num);
    CF_segs=zeros(1,segment_num);
    SF_segs=zeros(1,segment_num);
    IF_segs=zeros(1,segment_num);
    
    %每位患者的SM,SSD,SD,SK特征
    fft_SM_segs=zeros(1,segment_num);
    fft_SSD_segs=zeros(1,segment_num);
    fft_SD_segs=zeros(1,segment_num);
    fft_SK_segs=zeros(1,segment_num);

    %每位患者的小波包SVD,MAX,STD特征
    decompose_level=4;
    %WPD_Charact_type=3;
    WPD_Charact_num=power(2,decompose_level);
    wSVD_segs=zeros(segment_num,WPD_Charact_num);
    wMAX_segs=zeros(segment_num,WPD_Charact_num);
    wSTD_segs=zeros(segment_num,WPD_Charact_num);

    for n=1:segment_num
        segment_start_present=segment_start+(n-1)*(segment_len+segment_interval)
        segment_signal=x_signal(1,segment_start_present:(segment_start_present+segment_len-1));
        preprocessed_segment_signal=preprocess(segment_signal,fs);
        preprocessed_signal=[preprocessed_signal;preprocessed_segment_signal]; %每一行是一段10s的信号

        flag=0;
        [segment_ecg_filter,segment_Signal_filter_index,segment_Signal_filter]=Find_R_Peaks(preprocessed_segment_signal,fs,flag);
        [seg_Amp,seg_SD,seg_CF,seg_SF,seg_IF]=Seg_Time_Charact(segment_ecg_filter,segment_Signal_filter,segment_len);
        Amp_segs(n)=seg_Amp;
        SD_segs(n)=seg_SD;
        CF_segs(n)=seg_CF;
        SF_segs(n)=seg_SF;
        IF_segs(n)=seg_IF;

        [seg_fft_SM,seg_fft_SSD,seg_fft_SD,seg_fft_SK]=Seg_Freq_Charact(segment_ecg_filter,segment_len);
        fft_SM_segs(n)=seg_fft_SM;
        fft_SSD_segs(n)=seg_fft_SSD;
        fft_SD_segs(n)=seg_fft_SD;
        fft_SK_segs(n)=seg_fft_SK;
        
        wavename='db6';
        entropy_name='shannon';
        [seg_SVD,seg_MAX,seg_STD]=Seg_WPD_Charact(segment_ecg_filter,wavename,decompose_level,entropy_name);
        wSVD_segs(n,:)=seg_SVD;
        wMAX_segs(n,:)=seg_MAX;
        wSTD_segs(n,:)=seg_STD;
    end

    %% 时频域、小波域特征
    %每行一个患者
    %时域特征
    Amp=[Amp;Amp_segs];    %R峰振幅
    SD=[SD;SD_segs];     %标准差
    CF=[CF;CF_segs];     %峰值因子
    SF=[SF;SF_segs];     %形状因子
    IF=[IF;IF_segs];     %脉冲因子

    %频域特征
    fft_SM=[fft_SM;fft_SM_segs];
    fft_SSD=[fft_SSD;fft_SSD_segs];
    fft_SD=[fft_SD;fft_SD_segs];
    fft_SK=[fft_SK;fft_SK_segs];

    %小波域特征
    wSVD_segs=reshape(wSVD_segs,1,[]);
    wMAX_segs=reshape(wMAX_segs,1,[]);
    wSTD_segs=reshape(wSTD_segs,1,[]);
    wSVD=[wSVD;wSVD_segs];
    wMAX=[wMAX;wMAX_segs];
    wSTD=[wSTD;wSTD_segs];

end
%% 保存提取的特征
save('Amp.mat','Amp','-mat');
save('SD.mat','SD','-mat');
save('CF.mat','CF','-mat');
save('SF.mat','SF','-mat');
save('IF.mat','IF','-mat');
save('fft_SM.mat','fft_SM','-mat');
save('fft_SSD.mat','fft_SSD','-mat');
save('fft_SD.mat','fft_SD','-mat');
save('fft_SK.mat','fft_SK','-mat');
save('wSVD.mat','wSVD','-mat');
save('wMAX.mat','wMAX','-mat');
save('wSTD.mat','wSTD','-mat');

%% 归一化
Amp=importdata('Amp.mat');
SD=importdata('SD.mat');
CF=importdata('CF.mat');
SF=importdata('SF.mat');
IF=importdata('IF.mat');
fft_SM=importdata('fft_SM.mat');
fft_SSD=importdata('fft_SSD.mat');
fft_SD=importdata('fft_SD.mat');
fft_SK=importdata('fft_SK.mat');
wSVD=importdata('wSVD.mat');
wMAX=importdata('wMAX.mat');
wSTD=importdata('wSTD.mat');

sample_size=length(fileList);      %总的样本数量
charact_type=57;                   %总的特征数量 9+16*3
%分类标签
labels=importdata('Labels.mat');
labels=labels';
%dataset
dataset_row_num=sample_size;   %行数为患者数
dataset_col_num=charact_type*segment_num;  %列数为特征种类*信号段数
dataset=[Amp SD CF SF IF fft_SM fft_SK fft_SD fft_SSD wSTD wMAX wSVD];
normalized_dataset=zeros(dataset_row_num,dataset_col_num); %归一化的数据集
for k=1:dataset_col_num
    col_mean=mean(dataset(:,k)); %列均值
    col_max=max(dataset(:,k));   %列最大值
    col_min=min(dataset(:,k));   %列最小值
    normalized_dataset(:,k)=(dataset(:,k)-col_mean)./(col_max-col_min);
end
%% 作图
numRows=size(Amp,1);
numCols=size(Amp,2);
colors=jet(numRows);

figure;
hold on;

for n=1:numRows
    x=1:numCols;
    y=Amp(n,:);
    scatter(x,y,[],colors(n,:),'filled');
end

figure;
hold on;
for n=1:numRows
    x=1:numCols;
    y=fft_SK(n,:);
    scatter(x,y,[],colors(n,:),'filled');
end

figure;
hold on;
for n=1:numRows
    x=1:numCols;
    y=fft_SM(n,:);
    scatter(x,y,[],colors(n,:),'filled');
end
%% SVM
%从中抽取0.8作为训练集
test_set_num=floor(dataset_row_num*0.3);
x_test=zeros(test_set_num,dataset_col_num);
[x_test,test_index]=datasample(normalized_dataset,test_set_num,1,'Replace',false);
y_test=labels(test_index,:)
x_train=normalized_dataset;
x_train(test_index,:)=[];
y_train=labels;
y_train(test_index,:)=[];

% x_test=normalized_dataset;
% y_test=labels;

x_train=real(x_train);
x_test=real(x_test);
% 特征归一化
% label采用二进制编码
[y_predict,models]=MultiSvm(x_train,y_train,x_test);

figure;
cm=confusionchart(y_test,y_predict);

% [results,~]=confusionmat(y_test,y_predict);
% % 计算-1类的评价值
% c1_precise = results(1,1)/(results(1,1) + results(2,1));
% c1_recall = results(1,1)/(results(1,1) + results(1,2));
% c1_F1 = 2 * c1_precise * c1_recall/(c1_precise + c1_recall);
% 
% % 计算1类的评价值
% c2_precise = results(2,2)/(results(1,2) + results(2,2));
% c2_recall = results(2,2)/(results(2,1) + results(2,2));
% c2_F1 = 2 * c2_precise * c2_recall/(c2_precise + c2_recall);
%% 随机森林
%拆分为训练集和测试集
trainRatio=0.8;
[trainInd,~,testInd]=dividerand(size(normalized_dataset,1),trainRatio,0,1);
XTrain=normalized_dataset(trainInd,:);
YTrain=labels(trainInd);
% XTest=normalized_dataset(testInd,:);
% YTest=labels(testInd);
XTest=normalized_dataset;
YTest=labels;
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