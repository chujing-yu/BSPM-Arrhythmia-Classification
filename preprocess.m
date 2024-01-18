function signal=preprocess(x_signal,fs)
    %巴特沃斯低通滤波
    Wc=2*45/fs;           %截止频率 50Hz
    [b,a]=butter(4,Wc);
    ecg_signal=filter(b,a,x_signal);

    %% 小波变换＋去基线漂移
    %小波分解
    [c,l]=wavedec(ecg_signal,4,'db3');      %db3小波基函数，分解为4层
    %从[c,l]中提取第N层近似系数（CA）和细节系数（DA）
    ca4=appcoef(c,l,'db3',4);
    cd1=appcoef(c,l,'db3',1);
    cd2=appcoef(c,l,'db3',2);
    cd3=appcoef(c,l,'db3',3);
    cd4=appcoef(c,l,'db3',4);
    %使用stein的无偏似然估计原理选择各层的阈值
    %'rigrsure’为无偏似然估计阈值类型
    thr1=thselect(cd1,'rigrsure');
    thr2=thselect(cd2,'rigrsure');
    thr3=thselect(cd3,'rigrsure');
    thr4=thselect(cd4,'rigrsure');
    %各层的阈值
    TR=[thr1,thr2,thr3,thr4];
    %'s'为软阈值;'h'硬阈值。
    SORH='s';
    %去噪
    %XC为去噪后信号
    %[CXC,LXC]为 小波分解结构
    %PERF0和PERF2是恢复和压缩的范数百分比。
    %'lvd'为允许设置各层的阈值,
    %'gbl'为固定阈值。
    %3为阈值的长度
     [XC,CXC,LXC,PERF0,PERF2]=wdencmp('lvd',ecg_signal,'db3',4,TR,SORH);
    
    maxlevel=7;
    wavename='db6';
    [C,L]=wavedec(XC,maxlevel,wavename);
    
    A7=appcoef(C,L,wavename,7);
    D1=detcoef(C,L,1);
    D2=detcoef(C,L,2);
    D3=detcoef(C,L,3);
    D4=detcoef(C,L,4);
    D5=detcoef(C,L,5);
    D6=detcoef(C,L,6);
    D7=detcoef(C,L,7);
    
    D1=zeros(1,length(D1)); %去高频噪声
    D2=zeros(1,length(D2));
    D3=zeros(1,length(D3));
    A7=zeros(1,length(A7));
    C2=[A7,D7,D6,D5,D4,D3,D2,D1];
    signal=(waverec(C2,L,wavename));  %重构信号
    signal=signal(1000:end);
end 