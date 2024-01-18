function [seg_Amp,seg_SD,seg_CF,seg_SF,seg_IF]=Seg_Time_Charact(seg_ecg_filter,seg_Signal_filter,seg_len)
    %% 时域特征提取
    
    %每个信号段的R峰振幅
    seg_Amp=max(seg_Signal_filter); %取每段信号最大峰值为R峰幅值
    %每个信号段的均值
    seg_mean=sum(seg_ecg_filter,2)./seg_len;
    %每段信号的标准差
    seg_SD=std(seg_ecg_filter,0,2); 
    %峰值因子CF
    seg_RMS=rms(seg_ecg_filter,2)    %每段信号的均方根pieces_RMS
    seg_CF=abs(seg_Amp)/seg_RMS;
    %形状因子SF
    seg_abs_mean=sum(abs(seg_ecg_filter),2)/seg_len;
    seg_SF=seg_RMS/seg_abs_mean;
    %脉冲因子IF
    seg_IF=seg_Amp./seg_abs_mean;
end