function [seg_fft_SM,seg_fft_SSD,seg_fft_SD,seg_fft_SK]=Seg_Freq_Charact(seg_ecg_filter,seg_len)
    %% 频域特征提取
    seg_ecg_fft=fft(seg_ecg_filter);
    %频谱均值SM
    seg_fft_SM=sum(seg_ecg_fft,2)/seg_len;
    %频谱标准差SSD
    seg_fft_SSD=std(seg_ecg_fft,0,2);
    %频谱偏度SD
    seg_fft_SD=sum((seg_ecg_fft-seg_fft_SM),2)/power(seg_fft_SSD,1.5);
    seg_fft_SD=seg_fft_SD/seg_len;
    %频谱峭度SK
    seg_fft_SK=sum(power((seg_ecg_fft-seg_fft_SM),4),2)/power(seg_fft_SSD,2);
   seg_fft_SK=seg_fft_SK/seg_len;
end