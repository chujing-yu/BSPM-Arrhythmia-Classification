function [ecg_filter,Signal_filter_index,Signal_filter]=Find_R_Peaks(signal,fs,flag)
    %% R波定位Pan-Tompkins法
    Fs = fs;
    %ecg是输入的信号，Fs是抽样频率，plot是是否绘图，绘图为1，不绘图为0
    
    
    % 经过滤波求导平方均值后的信号和阈值
    Signal_mean = [];
    Signal_mean_index = [];
    Thresh_Signal_mean = 0;
    Thresh_Noise_mean = 0;
    
    % 经过滤波求导平方均值后的噪声
    Noise_mean = [];
    Noise_mean_index = [];
    
    % 经过带通滤波后的信号和阈值
    Signal_filter = [];
    Signal_filter_index = [];
    Thresh_Signal_filter = 0;
    Thresh_Noise_filter = 0;
    
    % 当前信号、当前噪声大小
    Signal_temp_mean = 0;
    Noise_temp_mean = 0;
    Signal_temp_filter = 0;
    Noise_temp_filter = 0;
    
    % 经过滤波求导平方均值后的信号和噪声的暂存
    Signal_mean_buf = [];
    Noise_mean_buf = [];
    
    
    % 经过带通滤波后的信号和噪声的暂存
    Signal_filter_buf = [];
    Noise_filter_buf = [];
    
    Thresh_Signal_mean_buf = [];
    Thresh_Signal_filter_buf = [];
    
    % 8次有规律的心率间隔
    regular_RR_interval = 0;
    % 8次心率间隔
    mean_interval = 0;
    % 心率间隔，regular_RR_interval有值就是regular_RR_interval，不然就是mean_interval
    RR_interval =  0;
    
    %<0.36Fs时判断是否是T波，如果确实是R波就是0，否则为1
    R_Flag = 0;
    
    %% 滤波1，先低通后高通
    if Fs == 200
        
        %低通，H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2
        
        b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
        a = [1 -2 1];
        
        h_l = filter(b,a, [1,zeros(1,12)]);
        
        ecg_l = conv(signal, h_l);
        ecg_l = ecg_l / max(abs(ecg_l));
        
        
        %高通 H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1))
        
        b = [-1 zeros(1,15) 32 -32 zeros(1,14) 1];
        a = [1 -1];
        
        h_h = filter(b,a, [1, zeros(1,32)]);
        ecg_filter = conv(ecg_l, h_h);
        ecg_filter = ecg_filter / max(abs(ecg_filter));
    
    %% 滤波2，直接buttord滤波，5-15Hz的数字带通滤波器
    else
        wp = [5 15] * 2 / Fs;
        
        [B,A] = butter(3, wp);%3阶滤波器
        
        ecg_filter = signal%filtfilt(B,A, signal);
        ecg_filter = ecg_filter / max(abs(ecg_filter));
        
    end
    
    %% 求导 H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))
        
        h_d = [-1 -2 0 2 1];
        ecg_deri = conv(h_d, ecg_filter);
        ecg_deri = ecg_deri / max(abs(ecg_deri));
    
    %% 平方
        
        ecg_square = ecg_deri .^ 2;
        
    %% 窗口均值化
        Window_width = 0.15 * Fs;
        
        ecg_mean = conv(ecg_square, ones(1, Window_width) / Window_width);
    
    %% 绘图
    if flag
        figure;
        if Fs == 200
        ax(1)=subplot(321);plot(signal);axis tight;title('Raw ECG Signal');
        ax(2)=subplot(322);plot(ecg_l);axis tight;title('After LowerPass');
        else
            ax(1)=subplot(3,2,[1,2]);plot(signal);axis tight;title('Raw ECG Signal');
        end
        
        ax(3) = subplot(323);plot(ecg_filter);axis tight;title('After band filter');
        ax(4) = subplot(324);plot(ecg_deri);axis tight;title('After Derivative');
        ax(5) = subplot(325);plot(ecg_square);axis tight;title('After Square');
        ax(6) = subplot(326);plot(ecg_mean);axis tight;title('After Mean');
        
    end
    
    
    %% 通过findpeaks找到局部的波峰
    
    [peaks, locs] = findpeaks(ecg_mean, 'MINPEAKDISTANCE', round(0.2 * Fs));
    
    
    %% 对阈值进行初始化，使用前2s的数据
    
    
    % 初始化信号阈值和噪声阈值，只使用平方并均值后信号前两秒
    %信号阈值为前两秒最大值的0.33,噪声阈值为前两秒平均值的0.5
    Thresh_Signal_mean = max(ecg_mean(1:2*Fs)) / 3;
    Thresh_Noise_mean = mean(ecg_mean(1:2*Fs)) / 2;
    
    Signal_temp_mean = Thresh_Signal_mean;
    Noise_temp_mean = Thresh_Noise_mean;
    
    % 初始化信号阈值和噪声阈值，只使用带通滤波后信号前两秒
    %信号阈值为前两秒最大值的0.33,噪声阈值为前两秒平均值的0.5
    Thresh_Signal_filter = max(ecg_filter(1:2*Fs)) / 3;
    Thresh_Noise_filter = mean(ecg_filter(1:2*Fs)) / 2;
    
    Signal_temp_filter = Thresh_Signal_filter;
    Noise_temp_filter = Thresh_Noise_filter;
    
    
    %% 对每个波峰进行处理
    for i = 1 : length(peaks)  
        
        %% 找到滤波后信号波峰的位置，index是相对值
        if locs(i) - round(0.15 * Fs) >= 1 && locs(i) <= length(ecg_filter)
            [peak_filter, index_filter] = max( ecg_filter(locs(i) - round(0.15 * Fs) : locs(i)));
            index_filter = locs(i) - round(0.15 * Fs) + index_filter - 1;
        else
            if i == 1 % i = 1是时候index是绝对值
                [peak_filter, index_filter] = max( ecg_filter(1 : locs(i)));
            else %index是相对值
                [peak_filter, index_filter] = max( ecg_filter(locs(i) - round(0.15 * Fs) : end));
                index_filter = locs(i) - round(0.15 * Fs) + index_filter - 1 ;
            end
        end
        %% 判断间隔，如果大于1.66要先处理被遗漏的R波，这样才能保持顺序
        
        %要有9个R波，才有8个间隔，才能判断是否大于1.66
        if length(Signal_mean) >= 9
            RR_diff = diff(Signal_mean_index(end - 8:end));
            mean_interval = mean(RR_diff);
            
            temp_interval = Signal_mean_index(end) - Signal_mean_index(end - 1);
            
            if(temp_interval >= 0.92 * mean_interval || temp_interval <= 1.16 * mean_interval)
                Thresh_Signal_mean = Thresh_Signal_mean / 2;
                Thresh_Signal_filter = Thresh_Signal_filter / 2;
            else
                regular_RR_interval = mean_interval;
            end
        end
        
        if regular_RR_interval
            RR_interval = regular_RR_interval;
        elseif mean_interval
            RR_interval = mean_interval;
        else
            RR_interval = 0;
        end
        
        % 心率间隔有值的话
        if RR_interval
            % 超过了1.66 RR，说明中间有漏检,不应期为0.2ms
            if locs(i) - Signal_mean_index(end) >= 1.66 * RR_interval
                [pk_temp, pk_temp_ind] = max( ecg_mean(Signal_mean_index(end) + round(0.2 * Fs) : locs(i) - round(0.2 * Fs)));
                
                pk_temp_ind = Signal_mean_index(end) + round(0.2 * Fs) + pk_temp_ind - 1;
                
                % 均值信号大于噪声峰值，加到buf里
                if pk_temp >= Thresh_Noise_mean
                    Signal_mean = [Signal_mean pk_temp];
                    Signal_mean_index = [Signal_mean_index pk_temp_ind];
                    
                    %找到滤波信号的位置
                    if pk_temp_ind <= length(ecg_filter)
                        [pk_filter_temp, pk_filter_temp_ind] = max(ecg_filter(pk_temp_ind - round(0.15 * Fs) :  pk_temp_ind));
                    else
                        [pk_filter_temp, pk_filter_temp_ind] = max(ecg_filter(pk_temp_ind - round(0.15 * Fs) :  length(ecg_filter)));
                    end
                    
                    pk_filter_temp_ind = ecg_filter(pk_temp_ind - round(0.15 * Fs) + pk_filter_temp_ind - 1);
                    
                    % 滤波信号也大于噪声阈值，加到buf里，同时sig使用新的更新策略
                    if pk_filter_temp >= Thresh_Noise_filter
                        Signal_filter = [Signal_filter pk_filter_temp];
                        Signal_filter_index = [Signal_filter_index pk_filter_temp_ind];
                        
                        Signal_temp_filter = 0.25 * pk_filter_temp + 0.75 * Signal_temp_filter;
                    end
                        
                    Signal_temp_mean = 0.25 * pk_temp + 0.75 * Signal_temp_mean;
    
                end
            end
        end
        
        %% 当前波峰大于信号阈值
        if peaks(i) >= Thresh_Signal_mean
            
            %要有三个波才能有两个间隔
            if(length(Signal_mean) >= 3)
                
                %如果间隔小于0.36Fs,判断是T波还是噪声
                if locs(i) - Signal_mean_index(end) < round(0.36 * Fs)
                    Slope_1 = mean(diff(ecg_mean(locs(i) - round(0.075 * Fs) : locs(i))));
                    Slope_2 = mean(diff(ecg_mean(Signal_mean_index(end) - round(0.075 * Fs) : Signal_mean_index(end))));
                    
                    %斜率太小了，当前波峰识别为噪声
                    if Slope_1 <= Slope_2 / 2
                        
                        Noise_mean = [Noise_mean peaks(i)];
                        Noise_mean_index = [Noise_mean_index locs(i)];
                        Noise_temp_mean = 0.25 * peaks(i) + 0.75 * Noise_temp_mean;
                        Noise_temp_filter = 0.25 * peak_filter + 0.75 * Noise_temp_filter;
                        R_Flag = 1;%是T波
                    
                    else%斜率够大，确实是T波
                        R_Flag = 0;%确实是T波
                    end
                end
            end
            
            if R_Flag == 0
    
                Signal_mean = [Signal_mean peaks(i)];
                Signal_mean_index = [Signal_mean_index locs(i)];
                Signal_temp_mean = 0.125 * peaks(i) + 0.875 * Signal_temp_mean;
    
                if peak_filter >= Thresh_Signal_filter
    
                    Signal_filter = [Signal_filter peak_filter];
                    Signal_filter_index = [Signal_filter_index index_filter];
                    Signal_temp_filter = 0.125 * peak_filter + 0.875 * Signal_temp_filter;
    
                end
            end
            
                        
        %% 当前波峰大于噪声阈值，但小于信号阈值,此时的波峰为噪声波峰
        elseif peaks(i) < Thresh_Signal_mean && peaks(i) >= Thresh_Noise_mean
            
            Noise_temp_mean = 0.125 * peaks(i) + 0.875 * Noise_temp_mean;
            Noise_temp_filter = 0.125 * peak_filter + 0.875 * Noise_temp_filter;
            
        %% 当前波峰小于噪声阈值,此时的波峰为噪声波峰
        else
            Noise_temp_mean = 0.125 * peaks(i) + 0.875 * Noise_temp_mean;
            Noise_temp_filter = 0.125 * peak_filter + 0.875 * Noise_temp_filter;
            
            Noise_mean = [Noise_mean peaks(i)];
            Noise_mean_index = [Noise_mean_index locs(i)];
            
        end
        
        %% 更新阈值参数
        Thresh_Signal_mean = Noise_temp_mean + 0.25 * (Signal_temp_mean - Noise_temp_mean);
        Thresh_Noise_mean = Thresh_Signal_mean / 2;
        
        Thresh_Signal_filter = Noise_temp_filter + 0.25 * (Signal_temp_filter - Noise_temp_filter);
        Thresh_Noise_filter = Thresh_Signal_filter / 2;
        
        %% 信号当前值存入buf
        Signal_mean_buf = [Signal_mean_buf Signal_temp_mean];
        Noise_mean_buf = [Noise_mean_buf Noise_temp_mean];
        
    
        % 经过带通滤波后的信号和噪声的暂存
        Signal_filter_buf = [Signal_filter_buf Signal_temp_filter];
        Noise_filter_buf = [Noise_filter_buf Noise_temp_filter];
        
        Thresh_Signal_mean_buf = [Thresh_Signal_mean_buf Thresh_Signal_mean];
        Thresh_Signal_filter_buf = [Thresh_Signal_filter_buf Thresh_Signal_filter];
        
        
        % 重置Flag
        R_Flag = 0;
        
    end
    
    
    %% 绘制结果图
    if flag
        figure;
        ax(1) = subplot(211);
        plot(ecg_filter);axis tight;title('Filter Signal');
        hold on,scatter(Signal_filter_index, Signal_filter);
        hold on,plot(locs,Noise_filter_buf,'LineWidth',2,'Linestyle','--','color','k');
        hold on,plot(locs,Signal_filter_buf,'LineWidth',2,'Linestyle','-.','color','r');
        hold on,plot(locs,Thresh_Signal_filter_buf,'LineWidth',2,'Linestyle','-.','color','g');
    
        ax(2) = subplot(212);
        plot(ecg_mean);axis tight;title('Mean Signal');
        hold on,scatter(Signal_mean_index, Signal_mean);
        hold on,plot(locs,Noise_mean_buf,'LineWidth',2,'Linestyle','--','color','k');
        hold on,plot(locs,Signal_mean_buf,'LineWidth',2,'Linestyle','-.','color','r');
        hold on,plot(locs,Thresh_Signal_mean_buf,'LineWidth',2,'Linestyle','-.','color','g');
    end
    %% 绘制结果图
    if flag
        figure;
        plot(ecg_filter);axis tight;title('R波定位结果');
        hold on,scatter(Signal_filter_index, Signal_filter);
    end
end