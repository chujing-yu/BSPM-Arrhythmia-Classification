function [seg_SVD,seg_MAX,seg_STD]=Seg_WPD_Charact(seg_ecg_filter,wavename,decompose_level,entropy_name)
    %小波包分解
    wpd_tree=wpdec(seg_ecg_filter,decompose_level,wavename,entropy_name);
    n=power(2,decompose_level)-1;
    wpd_coef_layer=[];   %最后一层的小波包分解系数
    for index=0:n
        coef=wpcoef(wpd_tree,[decompose_level index]);
        wpd_coef_layer=[wpd_coef_layer;coef];
    end

    seg_SVD=svd(wpd_coef_layer);
    seg_SVD=seg_SVD';
    seg_MAX=max(wpd_coef_layer,[],2);
    seg_MAX=seg_MAX';
    seg_STD=std(wpd_coef_layer,0,2);
    seg_STD=seg_STD';
   
   