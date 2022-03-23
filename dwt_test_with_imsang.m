clear
close all
clc

load('C:\Lee\DWT\seq11_testset')
load('C:\Lee\Zing Seqence 결정 관련\밸리데이션 데이터\1차\170526 AgaMatrix Z!ng 신규 알고리즘 검증 내부임상\Zing_validation')


bi=seq1_validation;

seq1_c_t23.data=data;
seq1_c_t23.index=index;

tic
for k=1:length(bi.data(1,:))

    bi1_data=bi.data(:,k);
    
    
    sample_leng=length(index);
    %% 구간 끊기 없이
    pair=zeros(1,sample_leng);
    for i=1:sample_leng
        [pair(i),mapp]=DTW([bi1_data(172:450,1);bi1_data(451:10:971,1);bi1_data(972:1210,1)],[seq1_c_t23.data(172:450,i);seq1_c_t23.data(451:10:971,i);seq1_c_t23.data(972:1210,i)]);
        total_map(:,:,i)=mapp;
    end
    
    A=[pair;seq1_c_t23.index';(1:1:sample_leng);repmat(k,1,sample_leng)]';
    vali_info=[k [bi.ysi(k) bi.hct(k) bi.temp(k)]];
    rankk=sortrows(A,1);
    mini30=rankk(1:10,:);
    
    f_vali_info=[f_vali_info;vali_info];
    f_rank_value=[f_rank_value;mini30];
    dtw_map(:,:,k)=total_map(:,:,mini30(1,5));
    
    disp(k)
    
end