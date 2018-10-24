clear all;
cd /export/home1/NoCsBack/thesisnl/r0649141/;
%% mkdir test0330; % <----- only once
%% {only once CREATE DIRECTORY
%mkdir('./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_BE_error_consider_permutation/workspace');
%mkdir('./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_BE_error_consider_permutation/matlab_output');
%mkdir('./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_BE_error_consider_permutation/error_per_VOC_TOP');
%mkdir('./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_BE_error_consider_permutation/convergence_plot');
%% }
addpath('./TENSORLAB')
% addpath('./Thesis/matlab_code/parameter_recover')
addpath('./NIck');
%% Calcurate Frobinean (error) distance
fro = nan(1,100);
fro_test = nan(1,100);
fro_modi = nan(1,100);

clear;
%TOP_list = [3,4,5,6,7,8,9,10];
TOP_list = [3,4,5,6,7,8,9,10];
VOC = 20;
%VOC_list = [20,30,40,50,60,70];
VOC_list = [60,70];
NT = 100;
j = 100;
%#documents
%DN = [1000,100,10];
DN = 1000;
i = 1000;
%#restriction for orthogonality---> dot product < 10^(-min)
%% min corresponding to the threshould for smallest singular value--ortho = [1.1, 0.9, 0.5, 0.1, 0.05, 0.01]
%% min=[1~5]--> ortho[i]< singular <ortho[i+1], min[6] --> ortho[5] < sing < ortho[6]
min = [1,2,3,4,5,6];
%min = 2;
%k = 2;
%count = 0;
N = repelem(j,i);
k = 6;
x = [1,2,3,4,5,6];
test_cond_gibbs = zeros(10,1);
test_time_gibbs = zeros(10,1);
test_we_gibbs = zeros(10,1);
test_cond_mean = zeros(10,1);
test_time_mean = zeros(10,1);
test_we_mean = zeros(10,1);
test_cond_median = zeros(10,1);
test_time_median = zeros(10,1);
test_we_median = zeros(10,1);

%% for loop wrt VOC, uncomment out for the following two lines
%for VOC = VOC_list
    %TOP = 5
%% for loop wrt TOP, comment out for the above two lines and use
%for TOP = TOP_list


count = 0;
for VOC = VOC_list %VOC=70
%for TOP = TOP_list% TOP = 3
%for ii = 2:2
    TOP = 5;
    pdir = repelem(round(1/TOP,4),TOP);
    %k = 6
    count = count + 1;
    %load(['./Error_Ten_0410/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat'])
    %load(['data_gene_0417_ranini_ITE300_TEN/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    load(['data_gene_0418_ite5000_TolFun24_TolX24_tendec/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %'out_gcond' ,'out_alpha' ,'out_word' ,'out_time' ,'out_rel_dif_alha'
    %load(['./Gibbs0408/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %load(['./Gibbs0409_BI100_ITE200/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
  
    test_cond_mean(count,1) = mean(out_gcond(:,1));
    test_time_mean(count,1) = mean(out_time(:,1));
    test_we_mean(count,1) = mean(out_word(:,1));
    test_cond_median(count,1) = median(out_gcond(:,1));
    test_time_median(count,1) = median(out_time(:,1));
    test_we_median(count,1) = median(out_word(:,1));
    
    %{
    test_cond_gibbs(count,1) = out_gcond_gibbs(1,1);
    test_time_gibbs(count,1) = out_time_gibbs(1,1);
    test_we_gibbs(count,1) = out_word_gibbs(1,1);
    %}
%end
    
    %end
end