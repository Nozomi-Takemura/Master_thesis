%% from left to right ex1-->94|| from top to bottom VOC = 20 --> 90
clear all;
cd /export/home1/NoCsBack/thesisnl/r0649141/;
%% mkdir test0330; % <----- only once
%% {only once CREATE DIRECTORY
mkdir('./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_theoretical_gcond/workspace');
mkdir('./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_theoretical_gcond/matlab_output');
mkdir('./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_theoretical_gcond/error_per_VOC_TOP');
mkdir('./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_theoretical_gcond/convergence_plot');
%% }
addpath('./TENSORLAB')
% addpath('./Thesis/matlab_code/parameter_recover')
addpath('./NIck');
% addpath('./package')
% addpath('./Thesis/NIck');
%VOC_list = [20,30,40,50,60,70];
VOC_list = [20,30,40,50,60,70];
TOP_list = [3,4,5,6,7,8,9,10];
%TOP_list = [9,10];
%mkdir(['./LDA/test0321/matlab_output/run' int2str(ii) '/']);
%#terms/doc
%NT = [1000,100,10];
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

%% Insert #rows == #runs for each comb:VOC x TOP x Orth 
for VOC = VOC_list %VOC = 10 VOC =30 VOC=50
    %count = count + 1;
    mkdir(['./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_theoretical_gcond/workspace/VOC' int2str(VOC) '/']);
    mkdir(['./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_theoretical_gcond/matlab_output/VOC' int2str(VOC) '/']);    
    for TOP = TOP_list% TOP = 3
        pdir = repelem(round(1/TOP,4),TOP);
        %k = 0
        count = 0;
        out_gcond_theory = zeros(1,length(min));
        for k = min%k =1
            count = count + 1;
            for ii = 1:5%ii =1 
                sum_difalpha = nan(1,1);
                sum_difword = nan(1,1);
                sum_trd_geo = nan(1,1);
                %Taking loop for different dataset
                %Condition number ?
                load(['./data_gene_0407/parameter_original_K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
                %% Calculate Gcond for theoretical topic-word distribution
                trd_geo_theory = trd_geomcond({phi_matrix_original,phi_matrix_original,phi_matrix_original});
                save(['./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_theoretical_gcond/workspace/VOC' int2str(VOC) '/run' int2str(ii) 'matlab_result_K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat'], '-regexp', '^(?!(out_gcond_theory)$).');
                out_gcond_theory(ii,count) = trd_geo_theory;
            end
        end
        save(['./data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_theoretical_gcond/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat'], 'out_gcond_theory');
    end
end