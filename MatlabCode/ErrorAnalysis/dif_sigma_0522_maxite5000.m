%% 1. cp all the data into tmp
%% 2. cp .m file || 3. cp TENSORLAB || 4. cp Nick
%% 3. move to /tmp/ and check path and run!
%% 4. mkdir test_date directory (e.g. test0330)
%% Main code-----------K_changes_0331_VOC_changes_ortho_changes.py---------------------------------------------
%{
outl = zeros(100,10);      %preallocate
for ii = 1:size(ppp,1)    %loop over 100 rows
  ppp(ii,:) = rand(1,10); %insert 10 random numbers
end
%}
clear all;
cd /export/home1/NoCsBack/thesisnl/r0649141/;
%% mkdir test0330; % <----- only once
%% {only once CREATE DIRECTORY
mkdir('./data_gene_0522_range_diff/workspace');
mkdir('./data_gene_0522_range_diff/matlab_output');
mkdir('./data_gene_0522_range_diff/error_per_VOC_TOP');
mkdir('./data_gene_0522_range_diff/convergence_plot');
%% }
addpath('./TENSORLAB')
% addpath('./Thesis/matlab_code/parameter_recover')
addpath('./NIck');
% addpath('./package')
% addpath('./Thesis/NIck');
%VOC_list = [20,30,40,50,60,70];
%VOC_list = [20,30,40,50,60,70];
%VOC_list = [70]
%% -----------------------IF NEED MAKE A LOOP FOR VOC--------------------------
VOC_list = [10];
%% -----------------------IF NEED MAKE A LOOP FOR TOP--------------------------
%TOP_list = [3,4,5,6,7,8,9,10];
TOP_list = [10];
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
%min = [1];
%k = 2;
%count = 0;
N = repelem(j,i);
%% from left to right ex1-->94|| from top to bottom VOC = 20 --> 90

%% Insert #rows == #runs for each comb:VOC x TOP x Orth 
for VOC = VOC_list %VOC = 10 VOC =30 VOC=50
    %count = count + 1;
    mkdir(['./data_gene_0522_range_diff/workspace/VOC' int2str(VOC) '/']);
    mkdir(['./data_gene_0522_range_diff/matlab_output/VOC' int2str(VOC) '/']);    
    for TOP = TOP_list% TOP = 3
        pdir = repelem(round(1/TOP,4),TOP);
        %k = 0
        count = 0;
        out_gcond = zeros(1,length(min));
        out_alpha = zeros(1,length(min));  
        out_word = zeros(1,length(min));
        out_time = zeros(1,length(min));
        out_rel_dif_alha = zeros(1,length(min));
        out_ite = zeros(1,length(min));
        out_sdf_nls_abserr = zeros(1,length(min));
        out_sdf_nls_relerr = zeros(1,length(min));
        for k = min%k =1
            count = count + 1;
            for ii = 2:5%ii =1 
                sum_difalpha = nan(1,1);
                sum_difword = nan(1,1);
                sum_trd_geo = nan(1,1);
                %Taking loop for different dataset

                %for i = DN
                 %   for j = NT
                  %      for k = min

                %count = count + 1;
                %% A = load(['./test_different_VOC_0328/observed_word_matrix_DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
                A = load(['./data_gene_0522_range_dif/observed_word_matrix__K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
                A = A.matrix;
                
                tic;
                [x,y,z,v,w] = para_recover_0131_nopause_maxite5000(A, TOP, sum(pdir), N);
                elapse_time = toc;
                fig = figure('visible','off');
                semilogy(v.relfval);
                legend(['K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho\_minus' num2str(k) 'index\_expe' num2str(ii)']);
                %title(['K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho\_minus' num2str(k) 'index\_expe' num2str(ii)'])
                xlabel('Iterations');
                ylabel('Log(refval)');
                saveas(fig,['./data_gene_0522_range_diff/convergence_plot/K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii)'],'epsc');
                
                
                %Condition number ?
                p = struct2cell(w.factors);
                %trd_geomcond(p)%Infinitity...
                trd_geo = trd_geomcond({x,x,x});
                %save(['matlab_result_run' num2str(j) '.mat']);
                %% load(['./test_different_VOC_0328/parameter_original_DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
                load(['./data_gene_0522_range_dif/parameter_original_K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
                %parameter_original_DN10NT10VOC10ortho_minus7.mat
                generated_para = phi_matrix_original;
                recovered_para = x;
                PW = [generated_para, recovered_para];
                
                %% error for alpha-prior
                %% sum of abso error 
                abs_dif_alpha = sum(abs(pdir - y));
                %% relative errro ?
                rel_dif_alpha = sum(abs(pdir - y))/norm(pdir,2);
                
                %% RELATIVE ERROR --- 'A \ B'-->  Solve systems of linear equations Ax = B for x for topic-word-matrix
%                 fro = nan(1,100);
%                 fro_test = nan(1,100);
%                 fro_modi = nan(1,100);
%                 clearvars -except ii fro fro_test fro_modi
%                 load(['./LDA/test0321/workspace/run' int2str(ii) '/matlab_result_DN1000NT100VOC10ortho_minus2.mat']);
%                 load(['./LDA/test0321/parameter_original_DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
                %% || A_true - A_computed * (A_computed \ A_true) ||_F^2 + || B_true - B_computed * (B_computed \ B_true) ||_F^2 + || C_true - C_computed * (C_computed \ C_true) ||_F^2
                %fro(ii) = sqrt((norm(phi_matrix_original-x*(x \ phi_matrix_original), 'fro'))^2*3);
                rel_dif_word = sqrt(3)*norm( phi_matrix_original - x * (x \ phi_matrix_original), 'fro');
                
                %% RELATIVE ERROR    
                %save to .txt
                %PWt = array2table([generated_para, recovered_para], 'VariableNames',{'generated_para_T1','generated_para_T2','recovered_para_T1', 'recovered_para_T2'});
                %PW_p = array2table([trd_geo, abs_dif_alpha, rel_dif_alpha, rel_dif_word], 'VariableNames',{'trd_geo', 'abs_dif_alpha', 'rel_difalpha', 'rel_dif_alpha, rel_dif_topic_word_dis'});
                clear A;
                clear generated_para;
                clear phi_matrix_original;
                clear recovered_para;
                save(['./data_gene_0522_range_diff/workspace/VOC' int2str(VOC) '/run' int2str(ii) 'matlab_result_K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat'], '-regexp', '^(?!(out_gcond|out_word|out_alpha)$).');
                %mkdir(['./LDA/test0321/matlab_output/run' int2str(ii) '/']);
                %writetable(PWt,['./data_gene_0522_range_diff/matlab_output/VOC' int2str(VOC) '/run' int2str(ii) 'Parameter_RC_vs_ORI_DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.txt'],'Delimiter','\t');
                %writetable(PW_p,['./data_gene_0522_range_diff/matlab_output/VOC' int2str(VOC) '/run' int2str(ii) 'hypara_geocond_summary_value_K' num2mstr(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.txt'],'Delimiter','\t');
                sum_difalpha(1) = abs_dif_alpha;
                sum_difword(1) = rel_dif_word;
                sum_trd_geo(1) = trd_geo;
                %open .txt
                            %FID=fopen('testtxt1.txt');
                            %datacell = textscan(FID, '%f%f%f%f', 'HeaderLines', 10, 'CollectOutput', 1);
                            %fclose(FID);
                            %B = datacell{1};
                    %disp(i)
                  %      end
                 %   end
                %end
                out_gcond(ii,count) = sum_trd_geo;
                out_alpha(ii,count) = sum_difalpha;  
                out_word(ii,count) = sum_difword;
                out_time(ii,count) = elapse_time;
                out_rel_dif_alha(ii, count)  = rel_dif_alpha;
                out_ite (ii,count) = v.iterations;
                out_sdf_nls_abserr = v.abserr;
                out_sdf_nls_relerr = v.relerr;
            end
        end
        save(['./data_gene_0522_range_diff/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat'], 'out_gcond' ,'out_alpha' ,'out_word' ,'out_time' ,'out_rel_dif_alha');
    end
end
