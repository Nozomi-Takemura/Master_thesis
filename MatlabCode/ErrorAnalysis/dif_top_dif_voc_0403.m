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
mkdir('./test0403/workspace');
mkdir('./test0403/matlab_output')
mkdir('./test0403/error_per_VOC_TOP')
%% }
addpath('./TENSORLAB')
% addpath('./Thesis/matlab_code/parameter_recover')
addpath('./NIck');
% addpath('./package')
% addpath('./Thesis/NIck');
VOC_list = [20,30,40,50,60,70,80,90];
TOP_list = [3,4,5,6,7];
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
%% from left to right ex1-->94|| from top to bottom VOC = 20 --> 90

%% Insert #rows == #runs for each comb:VOC x TOP x Orth 
for VOC = VOC_list %VOC = 10
    %count = count + 1;
    mkdir(['./test0403/workspace/VOC' int2str(VOC) '/']);
    mkdir(['./test0403/matlab_output/VOC' int2str(VOC) '/']);    
    for TOP = TOP_list% TOP = 3
        pdir = repelem(round(1/TOP,4),TOP);
        %k = 0
        count = 0;
        out_gcond = zeros(1,length(min));
        out_alpha = zeros(1,length(min));  
        out_word = zeros(1,length(min));
        out_time = zeros(1,length(min));
        out_rel_dif_alha = zeros(1,length(min));
        for k = min%k =1
            count = count + 1;
            for ii = 1:size(out_gcond,1)%ii =1 
                %clear all;
                % specific clear
                %clearvars -except j
                %path to tensorlab
                %clear all;
                %cd '~/Thesis/TENSORLAB';
                %addpath('~/Thesis/TENSORLAB')
                %addpath(pwd)
                %path to function para_recover
                %cd '~/Thesis/matlab_code/parameter_recover';
                %addpath(pwd);
                %addpath('~/Thesis/matlab_code/parameter_recover')
                %cd '~/PYTHON/LDA/test0201/ND1200NT1200/';
                %cd '~/PYTHON/LDA/test0201'/ND600NT600/;
                %cd '~/PYTHON/LDA/test0208/K2ND50TE50/'
                %cd '~/PYTHON/LDA/test0217/K2ND1000TE1000VOC10'
                %cd '~/PYTHON/LDA/test0217/K2ND1000TE100VOC10'
                %cd '~/PYTHON/LDA/test0225';
                %addpath('~\Documents\thesis\matlab_code\parameter_recover')
                %for trd_geo function
                %addpath('~/Thesis/NIck');
                %choose #words in documents * #distict words in corpus observed-word-matrix of unit vectors  
                %A = uiimport;
                %A = load('observed_word_matrix_DN1200NT1200.mat');
                %A = load('observed_word_matrix_DN600NT600.mat');
                %A = load('observed_word_matrix_DN1000NT100AR1.0.mat');

                %true_hyper_for dirichlet
                %pdir = [0.5, 0.5];
                sum_difalpha = nan(1,1);
                sum_difword = nan(1,1);
                sum_trd_geo = nan(1,1);
                %Taking loop for different dataset

                %for i = DN
                 %   for j = NT
                  %      for k = min

                %count = count + 1;
                %% A = load(['./test_different_VOC_0328/observed_word_matrix_DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
                A = load(['./r0649141_test0402/observed_word_matrix__K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
                A = A.matrix;
                %a = sum(sum(A));   %test 
                %disp(a)    %test
                %A = A.df;

                % A <= observed-word matrix of unit-vector, 2<= rank, 1 <= alpha_{0} of sum
                % of parameters of Dirichlet

                % A vector of # of terms contained in each document. ---% Specification of 
                % the vector N for #doc=600 & #terms in each doc =600 is shown below.
                %N` = [NT];
                tic;
                [x,y,z,v,w] = para_recover_0131_nopause_maxite300(A, TOP, sum(pdir), N);
                elapse_time = toc; 
                %Condition number ?
                p = struct2cell(w.factors);
                %trd_geomcond(p)%Infinitity...
                trd_geo = trd_geomcond({x,x,x});
                %save(['matlab_result_run' num2str(j) '.mat']);
                %% load(['./test_different_VOC_0328/parameter_original_DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
                load(['./r0649141_test0402/parameter_original_K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
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
                save(['./test0403/workspace/VOC' int2str(VOC) '/run' int2str(ii) 'matlab_result_K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat'], '-regexp', '^(?!(out_gcond|out_word|out_alpha)$).');
                %mkdir(['./LDA/test0321/matlab_output/run' int2str(ii) '/']);
                %writetable(PWt,['./test0403/matlab_output/VOC' int2str(VOC) '/run' int2str(ii) 'Parameter_RC_vs_ORI_DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.txt'],'Delimiter','\t');
                %writetable(PW_p,['./test0403/matlab_output/VOC' int2str(VOC) '/run' int2str(ii) 'hypara_geocond_summary_value_K' num2mstr(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.txt'],'Delimiter','\t');
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
            end
        end
        save(['./test0403/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat'], 'out_gcond' ,'out_alpha' ,'out_word' ,'out_time' ,'out_rel_dif_alha');
    end
end
save('./test0403/after90runs_times_9voc.mat');
%% Until here--------------------------------------------------------

%std for each column
gcond_std = std(out_gcond);
alpha_std = std(out_alpha);
word_std = std(out_word);
%average for each column
gcond_mean = mean(out_gcond);
alpha_mean = mean(out_alpha);
word_mean = mean(out_word);
%median for each column
gcond_median = median(out_gcond);
alpha_median = median(out_alpha);
word_median = median(out_word);

    NT = [1000,100,10];
    %#documents
    DN = [1000,100,10];
xaxis = [1, 2, 3, 4, 5, 6, 7];
plot(xaxis, sum_difalpha(1:7))
hold;
plot(xaxis, sum_difalpha(8:14), 'r');
plot(xaxis, sum_difalpha(15:21), 'g');
plot(xaxis, sum_difalpha(22:28), 'c');
plot(xaxis, sum_difalpha(29:35), 'm');
plot(xaxis, sum_difalpha(36:42), 'y');
plot(xaxis, sum_difalpha(43:49), 'k');
plot(xaxis, sum_difalpha(50:56), ':');
plot(xaxis, sum_difalpha(57:63), '--');
legend('sample size:1000000','sample size:100000', 'sample size:10000','sample size:100000', 'sample size:10000', 'sample size:1000','sample size:10000','sample size:1000','sample size:100');
title('Sum of absolute difference of "alpha" parameter vecotr ')
hold off;
% plot for conditional number 
plot(xaxis, sum_trd_geo(1:7))
hold;
plot(xaxis, sum_trd_geo(8:14), 'r');
plot(xaxis, sum_trd_geo(15:21), 'g');
plot(xaxis, sum_trd_geo(22:28), 'c');
plot(xaxis, sum_trd_geo(29:35), 'm');
plot(xaxis, sum_trd_geo(36:42), 'y');
plot(xaxis, sum_trd_geo(43:49), 'k');
plot(xaxis, sum_trd_geo(50:56), ':');
plot(xaxis, sum_trd_geo(57:63), '--');
legend('sample size:1000000','sample size:100000', 'sample size:10000','sample size:100000', 'sample size:10000', 'sample size:1000','sample size:10000','sample size:1000','sample size:100');
title('Conditional value ?')
hold off;

%plot for word-parameter matrix
plot(xaxis, sum_difword(1:7))
hold;
plot(xaxis, sum_difword(8:14), 'r');
plot(xaxis, sum_difword(15:21), 'g');
plot(xaxis, sum_difword(22:28), 'c');
plot(xaxis, sum_difword(29:35), 'm');
plot(xaxis, sum_difword(36:42), 'y');
plot(xaxis, sum_difword(43:49), 'k');
plot(xaxis, sum_difword(50:56), ':');
plot(xaxis, sum_difword(57:63), '--'); 
%legend('dot product < 10^(-1)','dot product < 10^(-2)', 'dot product < 10^(-3)','dot product < 10^(-4)', 'dot product < 10^(-5)', 'dot product < 10^(-6)','dot product < 10^(-7)');
legend('sample size:1000000','sample size:100000', 'sample size:10000','sample size:100000', 'sample size:10000', 'sample size:1000','sample size:10000','sample size:1000','sample size:100');
title('Sum of absolute difference of "word_topic_matrix" parameter vecotr ')
hold off;
