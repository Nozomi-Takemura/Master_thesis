%% Calcurate Frobinean (error) distance
fro = nan(1,100);
fro_test = nan(1,100);
fro_modi = nan(1,100);

clear;
%TOP_list = [3,4,5,6,7,8,9,10];
TOP_list = [3,4,5,6,7,8,9,10];
VOC = 20;
VOC_list = [20,30,40,50,60,70];
%VOC_list = [20,30,40,50];
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
%% k = 6 means that non-orthogonality --> 0.05< sigma <0.01
%% k = 1 means that non-orthogonality --> 1.1< sigma <0.7
k = 6;
%k = 1
x = [1,2,3,4,5,6];
test_abserr_median = zeros(10,1);
test_abserr_mean = zeros(10,1);
test_cond_mean = zeros(10,1);
test_cpderr_mean = zeros(10,1);
test_relerr_mean = zeros(10,1);
test_cond_median = zeros(10,1);
test_cpderr_median = zeros(10,1);
test_relerr_median = zeros(10,1);
%% ------------------------------------------------------------------------------- %%
%% THIS IS THE LOOP TO GET ERROS FOR DIFFERENT NUMBER OF TOPICS
%% for loop wrt VOC, uncomment out for the following two lines
%for VOC = VOC_list
    %TOP = 5
%% for loop wrt TOP, comment out for the above two lines and use
%% ATTENTION FOR EACH LOOP DON'T FORGET MOVING COUNT!!
count = 0;
for TOP = TOP_list

%% ATTENTION FOR EACH LOOP DON'T FORGET MOVING COUNT!!
%count = 0;
%for VOC = VOC_list
%for TOP = TOP_list% TOP = 3
%for ii = 2:2
    %TOP = 5;
    pdir = repelem(round(1/TOP,4),TOP);
    %k = 0
    count = count + 1;
    %load(['./Error_Ten_0410/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat'])
    %load(['data_gene_0417_ranini_ITE300_TEN/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %load(['data_gene_0418_ite_5000_tendec/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %% USE THIS FOR EMPIRICALL VALUE
    %load(['data_gene_0418_ite5000_TolFun24_TolX24_tendec/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %% USE THIS FOR THEORETICAL VALUE
    load(['data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_theoretical_gcond/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    load(['data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_BE_error_consider_permutation/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %'out_gcond' ,'out_alpha' ,'out_word' ,'out_time' ,'out_rel_dif_alha'
    %load(['./Gibbs0408/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %load(['./Gibbs0409_BI100_ITE200/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %% for empirical value
    %{
    test_cond_mean(count,1) = mean(out_gcond(:,1));
    test_time_mean(count,1) = mean(out_time(:,1));
    test_we_mean(count,1) = mean(out_word(:,1));
    test_cond_median(count,1) = median(out_gcond(:,1));
    test_time_median(count,1) = median(out_time(:,1));
    test_we_median(count,1) = median(out_word(:,1));
    %}
    %% FOR theoretical value
    test_cond_mean(count,1) = mean(out_gcond_theory(:,1));
    test_cond_median(count,1) = median(out_gcond_theory(:,1));
    test_abserr_median(count,1) = median(out_abserr(:,1));
    test_abserr_mean(count,1) = mean(out_abserr(:,1));
    test_cpderr_mean(count,1) = mean(out_cpderr(:,1));
    test_relerr_mean(count,1) = mean(out_relerr(:,1));
    test_cpderr_median(count,1) = median(out_cpderr(:,1));
    test_relerr_median(count,1) = median(out_relerr(:,1));
    %{
    test_cond_gibbs(count,1) = out_gcond_gibbs(1,1);
    test_time_gibbs(count,1) = out_time_gibbs(1,1);
    test_we_gibbs(count,1) = out_word_gibbs(1,1);
    %}
%end
    
    %end
end
%% ------------------------------------------------------------------------------- %%

%% ------------------------------------------------------------------------------- %%
%% THIS IS THE SETTING TO GET ERRORS FOR DIFFERENT DICTIONARY SIZE
%count = 0;
%for TOP = TOP_list

%% ATTENTION FOR EACH LOOP DON'T FORGET MOVING COUNT!!
count = 0;
for VOC = VOC_list
%for TOP = TOP_list% TOP = 3
%for ii = 2:2
    TOP = 5;
    pdir = repelem(round(1/TOP,4),TOP);
    %k = 0
    count = count + 1;
    %load(['./Error_Ten_0410/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat'])
    %load(['data_gene_0417_ranini_ITE300_TEN/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %load(['data_gene_0418_ite_5000_tendec/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %% USE THIS FOR EMPIRICALL VALUE
    %load(['data_gene_0418_ite5000_TolFun24_TolX24_tendec/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %% USE THIS FOR THEORETICAL VALUE
    load(['data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_theoretical_gcond/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    load(['data_gene_0418_ite5000_TolFun24_TolX24_tendec_0520_BE_error_consider_permutation/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %'out_gcond' ,'out_alpha' ,'out_word' ,'out_time' ,'out_rel_dif_alha'
    %load(['./Gibbs0408/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %load(['./Gibbs0409_BI100_ITE200/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat']);
    %% for empirical value
    %{
    test_cond_mean(count,1) = mean(out_gcond(:,1));
    test_time_mean(count,1) = mean(out_time(:,1));
    test_we_mean(count,1) = mean(out_word(:,1));
    test_cond_median(count,1) = median(out_gcond(:,1));
    test_time_median(count,1) = median(out_time(:,1));
    test_we_median(count,1) = median(out_word(:,1));
    %}
    %% FOR theoretical value
    test_cond_mean(count,1) = mean(out_gcond_theory(:,1));
    test_cond_median(count,1) = median(out_gcond_theory(:,1));
    test_abserr_median(count,1) = median(out_abserr(:,1));
    test_abserr_mean(count,1) = mean(out_abserr(:,1));
    test_cpderr_mean(count,1) = mean(out_cpderr(:,1));
    test_relerr_mean(count,1) = mean(out_relerr(:,1));
    test_cpderr_median(count,1) = median(out_cpderr(:,1));
    test_relerr_median(count,1) = median(out_relerr(:,1));
    %{
    test_cond_gibbs(count,1) = out_gcond_gibbs(1,1);
    test_time_gibbs(count,1) = out_time_gibbs(1,1);
    test_we_gibbs(count,1) = out_word_gibbs(1,1);
    %}
%end
    
    %end
end





%% for orthogonality changes - remove the last column
% load Error_summary_VOC20TOP-2to10-DN1000NT100VOC20ortho_minus6.mat—OC1.1-0,.7: Max ITE:300
mean(out_abserr, 1)'
mean(out_cpderr, 1)'
mean(out_relerr, 1)'

median(out_abserr, 1)'
median(out_cpderr, 1)'
median(out_relerr, 1)'

mean(out_gcond_theory, 1)'
median(out_gcond_theory, 1)'


        %'out_gcond_gibbs' ,'out_word_gibbs' ,'out_time_gibbs'
        scatter(log(out_gcond(ii,:)), log(out_word(ii,:)),'r','LineWidth',1)
        hold on;
        scatter(log(out_gcond_gibbs(ii,:)), log(out_word_gibbs(ii,:)),'b','LineWidth',1)
        xlabel('log(condition number)','FontSize',18)
        ylabel('$log\{\sum_{k=1}^{2}\{\sum_{n=1}^{10}|| (\phi_{true})_{n,k} - (\phi_{recovered})_{n,k} ||\}\}$','Interpreter','latex','FontSize',14)
        corrcoef(log(out_gcond(ii,:)), log(out_word(ii,:)))
        corrcoef(log(out_gcond_gibbs(ii,:)), log(out_word_gibbs(ii,:)))
        legend('Tensor decomposition approach', 'Collapsed gibbs sampling approach');
        title(['Log(condition number) vs log(relative errors): #topics =' num2str(TOP)] )
        hold off;
        pause;
        
        scatter(x,log(out_gcond(ii,:)),'r','LineWidth',1)
        hold on;
        scatter(x,log(out_gcond_gibbs(ii,:)),'b','LineWidth',1)
        xlabel('Orthogonality constrain for generated topic-word vectors: righter is less semi-orthogonal','FontSize',14)
        ylabel('log(condition number)')
        %ylabel('$log\{\sum_{k=1}^{2}\{\sum_{n=1}^{10}|| (\phi_{true})_{n,k} - (\phi_{recovered})_{n,k} ||\}\}$','Interpreter','latex','FontSize',18)
        legend('Tensor decomposition approach', 'Collapsed gibbs sampling approach');
        title(['Magnitude of orthogonarility constrain  vs log(Condition number): #topics =' num2str(TOP)])
        hold off;
        pause;
        
        scatter(x,log(out_word(ii,:)),'r','LineWidth',1)
        hold on;
        scatter(x,log(out_word_gibbs(ii,:)),'b','LineWidth',1)
        xlabel('Orthogonality constrain for generated topic-word vectors: righter is less semi-orthogonal','FontSize',14)
        %xlabel('log(condition number)','FontSize',18)
        ylabel('$log\{\sum_{k=1}^{2}\{\sum_{n=1}^{10}|| (\phi_{true})_{n,k} - (\phi_{recovered})_{n,k} ||\}\}$','Interpreter','latex','FontSize',18)
        %corrcoef(log(out_gcond(ii,:)), log(out_word(ii,:)))
        %corrcoef(log(out_gcond_gibbs(ii,:)), log(out_word_gibbs(ii,:)))
        legend('Tensor decomposition approach', 'Collapsed gibbs sampling approach');
        title(['Magnitude of orthogonarility constrain vs log(relative errors): #topics =' num2str(TOP)])
        hold off;
        pause;
        
        scatter(x,log(out_time(ii,:)),'r','LineWidth',1)
        hold on;
        scatter(x,log(out_time_gibbs(ii,:)),'b','LineWidth',1)
        xlabel('Orthogonality constrain for generated topic-word vectors: righter is less semi-orthogonal','FontSize',14)
        %xlabel('log(condition number)','FontSize',18)
        %ylabel('$log\{\sum_{k=1}^{2}\{\sum_{n=1}^{10}|| (\phi_{true})_{n,k} - (\phi_{recovered})_{n,k} ||\}\}$','Interpreter','latex','FontSize',18)
        ylabel('Computational time')
        %corrcoef(log(out_gcond(ii,:)), log(out_word(ii,:)))
        %corrcoef(log(out_gcond_gibbs(ii,:)), log(out_word_gibbs(ii,:)))
        legend('Tensor decomposition approach', 'Collapsed gibbs sampling approach');
        title(['Magnitude of orthogonarility constrain vs log(computational time): #topics =' num2str(TOP)])
        hold off;
        pause;
    end
end


 'out_gcond' ,'out_alpha' ,'out_word' ,'out_time' ,'out_rel_dif_alha'


 
 
 
 
 
 
 
 for ii = 1:100
    clearvars -except ii fro fro_test fro_modi
    load(['~/PYTHON/LDA/test0304/workspace/run' int2str(ii) '/matlab_result_DN1000NT100VOC10ortho_minus2.mat']);
    load(['~/PYTHON/LDA/test0304/parameter_original_DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
    %% || A_true - A_computed * (A_computed \ A_true) ||_F^2 + || B_true - B_computed * (B_computed \ B_true) ||_F^2 + || C_true - C_computed * (C_computed \ C_true) ||_F^2
    %fro(ii) = sqrt((norm(phi_matrix_original-x*(x \ phi_matrix_original), 'fro'))^2*3);
    fro(ii) = sqrt(3)*norm( phi_matrix_original - x * (x \ phi_matrix_original), 'fro');
    generated_para = phi_matrix_original;
    recovered_para = x;
    PW = [generated_para, recovered_para];
    dif1 = sum(abs(PW(:,1) - PW(:,3)));
    dif4 = sum(abs(PW(:,2) - PW(:,4)));
    fro_test(ii) = dif1 + dif4;
    dif2 = sum(abs(PW(:,1) - PW(:,4)));
    if dif1 > dif2
        recovered_para = [recovered_para(:,2),recovered_para(:,1)];
    end
    fro_modi(ii) = sqrt(3)*norm( phi_matrix_original - recovered_para * (recovered_para \ phi_matrix_original), 'fro');
    
end

%%  Remove outlier
clearvars -except fro fro_test fro_modi
load('~/PYTHON/LDA/test0304/after100runs.mat');
nonoutlier_fro = find(fro < 1);%index of outlier
%nonoutlier_fro = find(fro < 0.008);%index of outlier
sum(fro)
%% remove outlier---for "fro"
fro_nonout = fro(nonoutlier_fro);
fro_nonout = fro_nonout';

%% remove outliers for "fro_test" -- not fnorm bu just sum of absolute difference, ignoring permutation
nonoutlier_fro_test = find(fro_test < 1);%index of outlier
%% remove outlier---for "fro"
fro_test_nonout = fro_test(nonoutlier_fro_test);
fro_test_nonout = fro_test_nonout';

%% remove outlier-- for "fro_modi"
%clear nonoutlier_fro_modi fro_modi_nonout
nonoutlier_fro_modi = find(fro_modi < 1);%index of outlier
%% remove outlier---for "fro_modi"
fro_modi_nonout = fro_modi(nonoutlier_fro_modi);
fro_modi_nonout = fro_modi_nonout';


%% remove outliers -- for "out_word"
nonoutlier = find(out_word < 1);%index of outlier

word_nonout = out_word(nonoutlier);
gcond_nonout = out_gcond(nonoutlier);
alpha_nonout = out_alpha(nonoutlier);

%% check if "fro" and "out_word"'s outliers're same or not
isequal(nonoutlier,nonoutlier_fro)% return 1 == same
isequal(nonoutlier,nonoutlier_fro_modi)
%% correlation coefficient -- fro_test
corrcoef(log(out_gcond(nonoutlier_fro_test)), log(fro_test_nonout))
%{
    1.0000    0.8402
    0.8402    1.0000
%}
%% correlation coefficient -- fro
corrcoef(log(out_gcond(nonoutlier_fro)), log(fro_nonout))
%{
    1.0000    0.1900
    0.1900    1.0000
%}
%% correlation coefficient -- fro_modi
corrcoef(log(out_gcond(nonoutlier_fro_modi)), log(fro_modi_nonout))
%{
    1.0000    0.1900
    0.1900    1.0000
%}
%% plot log_fro vs log_gcond
scatter(log(out_gcond(nonoutlier_fro)), log(fro_nonout))
hold on;
xlabel('log(condition number)','FontSize',18)
ylabel('$log\{\sum_{k=1}^{2}\{\sum_{n=1}^{10}|| (\phi_{true})_{n,k} - (\phi_{recovered})_{n,k} ||\}\}$','Interpreter','latex','FontSize',18)
corrcoef(log(gcond_nonout), log(word_nonout))
%{
    1.0000    0.8050
    0.8050    1.0000
%}
dim = [.5 .09 .7 .3];
str = 'correlation coefficent: 0.8050';
a = annotation('textbox',dim,'String',str,'FitBoxToText','on');
a.FontSize = 14;