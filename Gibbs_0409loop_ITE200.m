%% Data modification
clear
cd /export/home1/NoCsBack/thesisnl/r0649141/;
%% mkdir('Gibbs0409_BI100_ITE200') %only once CREATE DIRECTORY
mkdir('./Gibbs0409_BI100_ITE200/workspace_Gibbs_BI100_ITE200');
mkdir('./Gibbs0409_BI100_ITE200/matlab_output_Gibbs')
%% create directory only once
%% add path
%addpath('./LDA/LDA_Gibbs')
%addpath('/home/nozomi/PYTHON/TenDecLDA/Thesis/LDA_Gibbs/')
%addpath('/home/nozomi/PYTHON/TenDecLDA/Thesis/NIck')
%addpath('./LDA/test0321/')
addpath('/home/r0649141/matlab_code/')
addpath('/home/r0649141/NIck/')
%addpath('/home/r0649141/test0327/')
%% load('./LDA/test0321/after100runs_0321.mat')
%VOC = 10;
NT = 100;
j = 100;
DN = 1000;
i = 1000;
%min = 2;
%k = 2;
count = 0;%count = 1;
N = repelem(j,i);
%VOC_list = [20,30,40,50,60,70,80,90,100];
VOC_list = [20,30,40,50,60,70];
TOP_list = [3,4,5,6,7,8,9,10];
min = [1,2,3,4,5,6];
TOT_RUN = 1;
mkdir('./Gibbs0409_BI100_ITE200/error_per_VOC_TOP')

for VOC = VOC_list%VOC=20
    count = count + 1;
    mkdir(['./Gibbs0409_BI100_ITE200/workspace_Gibbs_BI100_ITE200/VOC' int2str(VOC) '/']);
    mkdir(['./Gibbs0409_BI100_ITE200/matlab_output_Gibbs/VOC' int2str(VOC) '/']);

    for TOP = TOP_list% TOP = 3
        pdir = repelem(round(1/TOP,4),TOP);
        %k = 0
        count = 0;
        out_gcond_gibbs = zeros(TOT_RUN,length(min));
        out_alpha_gibbs = zeros(TOT_RUN,length(min));  
        out_word_gibbs = zeros(TOT_RUN,length(min));
        out_time_gibbs = zeros(TOT_RUN,length(min));
        out_rel_dif_alha = zeros(TOT_RUN,length(min));
        for k = min%k =1
            count = count + 1;
            for ii = 1:TOT_RUN%ii =1 

		        A = load(['./data_gene_0407/observed_word_matrix__K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
		        A = A.matrix;

		        WS = zeros(1,size(A,1));
		        DS = zeros(1,size(A,1));
		        for p = 1:size(A,1)
		            WS(p) = find(A(p,:));
		        end

		        sum_tmp = 0;
		        t = 0;
		        doc = 0;
		        for t = N
		            doc = doc + 1;
		            doc_index = 0;
		            while doc_index < t
		                doc_index = doc_index + 1;
		                if doc == 1
		                    DS(doc_index) = doc;
		                else
		                    DS(sum_tmp + doc_index) = doc;
		                end
		            end
		            sum_tmp = sum_tmp + N(doc);    
		        end

		        tic
		        T=TOP;
		        alpha=round(1/T, 4);
		        beta = round(1/VOC, 4);

		        BURNIN=100;
		        ITER=200; LAG=1;
		        %BURNIN=1;
		        %ITER=2; LAG=1;
		        %BURNIN=150;
		        %ITER=500; LAG=20;
		        %BURNIN=5;
		        %ITER=50; LAG=2;
		        [Phi,Theta,LL,LLall]=LDA_Gibbs_nozomi_modified(WS',DS',T,alpha,beta,ITER,BURNIN,LAG,1,VOC);
		        elapse_time = toc
		        display('Done!')
		        
		        %% Condition number
		        trd_geo_Gibbs = trd_geomcond({Phi,Phi,Phi});

		        %% Comparison to Gibbs
		        phi_matrix_original = load(['./data_gene_0407/parameter_original_K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'index_expe' num2str(ii) '.mat']);
		        generated_para = phi_matrix_original.phi_matrix_original;
		        recovered_para_gibbs = Phi;
		        PW_gibbs = [generated_para, recovered_para_gibbs];
		        rel_dif_word_gibbs = sqrt(3)*norm( generated_para - recovered_para_gibbs * (recovered_para_gibbs \ generated_para), 'fro');
		        

		        %% save data
		        clear A;
		        clear generated_para;
		        clear phi_matrix_original;
		        clear recovered_para;
		        clear DS WS
		        save(['./Gibbs0409_BI100_ITE200/matlab_output_Gibbs/VOC' int2str(VOC) '/run' int2str(ii) 'matlab_result_DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) 'BURNIN' num2str(BURNIN) 'ITER' num2str(ITER) 'LAG' num2str(LAG) '.mat'], 'Phi','Theta','LL','LLall', 'elapse_time', 'rel_dif_word_gibbs','BURNIN', 'ITER', 'LAG', 'alpha', 'trd_geo_Gibbs', 'beta');
		        save(['./Gibbs0409_BI100_ITE200/workspace_Gibbs_BI100_ITE200/VOC' int2str(VOC) '/run' int2str(ii) 'matlab_result_DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat'], '-regexp', '^(?!(out_gcond_gibbs|out_word_gibbs)$).');
		        

		        out_gcond_gibbs(ii,count) = trd_geo_Gibbs;
		        out_word_gibbs(ii,count) = rel_dif_word_gibbs;
		        out_time_gibbs(ii,count) = elapse_time;
		    end
	    end
	    save(['./Gibbs0409_BI100_ITE200/error_per_VOC_TOP/Error_summary_VOC' int2str(VOC) 'TOP' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho_minus' num2str(k) '.mat'], 'out_gcond_gibbs' ,'out_word_gibbs' ,'out_time_gibbs');
    end
end
%save('./Gibbs0409_BI100_ITE200/after2_5runs.mat');
