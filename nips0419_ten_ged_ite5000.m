cd /export/home1/NoCsBack/thesisnl/r0649141/
addpath ./topictoolbox/
addpath ./Thesis/PauseMATLAB/
addpath('./TENSORLAB')
addpath ./topictoolbox_tendec/
addpath /export/home1/NoCsBack/thesisnl/r0649141/
%load('bagofwords_nips.mat')
%load 'words_nips'

%% MAKE DIR TO WHICH YOU SAVE THE RESULTS
%mkdir nips0417_random_initialization
mkdir nips0419_ged_ite5000_tolx30 
%% CD TO THAT DIRECTOLY
cd nips0419_ged_ite5000_tolx30
%cd nips0417_random_initialization/

%% COMPILE CPP FILE USING MEX
mex /export/home1/NoCsBack/thesisnl/r0649141/topictoolbox/GibbsSamplerLDA.cpp
%A = sparse(size(DS,2) ,size(WO,1));
%ND = size(unique(DS),2);
%NDIND = zeros(1,ND);
load('bagofwords_nips.mat')
load 'words_nips'
VOC = zeros(1,size(WO,1));
if (isequal(DS,sort(DS)) ~= 1)
    error('The input arrary of Document indies is incorrect')
end
for i=1:size(DS,2)
    VOC(1,WS(i)) = VOC(1,WS(i)) + 1;
    %count = count + 1;
    %A(i, WS(i)) = 1;
    %NDIND(1,DS(i)) = NDIND(1,DS(i)) + 1;
end
vloop = [10,20,30,40,50,60,70,80,90,100];
for vs = vloop%vs=100
    [ temp1 , temp2 ] = sort(-VOC);
    WS100 = WS(ismember(WS,temp2(1:vs)));
    DS100 = DS(ismember(WS,temp2(1:vs)));
    WO100 = WO(temp2(1:vs));

    %% create matrix data for TenDec approach
    A = zeros(size(DS100,2),size(WO100,1));
    ND100 = size(unique(DS100),2);
    NDIND100 = zeros(1,ND100);
    for i=1:size(DS100,2)
        %VOC(1,WS(i)) = VOC(1,WS(i)) + 1;
        A(i, WS100(i)) = 1;
        NDIND100(1,DS100(i)) = NDIND100(1,DS100(i)) + 1;
    end    

    %% Tensor decomposition
    TOP = 2; %TOP=5
    pdir = repelem(2/TOP, TOP);%pdir = repelem(5/TOP, TOP);
    tic;
    [x,y,z,v,w] = para_recover_0131_nopause_maxite5000_tolx30(A, TOP, sum(pdir), NDIND100);
    elapse_time = toc;
    fig = figure('visible','off');
    semilogy(v.relfval);
    legend(['#TOP=' num2str(TOP) ', #VOC=' num2str(length(WO100))]);
    %title(['K' num2str(TOP) 'DN' num2str(i) 'NT' num2str(j) 'VOC' num2str(VOC) 'ortho\_minus' num2str(k) 'index\_expe' num2str(ii)'])
    xlabel('Iterations');
    ylabel('Log(refval)');
    saveas(fig,['./convergence_plot_TOP' num2str(TOP) 'VOC' num2str(length(WO100))],'epsc');
    %[S1] = WriteTopics_TenDec( x, WO100 , 7 , 0.7 );

    %Set the number of topics

    T=2;%T = 5
    %Set the hyperparameters

    BETA=0.01;
    ALPHA=2/T;%ALPHA = 5/T
    %The number of iterations

    N = 100;
    %The random seed

    SEED = 3;
    %What output to show (0=no output; 1=iterations; 2=all output)

    OUTPUT = 1;
    %This function might need a few minutes to finish

    tic
    [ WP,DP,Z ] = GibbsSamplerLDA( WS100 , DS100 , TOP , N , ALPHA , BETA , SEED , OUTPUT );
    time_gibbs = toc

    diary(['nips_gibbs_vs_tendec_VOC' num2str(vs) 'TOP' num2str(TOP) '.txt']);    
    diary on ;
    [Gibbs_sampling] = WriteTopics( WP , BETA , WO100 , 7 , 0.7 )
    [Tensor_decomposition] = WriteTopics_TenDec( x, WO100 , 7 , 0.7 )
    diary off ;
    save(['nips_gibbs_vs_tendec_VOC' num2str(vs) 'TOP' num2str(TOP) '.mat'],'x', 'y', 'z', 'v', 'w', 'WP', 'DP', 'Z', 'elapse_time', 'time_gibbs')
end
