function  [x,y,z,v,w]  = para_recover_0131_nopause_maxite300(matrix, rank, alpha0, num_term_vector)
A = matrix;
%{
ichi=input('Is your input {#words in all documents(row) * #distict words in corpus(column)) matrix of unit-vector ?(yes/no)\n' ,'s');
if strcmpi(ichi, 'yes')
A = matrix;
elseif strcmpi(ichi, 'no')
Ni = input('Is your input {#distict words in corpus(row) * #words in all documents(column)) matrix of unit-vector ?(yes/no)\n' ,'s');
    if strcmpi(Ni, 'yes')
        A = matrix.';
    elseif strcmpi(Ni, 'no')
        error('Please make an appropriate input of word matrix')
    else
        error('Please type "yes" or "no"');
    end
else
error('Please type "yes" or "no"');
end
%}

%A = matrix;
R = rank; % # (unobserable) topics
a = alpha0; % sum of all elements of the alpha prior vector for Dirichlet dis.
N = num_term_vector;
ndoc = length(num_term_vector);
%E(x)
E_M1 = mean(A,1);
%alpha = 1

%number of terms in each doc--list
% e.g. 100 doc and 100 terms---> ex; N = 20
%N = [100]; 
%N = repelem(N,100);


%N = [20];
%N = repelem(20,50);

% #documents --> ndoc = 100;
% # vocabulary --> nvoc = 1000;
%ndoc = 50;
%create #documet x #words in dictionary matrix by adding up all counts in a
%doc for the same word.


% looping for each doc with #terms each docment have
start = 1;
stop = 0;
nvoc = size(A,2);
%ndoc = size()
freq = nan(ndoc, nvoc);
count_doc = 0;
m2_den = 0;
m3_den = 0;
for i = N
    count_doc = count_doc + 1;
    stop = stop + i;
    freq(count_doc, :) = sum(A(start:stop,:));
    start = stop + 1;
    %disp(start);
    %disp(stop);
    m2_den = m2_den + i*(i-1);
    m3_den = m3_den + i*(i-1)*(i-2); % denominator for M3 ()
end
disp('If any NaN is outputted in following, error is occuring and script is wrong. Push any key')
%pause
freq(isnan(freq)==1) % check if error (NaN still exists) occurs
disp('Push any key')
%pause;


test = nchoosek(1:size(freq,2),2);
m2_til = nan(size(freq,2));
for i= 1:size(test,1)
    %FOR PAIR-WORDS
    test2 = 0;
    for j = 1:size(freq,1)
        one = freq(j,test(i,1));
        two = freq(j,test(i,2));
        test2 = test2 + one*two;
    end
    m2_til(test(i,1),test(i,2)) = test2;
    m2_til(test(i,2),test(i,1)) = test2;
end
% calculate diagonal element
dia_ele = size(zeros(1,size(freq,2)));
for i = 1:size(freq,2)
    test2 = 0;
    for j = 1:size(freq,1)
        one = freq(j,i);
        two = freq(j,i) - 1;
        test2 = test2 + one*two;
    end
    dia_ele(i) = test2;
end
%diago = diag(dia_ele);
% replace diagonal element 
m2_til(logical(eye(size(m2_til)))) = dia_ele;
% by dividing denominator term, it is done.
disp('If any NaN is outputted in following, error is occuring and script is wrong. Push any key')
%pause
m2_til(isnan(m2_til)==1) % check if error (NaN still exists) occurs
disp('Push any key')
%pause;

E_M2 = m2_til/m2_den;

%test3 = nan([1000,1000,1000])
test3 = nan([size(freq,2),size(freq,2),size(freq,2)]);
TRIPLE = nchoosek(1:size(freq,2),3);
%E3_count = zeros(size(A,2),size(A,2),size(A,2));
%Calculate frequency cubic for (X1,X2,X3)
%non symmetric element
%test3 = zeros(size(freq,2),size(freq,2),size(freq,2));
%M_12 = nan([size(freq,2),size(freq,2),size(freq,2)]); % create suppliment tensor for definition of M3 (pp.8)
for i= 1:size(TRIPLE,1)
    %FOR TRIPLER
    test2 = 0;
%    test1 = 0;
    for j = 1:size(freq,1)
        one = freq(j,TRIPLE(i,1));
        two = freq(j,TRIPLE(i,2));
        three = freq(j,TRIPLE(i,3));
        test2 = test2 + one*two*three;
%        test1 = E_M2(TRIPLE(i,1),TRIPLE(i,2))*E_M1(TRIPLE(i,3)) + E_M2(TRIPLE(i,2),TRIPLE(i,3))*E_M1(TRIPLE(i,1)) + E_M2(TRIPLE(i,3), TRIPLE(i,1))*E_M1(TRIPLE(i,2));
    end
    test3(TRIPLE(i,1),TRIPLE(i,2),TRIPLE(i,3)) = test2 ; test3(TRIPLE(i,1),TRIPLE(i,3),TRIPLE(i,2)) = test2;
    test3(TRIPLE(i,2),TRIPLE(i,1),TRIPLE(i,3)) = test2 ; test3(TRIPLE(i,2),TRIPLE(i,3),TRIPLE(i,1)) = test2;
    test3(TRIPLE(i,3),TRIPLE(i,1),TRIPLE(i,2)) = test2 ; test3(TRIPLE(i,3),TRIPLE(i,2),TRIPLE(i,1)) = test2;
%    M_12(TRIPLE(i,1),TRIPLE(i,2),TRIPLE(i,3)) = test1;
end

%diagonal element
%dia_ele = size(zeros(1,size(freq,2)));
for i = 1:size(freq,2)
    test2 = 0;
    for j = 1:size(freq,1)
        one = freq(j,i);
        two = freq(j,i) - 1;
        three = freq(j,i) - 2;
        test2 = test2 + one*two*three;
    end
    test3(i,i,i) = test2;
end

%semi-symmetics element
        

for i= 1:size(test,1)
    %FOR PAIR-WORDS
    test2 = 0;
    test4 = 0;
    % consider two cases: 'one occurs twice' and 'two occurs twice'
    for j = 1:size(freq,1)
        one = freq(j,test(i,1));
        two = freq(j,test(i,2));
        test2 = test2 + one*two*(two-1);
        test4 = test4 + one*(one - 1)*two;
    end
    test3(test(i,1),test(i,1),test(i,2)) = test4; test3(test(i,1),test(i,2),test(i,1)) = test4; test3(test(i,2),test(i,1),test(i,1)) = test4;
    test3(test(i,2),test(i,2),test(i,1)) = test2; test3(test(i,2),test(i,1),test(i,2)) = test2; test3(test(i,1),test(i,2),test(i,2)) = test2;
end
disp('If any  is outputted, error is occuring and script is wrong')
test3(isnan(test3)==1) % check if error (NaN still exists) occurs
%pause;

% obtained tensor is divided by denominator and done.
E_M3 = test3/m3_den;

%M1 o M1 (cross moment for 1st morment)
outerproduct = outprod(E_M1,E_M1);
%Definition of 2nd moment
M2 = E_M2 - (a/(a+1))*outerproduct; % a = 1 (sum of all elements of alpha_vector = 1)
%Definition of 3rd moment
    %Definition of M12 moment (tensor)
M_12 = nan([size(freq,2),size(freq,2),size(freq,2)]); % create suppliment tensor for definition of M3 (pp.8)
for h = 1:size(freq,2)
    for l = 1:size(freq,2)
        for m = 1:size(freq,2)           
            M_12(h,l,m) = E_M2(h,l)*E_M1(m) + E_M2(l,m)*E_M1(h) + E_M2(m,h)*E_M1(l);
            %disp([h,l,m])
            %pause
        end
    end
end
disp('If any  is outputted, error is occuring and script is wrong')
M_12(isnan(M_12)==1)% check if error (NaN still exists) occurs
%pause;
    %Definition of 3rd moment--final--
M3 = E_M3 - (a/(a+2))*M_12 + (2*a^(2)/(a+2)/(a+1))*outprod(E_M1,E_M1,E_M1);

%%                      %%
%% Recover parameter mu %%
%%                      %%

T = M3;
%[tr1, tr2] = trd_geomcond(cpd(T, R));
% Create model
model = struct;
model.variables.a = randn(size(T,1), R);
model.variables.extra = randn(1, R);
model.factors.A = {'a',@struct_nonneg};
model.factors.C = {'extra',@struct_nonneg};
model.factorizations.symm.data = T;
model.factorizations.symm.cpd  = {'A','A','A','C'};
sdf_check(model, 'print')
sdf_check(model);
disp('Please press a key if no errors are made !')% Press a key here.You can see the message 'Paused: Press any key' in        % the lower left corner of MATLAB window.
%pause;
%options.Display = 5; % View convergence progress every 5 iterations.
%% option statement
options.MaxIter = 5000;
%options.TolFun  = 1e-20 ;
%options.TolX    = 1e-20 ;
%options.TolFun  = 1e-22 ; around 3000 ites
%options.TolX    = 1e-22 ;
options.TolFun  = 1e-30 ;
options.TolX    = 1e-30 ;
%
[sol,out] = sdf_nls(model,options);
%[tr1, tr2] = trd_geomcond(sdf_nls(model,options));
solvar = sol.variables;
%factor matrices
solfac = sol.factors;
%unnormalized mu
topic_word_pro = zeros(size(A,2),R);
alpha      = zeros(1,R);
for i=1:R
    alpha(1,i) = solfac.C(1,i)*norm(solfac.A(:,i),1)*norm(solfac.A(:,i),1)*norm(solfac.A(:,i),1)*(a*(a+1)*(a+2))/2;
    topic_word_pro(:,i) = solfac.A(:,i)/ norm(solfac.A(:,i),1);
end

x = topic_word_pro;
y = alpha;
z = sum(alpha);
v = out;
w = sol;
%r = solvar;
fprintf('Z should be 1. Please check... \n ');
fprintf('Type x to see word * topic matrix. \n Type y to see vectors of alpha_{k}. \n Type z to see alpha_{0} \n ');
end