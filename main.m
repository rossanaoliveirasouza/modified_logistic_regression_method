%% ============================
%  BASE 1
%  ============================

load("base1/sem_sepse.mat");
load("base1/com_sepse.mat");

% Converter tabelas em arrays
pss = table2array(pacientes_sem_sepse);
pcs = table2array(pacientes_com_sepse);
genes = pacientes_com_sepse.(1);
genes = string(genes);

% Obter tamanhos
[linhas_s, colunas_s] = size(pss);
[linhas_c, colunas_c] = size(pcs);

% Inicializar arrays para armazenar os dados convertidos
Pss = zeros(linhas_s, colunas_s);
Pcs = zeros(linhas_c, colunas_c);

% Converter pss de char para numérico
for i = 1:linhas_s
    for j = 1:colunas_s
        aux = pss{i,j};
        aux = char(aux);
        Pss(i,j) = str2double(aux); % Usar str2double para tratar NaN corretamente
    end
end

% Converter pcs de char para numérico
for i = 1:linhas_c
    for j = 1:colunas_c
        aux = pcs{i,j};
        aux = char(aux);
        Pcs(i,j) = str2double(aux);
    end
end

% Remover a coluna de genes dos arrays, já está mapeado na variável genes
Pss = [Pss(:,2:end)];
Pcs = [Pcs(:,2:end)];

% Unir os dados de sem sepse e com sepse. Sabe-se que os 36 primeiros não tem sepse. E que os últimos 34 pacientes tem sepse.
A = [Pss,Pcs]; 

% Faz a transposta de A 
A = A'; % A possui 70 linhas e 8519 colunas

num_pacientes_sem_sepse = size(Pss, 2);
num_pacientes_com_sepse = size(Pcs, 2); 

classificador = [zeros(1, num_pacientes_sem_sepse), ones(1, num_pacientes_com_sepse)];

% Número de folds
K = 2;

fprintf('\n====================================\n');
fprintf('RESULTADOS DA BASE 1\n');
fprintf('====================================\n');

fprintf('\n====================================\n');
fprintf('MEU MÉTODO\n');
fprintf('====================================\n');
crossvalidation(A, classificador, genes, K);

fprintf('\n====================================\n');
fprintf('ADHOC\n');
fprintf('====================================\n');

k0 = 12;
k1 = 12;
%crossvalidation_adhoc(A, classificador, genes, K, k0, k1);

fprintf('\n====================================\n');
fprintf('LASSO LOGISTIC REGRESSION\n');
fprintf('====================================\n');
%crossvalidation_logistic_baselines(A, classificador, genes, K, "lasso", 20);

fprintf('\n====================================\n');
fprintf('ELASTIC NET LOGISTIC REGRESSION\n');
fprintf('====================================\n');
%crossvalidation_logistic_baselines(A, classificador, genes, K, "elasticnet", 20);

%% ============================
%  BASE 2
%  ============================

clear pacientes_sem_sepse pacientes_com_sepse pss pcs Pss Pcs A classificador genes

load("base2/sem_sepse.mat");
load("base2/com_sepse.mat");

pss = table2array(pacientes_sem_sepse_base2);
pcs = table2array(pacientes_com_sepse_base2);
genes = pacientes_com_sepse_base2.(1);
genes = string(genes);

[linhas_s, colunas_s] = size(pss);
[linhas_c, colunas_c] = size(pcs);

Pss = zeros(linhas_s, colunas_s);
Pcs = zeros(linhas_c, colunas_c);

for i = 1:linhas_s
    for j = 1:colunas_s
        aux = pss{i,j};
        aux = char(aux);
        Pss(i,j) = str2double(aux);
    end
end

for i = 1:linhas_c
    for j = 1:colunas_c
        aux = pcs{i,j};
        aux = char(aux);
        Pcs(i,j) = str2double(aux);
    end
end

Pss = Pss(:,2:end);
Pcs = Pcs(:,2:end);

A = [Pss, Pcs];
A = A';

num_pacientes_sem_sepse = size(Pss, 2);
num_pacientes_com_sepse = size(Pcs, 2);

classificador = [zeros(1, num_pacientes_sem_sepse), ones(1, num_pacientes_com_sepse)];

% Número de folds
K = 2;

fprintf('\n====================================\n');
fprintf('RESULTADOS DA BASE 2\n');
fprintf('====================================\n');

fprintf('\n====================================\n');
fprintf('MEU MÉTODO\n');
fprintf('====================================\n');
crossvalidation(A, classificador, genes, K);

fprintf('\n====================================\n');
fprintf('ADHOC\n');
fprintf('====================================\n');

k0 = 8;
k1 = 8;
%crossvalidation_adhoc(A, classificador, genes, K, k0, k1);


fprintf('\n====================================\n');
fprintf('LASSO LOGISTIC REGRESSION\n');
fprintf('====================================\n');
%crossvalidation_logistic_baselines(A, classificador, genes, K, "lasso", 20);

fprintf('\n====================================\n');
fprintf('ELASTIC NET LOGISTIC REGRESSION\n');
fprintf('====================================\n');
%crossvalidation_logistic_baselines(A, classificador, genes, K, "elasticnet", 20);

%% ============================
%  BASE 3
%  ============================

clear pacientes_sem_sepse pacientes_com_sepse pss pcs Pss Pcs A classificador genes

load("base3/sem_sepse.mat");
load("base3/com_sepse.mat");

pss = table2array(pacientes_sem_sepse_base3);
pcs = table2array(pacientes_com_sepse_base3);
genes = pacientes_com_sepse_base3.(1);
genes = string(genes);

[linhas_s, colunas_s] = size(pss);
[linhas_c, colunas_c] = size(pcs);

Pss = zeros(linhas_s, colunas_s);
Pcs = zeros(linhas_c, colunas_c);

for i = 1:linhas_s
    for j = 1:colunas_s
        aux = pss{i,j};
        aux = char(aux);
        Pss(i,j) = str2double(aux);
    end
end

for i = 1:linhas_c
    for j = 1:colunas_c
        aux = pcs{i,j};
        aux = char(aux);
        Pcs(i,j) = str2double(aux);
    end
end

Pss = Pss(:,2:end);
Pcs = Pcs(:,2:end);

A = [Pss, Pcs];
A = A';

num_pacientes_sem_sepse = size(Pss, 2);
num_pacientes_com_sepse = size(Pcs, 2);

classificador = [zeros(1, num_pacientes_sem_sepse), ones(1, num_pacientes_com_sepse)];

% Número de folds
K = 2;

fprintf('\n====================================\n');
fprintf('RESULTADOS DA BASE 3\n');
fprintf('====================================\n');

fprintf('\n====================================\n');
fprintf('MEU MÉTODO\n');
fprintf('====================================\n');
crossvalidation(A, classificador, genes, K);

fprintf('\n====================================\n');
fprintf('ADHOC\n');
fprintf('====================================\n');

k0 = 30;
k1 = 30;
%crossvalidation_adhoc(A, classificador, genes, K, k0, k1);


fprintf('\n====================================\n');
fprintf('LASSO LOGISTIC REGRESSION\n');
fprintf('====================================\n');
%crossvalidation_logistic_baselines(A, classificador, genes, K, "lasso", 20);

fprintf('\n====================================\n');
fprintf('ELASTIC NET LOGISTIC REGRESSION\n');
fprintf('====================================\n');
%crossvalidation_logistic_baselines(A, classificador, genes, K, "elasticnet", 20);