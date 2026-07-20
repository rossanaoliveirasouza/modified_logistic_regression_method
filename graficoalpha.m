% Carregar dados
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

% Unir os dados de sem sepse e com sepse. Sabe-se que os primeiros 36 pacientes tem sepse. E que os 34 últimos não tem.
A = [Pss,Pcs]; % A possui 70 linhas e 8519 colunas

% Faz a transposta de A 
A = A';

num_pacientes_sem_sepse = size(Pss, 2); % Pega o número de colunas
num_pacientes_com_sepse = size(Pcs, 2); % Pega o número de colunas

b(1:num_pacientes_sem_sepse) = log(0.00001 / (1 - 0.00001));
b(num_pacientes_sem_sepse+1:end) = log(0.999999 / (1 - 0.99999));
[alpha,x] = resolve(A,b);

[valores, pos] = sort(alpha);

dados = table();
dados.gene = genes;
dados.alpha = alpha;
writetable(dados, 'valor_logistica_sem_treino.xlsx');

% Gráfico das posições originais versus os valores ordenados para
% visualizar alpha

figure;
plot(pos, valores, '.'); % Plota os pontos normais
hold on; % Mantém o gráfico para adicionar os destaques

% Destaca os 10 últimos pontos com estrelas vermelhas
plot(pos(end-9:end), valores(end-9:end), 'b*', 'MarkerSize', 10, 'LineWidth', 2);

% Destaca os 10 últimos pontos com estrelas vermelhas
plot(pos(1:10), valores(1:10), 'm*', 'MarkerSize', 10, 'LineWidth', 2);

% Adiciona título e rótulos
title('(GSE12624) Ordered Alpha Values for Modified Logistic Regression', 'FontSize', 16); %base 1
%title('(GSM333448) Ordered Alpha Values for Modified Logistic Regression', 'FontSize', 16); % base 2
xlabel('Gene ID', 'FontSize', 14);
ylabel('Alpha Values', 'FontSize', 14);
grid on;



