function [] = crossvalidation(A, classificador, genes, K)

% Iniciar contagem do tempo total da validação cruzada
tempo_total_inicio = tic;

%% Dividir os pacientes em K grupos estratificados. 

% Garantir que o classificador seja vetor coluna
classificador = classificador(:);

% Número de genes mais importantes selecionados em cada fold
num_genes_selecionados = 20;

% Valor pequeno para transformar classes em log-odds
eps_prob = 1e-5;

% Criar partição estratificada
cv = cvpartition(classificador, 'KFold', K);

% Inicializar vetores para guardar os resultados de todos os pacientes
scores_todos = zeros(size(A,1), 1);
probabilidades_todos = zeros(size(A,1), 1);
predicoes_todos = zeros(size(A,1), 1);

% Inicializar células para guardar informações por fold
genes_selecionados_por_fold = cell(K, 1);
indices_genes_por_fold = cell(K, 1);
alphas_por_fold = cell(K, 1);

% Inicializar vetores para armazenar tempos computacionais por fold
tempo_fold = zeros(K, 1);
tempo_alpha = zeros(K, 1);
tempo_selecao = zeros(K, 1);
tempo_refit = zeros(K, 1);
tempo_predicao = zeros(K, 1);


%% Iniciar a cross validação 

for k = 1:K

    % Iniciar contagem do tempo total deste fold
    tempo_fold_inicio = tic;

   % fprintf('\n=============================\n');
   % fprintf('Fold %d de %d\n', k, K);
   % fprintf('=============================\n');

    idx_treino = training(cv, k);
    idx_teste  = test(cv, k);

    %% Separar dados
    A_treino = A(idx_treino, :);
    A_teste  = A(idx_teste, :);

    y_treino = classificador(idx_treino);
    y_teste  = classificador(idx_teste);

    % Inicializar vetor b_treino
    b_treino = zeros(length(y_treino), 1); 
    b_treino(y_treino == 0) = log(eps_prob / (1 - eps_prob)); % Pacientes sem sepse recebem log-odds próximo de probabilidade 0
    b_treino(y_treino == 1) = log((1 - eps_prob) / eps_prob); % Pacientes com sepse recebem log-odds próximo de probabilidade 1

    % Iniciar contagem do tempo para calcular alpha
    tempo_alpha_inicio = tic;

    % Calcular alpha usando apenas A_treino e b_treino
    [alpha, x] = resolve(A_treino, b_treino);

    % Armazenar tempo para calcular alpha
    tempo_alpha(k) = toc(tempo_alpha_inicio);

    % Iniciar contagem do tempo de seleção de genes
    tempo_selecao_inicio = tic;

    %% Ranquear os genes por |alpha_i|. 
    
    % %Importância de cada gene é o módulo do coeficiente alpha
    importancia_genes = abs(alpha);

    % Ordenar genes do mais importante para o menos importante
    [~, ordem_genes] = sort(importancia_genes, 'descend');

    % Selecionar os índices dos genes mais importantes
    indices_genes_selecionados = ordem_genes(1:num_genes_selecionados);

    % Obter os nomes dos genes selecionados
    genes_selecionados = genes(indices_genes_selecionados);

    % Armazenar tempo de seleção de genes
    tempo_selecao(k) = toc(tempo_selecao_inicio);

    %% Aplicar esses genes nos pacientes de teste

    % Manter apenas os genes selecionados no treino
    A_treino_sel = A_treino(:, indices_genes_selecionados);

    % Manter os mesmos genes selecionados no teste
    A_teste_sel = A_teste(:, indices_genes_selecionados);

    % Iniciar contagem do tempo de refit com genes selecionados
    tempo_refit_inicio = tic;

    % Recalcular alpha usando apenas os genes selecionados no treino
    [alpha_sel, ~] = resolve(A_treino_sel, b_treino);

    % Armazenar tempo de refit
    tempo_refit(k) = toc(tempo_refit_inicio);

    % Iniciar contagem do tempo de predição
    tempo_predicao_inicio = tic;

    %% Calcular score e probabilidade dos pacientes de teste

    % Score linear no conjunto de teste
    score_test = A_teste_sel * alpha_sel;

    % Converter score em probabilidade pela sigmoide
    prob_test = 1 ./ (1 + exp(-score_test));

    % Predição final usando limiar 0.5
    pred_test = prob_test >= 0.5;

    % Armazenar tempo de predição
    tempo_predicao(k) = toc(tempo_predicao_inicio);

    % Guardar os resultados deste fold

    % Guardar scores nas posições originais dos pacientes de teste
    scores_todos(idx_teste) = score_test;

    % Guardar probabilidades nas posições originais dos pacientes de teste
    probabilidades_todos(idx_teste) = prob_test;

    % Guardar predições nas posições originais dos pacientes de teste
    predicoes_todos(idx_teste) = pred_test;

    % Guardar genes selecionados neste fold
    genes_selecionados_por_fold{k} = genes_selecionados;

    % Guardar índices dos genes selecionados neste fold
    indices_genes_por_fold{k} = indices_genes_selecionados;

    % Guardar alpha final deste fold
    alphas_por_fold{k} = alpha_sel;

    %% Mostrar resumo do fold

    fprintf('Treino: %d pacientes sem sepse, %d com sepse\n', ...
        sum(y_treino == 0), sum(y_treino == 1));

    fprintf('Teste:  %d pacientes sem sepse, %d com sepse\n', ...
        sum(y_teste == 0), sum(y_teste == 1));

    %fprintf('Primeiros 10 genes selecionados neste fold:\n');
    %disp(genes_selecionados(1:min(10, num_genes_selecionados)));

    % Armazenar tempo total deste fold
    tempo_fold(k) = toc(tempo_fold_inicio);

end

% Armazenar tempo total da validação cruzada
tempo_total_cv = toc(tempo_total_inicio);

%% ============================================================
%  Avaliação final juntando as predições de todos os folds
%  ============================================================

% Garantir que os vetores estejam como coluna
y_real = classificador(:);
y_pred = predicoes_todos(:);
y_prob = probabilidades_todos(:);

% Calcular matriz de confusão
% Linha = classe real
% Coluna = classe prevista
matriz_confusao = confusionmat(y_real, y_pred);

% Extrair valores da matriz de confusão
% Considerando:
% classe 0 = sem sepse
% classe 1 = com sepse
TN = matriz_confusao(1,1);
FP = matriz_confusao(1,2);
FN = matriz_confusao(2,1);
TP = matriz_confusao(2,2);

% Calcular acurácia
acuracia = (TP + TN) / (TP + TN + FP + FN);

% Calcular sensibilidade
% Também chamada de recall ou taxa de verdadeiros positivos
sensibilidade = TP / (TP + FN);

% Calcular especificidade
% Taxa de verdadeiros negativos
especificidade = TN / (TN + FP);

% Calcular precisão
precisao = TP / (TP + FP);

% Calcular F1-score
F1 = 2 * (precisao * sensibilidade) / (precisao + sensibilidade);

% Calcular curva ROC e AUC
% A classe positiva é 1, ou seja, com sepse
[Xroc, Yroc, limiares, AUC] = perfcurve(y_real, y_prob, 1);

% Exibir resultados
fprintf('Acurácia:        %.4f\n', acuracia);
fprintf('Sensibilidade:   %.4f\n', sensibilidade);
fprintf('Especificidade:  %.4f\n', especificidade);
fprintf('Precisão:        %.4f\n', precisao);
fprintf('F1-score:        %.4f\n', F1);
fprintf('AUC:             %.4f\n', AUC);

fprintf('\nMatriz de confusão:\n');
disp(matriz_confusao);

%% ============================================================
%  Tempos computacionais
%  ============================================================

fprintf('\n============================================\n');
fprintf('TEMPOS COMPUTACIONAIS\n');
fprintf('============================================\n');

fprintf('Tempo total CV:              %.6f segundos\n', tempo_total_cv);
fprintf('Tempo médio por fold:        %.6f segundos\n', mean(tempo_fold));
fprintf('Desvio-padrão por fold:      %.6f segundos\n', std(tempo_fold));

fprintf('\nTempo médio alpha inicial:    %.6f segundos\n', mean(tempo_alpha));
fprintf('Tempo médio seleção:         %.6f segundos\n', mean(tempo_selecao));
fprintf('Tempo médio refit:           %.6f segundos\n', mean(tempo_refit));
fprintf('Tempo médio predição:        %.6f segundos\n', mean(tempo_predicao));

fprintf('\nTempo total alpha inicial:    %.6f segundos\n', sum(tempo_alpha));
fprintf('Tempo total seleção:         %.6f segundos\n', sum(tempo_selecao));
fprintf('Tempo total refit:           %.6f segundos\n', sum(tempo_refit));
fprintf('Tempo total predição:        %.6f segundos\n', sum(tempo_predicao));

%% Plotar curva ROC

%figure;
%plot(Xroc, Yroc, 'LineWidth', 2);
%xlabel('Taxa de Falsos Positivos');
%ylabel('Taxa de Verdadeiros Positivos');
%title(['Curva ROC - AUC = ', num2str(AUC, '%.4f')]);
%grid on;

%% ============================================================
%  Ver genes selecionados com maior frequência nos folds
%  ============================================================

% Juntar todos os índices de genes selecionados em todos os folds
todos_indices_genes = [];

for k = 1:K
    todos_indices_genes = [todos_indices_genes; indices_genes_por_fold{k}(:)];
end

% Encontrar genes únicos selecionados
genes_unicos = unique(todos_indices_genes);

% Contar quantas vezes cada gene apareceu
frequencia_genes = zeros(length(genes_unicos), 1);

for i = 1:length(genes_unicos)
    frequencia_genes(i) = sum(todos_indices_genes == genes_unicos(i));
end

% Ordenar genes por frequência, do maior para o menor
%[frequencia_ordenada, ordem_freq] = sort(frequencia_genes, 'descend');

% Obter os índices dos genes mais frequentes
%indices_genes_frequentes = genes_unicos(ordem_freq);

% Obter os nomes dos genes mais frequentes
%nomes_genes_frequentes = genes(indices_genes_frequentes);

% Criar tabela com os genes mais frequentes
%tabela_genes_frequentes = table( ...
%    nomes_genes_frequentes, ...
%    indices_genes_frequentes, ...
%    frequencia_ordenada, ...
%    'VariableNames', {'Gene', 'Indice', 'Frequencia'} ...
%);

% Exibir tabela completa
%fprintf('\n============================================\n');
%fprintf('GENES SELECIONADOS COM MAIOR FREQUÊNCIA\n');
%fprintf('============================================\n');

%disp(tabela_genes_frequentes);

% Exibir apenas os 20 genes mais frequentes
%topN = min(20, height(tabela_genes_frequentes));

%fprintf('\nTop %d genes mais frequentes:\n', topN);
%disp(tabela_genes_frequentes(1:topN, :));

%% Plotar frequência dos genes mais selecionados

%figure;
%bar(frequencia_ordenada(1:topN));
%xticks(1:topN);
%xticklabels(nomes_genes_frequentes(1:topN));
%xtickangle(45);
%ylabel('Frequência nos folds');
%title('Genes selecionados com maior frequência');
%grid on;