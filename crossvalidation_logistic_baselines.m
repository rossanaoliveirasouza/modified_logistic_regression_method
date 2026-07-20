function [] = crossvalidation_logistic_baselines(A, classificador, genes, K, metodo, num_features)

% A: matriz pacientes x genes
% classificador: vetor 0/1, onde 0 = sem sepse e 1 = com sepse
% genes: vetor com os nomes dos genes
% K: número de folds
% metodo: "classic", "lasso" ou "elasticnet"
% num_features: número de genes usados na logística clássica após filtragem

tempo_total_inicio = tic;

classificador = classificador(:);
genes = genes(:);

cv = cvpartition(classificador, 'KFold', K);

probabilidades_todos = zeros(size(A,1), 1);
predicoes_todos = zeros(size(A,1), 1);
num_features_por_fold = zeros(K, 1);

% Guardar genes selecionados em cada fold
indices_genes_por_fold = cell(K, 1);
genes_selecionados_por_fold = cell(K, 1);

% Guardar tempos por fold
tempo_fold = zeros(K, 1);
tempo_preprocessamento = zeros(K, 1);
tempo_selecao = zeros(K, 1);
tempo_treino = zeros(K, 1);
tempo_predicao = zeros(K, 1);

fprintf('\nMétodo: %s\n', metodo);
fprintf('====================================\n');

for k = 1:K

    tempo_fold_inicio = tic;

    idx_treino = training(cv, k);
    idx_teste  = test(cv, k);

    A_treino = A(idx_treino, :);
    A_teste  = A(idx_teste, :);

    y_treino = classificador(idx_treino);
    y_teste  = classificador(idx_teste);

    fprintf('Treino: %d pacientes sem sepse, %d com sepse\n', ...
        sum(y_treino == 0), sum(y_treino == 1));

    fprintf('Teste:  %d pacientes sem sepse, %d com sepse\n', ...
        sum(y_teste == 0), sum(y_teste == 1));

    %% Preprocessamento: padronização usando apenas treino
    t = tic;

    media_treino = mean(A_treino, 1);
    desvio_treino = std(A_treino, 0, 1);
    desvio_treino(desvio_treino == 0) = 1;

    A_treino_z = (A_treino - media_treino) ./ desvio_treino;
    A_teste_z  = (A_teste  - media_treino) ./ desvio_treino;

    tempo_preprocessamento(k) = toc(t);

    switch metodo

        case "classic"

            %% Seleção de genes por filtragem univariada
            t = tic;

            scores_features = zeros(1, size(A_treino_z, 2));

            for j = 1:size(A_treino_z, 2)
                grupo0 = A_treino_z(y_treino == 0, j);
                grupo1 = A_treino_z(y_treino == 1, j);

                scores_features(j) = abs(mean(grupo1) - mean(grupo0));
            end

            [~, ordem] = sort(scores_features, 'descend');

            num_sel = min(num_features, size(A_treino_z, 2));
            idx_genes = ordem(1:num_sel);

            Xtr = A_treino_z(:, idx_genes);
            Xte = A_teste_z(:, idx_genes);

            num_features_por_fold(k) = num_sel;

            tempo_selecao(k) = toc(t);

            %% Treino
            t = tic;

            try
                modelo = fitglm(Xtr, y_treino, ...
                    'Distribution', 'binomial', ...
                    'Link', 'logit');

                treino_falhou = false;

            catch
                warning('Regressão logística clássica falhou neste fold. Usando probabilidade média do treino.');

                prob_media = mean(y_treino);
                treino_falhou = true;
            end

            tempo_treino(k) = toc(t);

            %% Predição
            t = tic;

            if treino_falhou
                prob_test = prob_media * ones(size(y_teste));
            else
                prob_test = predict(modelo, Xte);
            end

            tempo_predicao(k) = toc(t);

        case "lasso"

            %% Treino com seleção embutida pelo LASSO
            t = tic;

            cv_interno = min(5, sum(y_treino == 0));
            cv_interno = min(cv_interno, sum(y_treino == 1));

            if cv_interno < 2
                cv_interno = 2;
            end

            [B, FitInfo] = lassoglm(A_treino_z, y_treino, ...
                'binomial', ...
                'Alpha', 1, ...
                'CV', cv_interno);

            idxLambda = FitInfo.IndexMinDeviance;

            coef = B(:, idxLambda);
            intercepto = FitInfo.Intercept(idxLambda);

            tempo_treino(k) = toc(t);

            %% Seleção de genes
            t = tic;

            idx_genes = find(coef ~= 0);
            num_features_por_fold(k) = length(idx_genes);

            tempo_selecao(k) = toc(t);

            %% Predição
            t = tic;

            prob_test = 1 ./ (1 + exp(-(A_teste_z * coef + intercepto)));

            tempo_predicao(k) = toc(t);

        case "elasticnet"

            %% Treino com seleção embutida pelo Elastic Net
            t = tic;

            cv_interno = min(5, sum(y_treino == 0));
            cv_interno = min(cv_interno, sum(y_treino == 1));

            if cv_interno < 2
                cv_interno = 2;
            end

            [B, FitInfo] = lassoglm(A_treino_z, y_treino, ...
                'binomial', ...
                'Alpha', 0.5, ...
                'CV', cv_interno);

            idxLambda = FitInfo.IndexMinDeviance;

            coef = B(:, idxLambda);
            intercepto = FitInfo.Intercept(idxLambda);

            tempo_treino(k) = toc(t);

            %% Seleção de genes
            t = tic;

            idx_genes = find(coef ~= 0);
            num_features_por_fold(k) = length(idx_genes);

            tempo_selecao(k) = toc(t);

            %% Predição
            t = tic;

            prob_test = 1 ./ (1 + exp(-(A_teste_z * coef + intercepto)));

            tempo_predicao(k) = toc(t);

        otherwise
            error('Método inválido. Use "classic", "lasso" ou "elasticnet".');
    end

    % Guardar genes selecionados no fold
    indices_genes_por_fold{k} = idx_genes(:);
    genes_selecionados_por_fold{k} = genes(idx_genes);

    pred_test = prob_test >= 0.5;

    probabilidades_todos(idx_teste) = prob_test;
    predicoes_todos(idx_teste) = pred_test;

    tempo_fold(k) = toc(tempo_fold_inicio);

end

tempo_total_cv = toc(tempo_total_inicio);

%% Avaliação final

y_real = classificador(:);
y_pred = predicoes_todos(:);
y_prob = probabilidades_todos(:);

matriz_confusao = confusionmat(y_real, y_pred);

TN = matriz_confusao(1,1);
FP = matriz_confusao(1,2);
FN = matriz_confusao(2,1);
TP = matriz_confusao(2,2);

acuracia = (TP + TN) / (TP + TN + FP + FN);

if (TP + FN) == 0
    sensibilidade = 0;
else
    sensibilidade = TP / (TP + FN);
end

if (TN + FP) == 0
    especificidade = 0;
else
    especificidade = TN / (TN + FP);
end

if (TP + FP) == 0
    precisao = 0;
else
    precisao = TP / (TP + FP);
end

if (precisao + sensibilidade) == 0
    F1 = 0;
else
    F1 = 2 * (precisao * sensibilidade) / (precisao + sensibilidade);
end

[~, ~, ~, AUC] = perfcurve(y_real, y_prob, 1);

fprintf('Acurácia:        %.4f\n', acuracia);
fprintf('Sensibilidade:   %.4f\n', sensibilidade);
fprintf('Especificidade:  %.4f\n', especificidade);
fprintf('Precisão:        %.4f\n', precisao);
fprintf('F1-score:        %.4f\n', F1);
fprintf('AUC:             %.4f\n', AUC);
fprintf('# Features médio: %.2f\n', mean(num_features_por_fold));

fprintf('\nMatriz de confusão:\n');
disp(matriz_confusao);

%% Tempos computacionais

fprintf('\n====================================\n');
fprintf('TEMPOS COMPUTACIONAIS\n');
fprintf('====================================\n');

fprintf('Tempo total CV:              %.6f segundos\n', tempo_total_cv);
fprintf('Tempo médio por fold:        %.6f segundos\n', mean(tempo_fold));
fprintf('Desvio-padrão por fold:      %.6f segundos\n', std(tempo_fold));

fprintf('\nTempo médio preprocessamento: %.6f segundos\n', mean(tempo_preprocessamento));
fprintf('Tempo médio seleção:         %.6f segundos\n', mean(tempo_selecao));
fprintf('Tempo médio treino:          %.6f segundos\n', mean(tempo_treino));
fprintf('Tempo médio predição:        %.6f segundos\n', mean(tempo_predicao));

fprintf('\nTempo total preprocessamento: %.6f segundos\n', sum(tempo_preprocessamento));
fprintf('Tempo total seleção:         %.6f segundos\n', sum(tempo_selecao));
fprintf('Tempo total treino:          %.6f segundos\n', sum(tempo_treino));
fprintf('Tempo total predição:        %.6f segundos\n', sum(tempo_predicao));

%% Genes selecionados por fold

%fprintf('\n====================================\n');
%fprintf('GENES SELECIONADOS POR FOLD\n');
%fprintf('====================================\n');

for k = 1:K
    %fprintf('\nFold %d - %d genes selecionados:\n', k, length(indices_genes_por_fold{k}));
  %  disp(genes_selecionados_por_fold{k});
end

%% Frequência dos genes selecionados

todos_indices_genes = [];

for k = 1:K
    todos_indices_genes = [todos_indices_genes; indices_genes_por_fold{k}(:)];
end

genes_unicos = unique(todos_indices_genes);

frequencia_genes = zeros(length(genes_unicos), 1);

for i = 1:length(genes_unicos)
    frequencia_genes(i) = sum(todos_indices_genes == genes_unicos(i));
end

[frequencia_ordenada, ordem_freq] = sort(frequencia_genes, 'descend');

indices_genes_frequentes = genes_unicos(ordem_freq);
nomes_genes_frequentes = genes(indices_genes_frequentes);

tabela_genes_frequentes = table( ...
    nomes_genes_frequentes, ...
    indices_genes_frequentes, ...
    frequencia_ordenada, ...
    'VariableNames', {'Gene', 'Indice', 'Frequencia'} ...
);

fprintf('\n====================================\n');
fprintf('FREQUÊNCIA DOS GENES SELECIONADOS\n');
fprintf('====================================\n');

%disp(tabela_genes_frequentes);

end