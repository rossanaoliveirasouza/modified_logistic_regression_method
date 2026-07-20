function resultados = crossvalidation_adhoc(A, classificador, genes, K, k0, k1)

tempo_total_inicio = tic;

% Garantir que o classificador seja vetor coluna
classificador = classificador(:);

% Criar partição estratificada
cv = cvpartition(classificador, 'KFold', K);

% Inicializar vetores para guardar os resultados de todos os pacientes
probabilidades_todos = zeros(size(A,1), 1);
predicoes_todos = zeros(size(A,1), 1);

% Tempos
tempo_fold = zeros(K, 1);
tempo_predicao = zeros(K, 1);

for k = 1:K

    tempo_fold_inicio = tic;

    idx_treino = training(cv, k);
    idx_teste  = test(cv, k);

    % Separar dados
    A_treino = A(idx_treino, :);
    A_teste  = A(idx_teste, :);

    y_treino = classificador(idx_treino);
    y_teste  = classificador(idx_teste);

    % Separar treino por classe
    A0 = A_treino(y_treino == 0, :);
    A1 = A_treino(y_treino == 1, :);

    % Garantir que k0 e k1 não sejam maiores que o número de amostras disponíveis
    k0_fold = min(k0, size(A0, 1));
    k1_fold = min(k1, size(A1, 1));

    fprintf('Treino: %d pacientes sem sepse, %d com sepse\n', ...
        sum(y_treino == 0), sum(y_treino == 1));

    fprintf('Teste:  %d pacientes sem sepse, %d com sepse\n', ...
        sum(y_teste == 0), sum(y_teste == 1));

    tempo_predicao_inicio = tic;

    % Predizer cada paciente do conjunto de teste
    indices_teste = find(idx_teste);

    for i = 1:length(indices_teste)

        q = A_teste(i, :);

        [~, p] = adhoc(q, A0, A1, k0_fold, k1_fold);

        probabilidades_todos(indices_teste(i)) = p;
        predicoes_todos(indices_teste(i)) = p >= 0.5;

    end

    tempo_predicao(k) = toc(tempo_predicao_inicio);
    tempo_fold(k) = toc(tempo_fold_inicio);

end

tempo_total_cv = toc(tempo_total_inicio);

%% Avaliação final

y_real = classificador(:);
y_pred = predicoes_todos(:);
y_prob = probabilidades_todos(:);

matriz_confusao = confusionmat(y_real, y_pred, 'Order', [0 1]);

TN = matriz_confusao(1,1);
FP = matriz_confusao(1,2);
FN = matriz_confusao(2,1);
TP = matriz_confusao(2,2);

acuracia = (TP + TN) / max((TP + TN + FP + FN), eps);
sensibilidade = TP / max((TP + FN), eps);
especificidade = TN / max((TN + FP), eps);
precisao = TP / max((TP + FP), eps);
F1 = 2 * (precisao * sensibilidade) / max((precisao + sensibilidade), eps);

[Xroc, Yroc, limiares, AUC] = perfcurve(y_real, y_prob, 1);

%% Guardar resultados

resultados.acuracia = acuracia;
resultados.sensibilidade = sensibilidade;
resultados.especificidade = especificidade;
resultados.precisao = precisao;
resultados.F1 = F1;
resultados.AUC = AUC;
resultados.matriz_confusao = matriz_confusao;
resultados.y_real = y_real;
resultados.y_pred = y_pred;
resultados.y_prob = y_prob;
resultados.k0 = k0;
resultados.k1 = k1;
resultados.tempo_total_cv = tempo_total_cv;
resultados.tempo_medio_fold = mean(tempo_fold);
resultados.tempo_medio_predicao = mean(tempo_predicao);

%% Exibir resultados

fprintf('Acurácia:        %.4f\n', acuracia);
fprintf('Sensibilidade:   %.4f\n', sensibilidade);
fprintf('Especificidade:  %.4f\n', especificidade);
fprintf('Precisão:        %.4f\n', precisao);
fprintf('F1-score:        %.4f\n', F1);
fprintf('AUC:             %.4f\n', AUC);
fprintf('probabilidade:   %.4f\n', y_prob);

fprintf('\nMatriz de confusão:\n');
disp(matriz_confusao);

%fprintf('\n============================================\n');
%fprintf('TEMPOS COMPUTACIONAIS - ADHOC\n');
%fprintf('============================================\n');

%fprintf('Tempo total CV:              %.6f segundos\n', tempo_total_cv);
%fprintf('Tempo médio por fold:        %.6f segundos\n', mean(tempo_fold));
%fprintf('Tempo médio predição:        %.6f segundos\n', mean(tempo_predicao));

%figure;
%plot(Xroc, Yroc, 'LineWidth', 2);
%xlabel('Taxa de Falsos Positivos');
%ylabel('Taxa de Verdadeiros Positivos');
%title(['Curva ROC - ADHOC - AUC = ', num2str(AUC, '%.4f')]);
%grid on;

end