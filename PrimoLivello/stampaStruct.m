function [] = stampaStruct(steps)
% STAMPASTRUCT Stampa uno struct array usando una table.
%   Usa i fieldnames dello struct array come nomi delle colonne. usa la
%   funzione stampaInLinea per trasformare il valore in una stringa,
%   poi usa categorical per fare in modo che non spuntino le virgolette
%   nella stampa (le stringhe vengono stampate con "" attorno, quindi
%   si converte la stringa a un oggetto categorical, che è una cosa
%   simile a un enum in C, e non fa spuntare le ""). La funzione
%   stampaInLinea gestisce i casi in cui la variabile è NaN, i casi in
%   cui è un vettore riga o colonna, e i casi in cui è una matrice.
    fn = fieldnames(steps);
    stringSteps(length(steps)) = struct;
    for k = 1:length(steps)
        for h = 1:length(fn)
            field = fn{h};
            stringSteps(k).(field) = categorical(stampaInLinea(steps(k).(field)));
        end
    end
    tab = struct2table(stringSteps);
    disp(tab);
end

