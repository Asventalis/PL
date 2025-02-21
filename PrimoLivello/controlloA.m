function [] = controlloA(A)
% CONTROLLOA Lancia un errore se la matrice A non può essere una 
% matrice dei vincoli per un problema di programmazione lineare in
% formato primale standard.
%   Controlla le seguenti cose:
%    - A è passato come argomento
%    - A è non vuoto
%    - A è numerico o simbolico
%    - Il numero di vincoli deve essere almeno la dimensione dello
%      spazio in cui siamo perché altrimenti il problema non avrebbe
%      vertici e sarebbe dunque banale.
%    - A deve essere una matrice non singolare, altrimenti tutte le
%      basi avrebbero matrice di base non invertibile, ovvero non
%      esisterebbero vertici.
%   PARAMETRI
%   A: matrice dei vincoli
    if ~exist("A","var")
        error("Errore. La matrice dei vincoli A va passata come argomento.");
    end
    if isempty(A)
        error("Errore. La matrice A è vuota.");
    end
    controlloNumSym(A);
    [nVincoli, nVariabili] = size(A);
    if nVincoli < nVariabili
        error("Errore. Il numero dei vincoli deve essere maggiore o uguale del numero delle variabili.");
    end
    if rank(A) ~= nVariabili
        error("Errore. A non ha rango massimo. Ciò vuol dire che non esistono basi per A.")
    end
end

