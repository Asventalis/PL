function [] = controlloPrimale(c,A,b)
    %CONTROLLOPRIMALE Lancia un errore se il problema lineare c'*x con
    %vincoli A*x<=b non è in formato primale standard.
    %   Controlla che i vettori c, A e b possano formare un problema di
    %   programmazione lineare in formato primale standard. Per fare ciò, c
    %   e b dovranno essere vettori colonna, e A dovrà avere come numero di
    %   vincoli (righe) la lunghezza di b e come numero di variabili
    %   (colonne) la lunghezza di c. Ovviamente tutte e tre le variabili
    %   dovranno essere degli array non vuoti di tipo numerico o simbolico.
    %   PARAMETRI
    %   c: vettore della funzione obiettivo
    %   A: matrice dei vincoli
    %   b: vettore dei termini noti
    
    if ~exist("c","var") || ~exist("A","var") || ~exist("b","var")
        error("Errore. Tutti gli argomenti vanno passati.");
    end
    controlloAb(A,b);
    controlloAc(A,c);
end

