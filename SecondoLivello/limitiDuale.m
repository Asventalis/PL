function [lb,ub] = limitiDuale(A,c)
% LIMITIDUALE Calcola l'estremo inferiore e l'estremo superiore che
% ogni variabile yi può assumere rispettando i vincoli y'*A=c'.
%   Usa l'algoritmo del simplesso duale considerando come b i vettori
%   della base canonica e anche i loro opposti, per trovare i valori
%   minimi e massimi assumibili da ogni coordinata del vettore y.
%   PARAMETRI
%   A: matrice dei vincoli
%   c: funzione obiettivo
%   OUTPUT
%   lb: vettore dei limiti inferiori
%   ub: vettore dei limiti superiori
%   ESEMPIO
%   c = [4; 5; 2];
%   A = [0 0.6 0.8;
%       -1 2 0;
%       1 0 -1;
%       -1 0 0;
%       0 -1 0;
%       0 0 -1];
%   b = [500; 0; 0; 0; 0; 0];
%   [lb, ub] = limitiDuale(A,c);
%   xSym = sym("x", [6 1]);
%   expr = lb <= xSym & xSym <= ub
    if ~exist("A", "var") || ~exist("c", "var")
        error("Errore. A e c devono essere tutti e due presenti.");
    end
    controlloAc(A,c);
    [nVariabili,~] = size(A);
    lb = sym(zeros(nVariabili,1));
    ub = sym(zeros(nVariabili,1));
    versori = eye(nVariabili);
    for i = 1:nVariabili
        b = versori(:,i);
        res = simplessoDuale(c,A,b);
        lb(i) = res.y(i);
        b = -versori(:,i);
        res = simplessoDuale(c,A,b);
        ub(i) = res.y(i);
    end
end

