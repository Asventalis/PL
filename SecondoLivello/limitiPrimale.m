function [lb,ub] = limitiPrimale(A,b)
% LIMITIPRIMALE Calcola l'estremo inferiore e l'estremo superiore che
% ogni variabile xi può assumere rispettando i vincoli A*x<=b.
%   Usa l'algoritmo del simplesso primale considerando come c i vettori
%   della base canonica e anche i loro opposti, per trovare i valori
%   massimi e minimi assumibili da ogni coordinata del vettore x.
%   PARAMETRI
%   A: matrice dei vincoli
%   b: vettore dei termini noti
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
%   [lb, ub] = limitiPrimale(A,b);
%   xSym = sym("x", [3 1]);
%   expr = lb <= xSym & xSym <= ub
    if ~exist("A", "var") || ~exist("b", "var")
        error("Errore. A e b devono essere tutti e due presenti.");
    end
    controlloAb(A,b);
    [~,nVariabili] = size(A);
    lb = sym(zeros(nVariabili,1));
    ub = sym(zeros(nVariabili,1));
    versori = eye(nVariabili);
    for i = 1:nVariabili
        c = -versori(:,i);
        res = simplessoPrimale(c,A,b);
        lb(i) = res.x(i);
        c = versori(:,i);
        res = simplessoPrimale(c,A,b);
        ub(i) = res.x(i);
    end
end

