function [ATagli,bTagli,info]= gomoryDuale(c,A,b)
    %GOMORYDUALE Trova i tagli di Gomory per un problema in forma duale
    %standard.
    %   Il problema deve essere nella forma min c'*x che rispetta i vincoli
    %   Ax = b e x >= 0. L'algoritmo procede trovando la soluzione ottimale
    %   del problema, con la relativa base ottimale. Usa le informazioni
    %   per calcolare la matrice ATilde e calcola tutti i piani di taglio
    %   possibili seguendo le regole dell'algoritmo.
    %   PARAMETRI
    %   c: vettore colonna della funzione obiettivo
    %   A: matrice dei vincoli
    %   b: vettore colonna dei termini noti
    %   OUTPUT
    %   ATagli: matrice dei vincoli dei piani di taglio
    %   bTagli: vettore dei termini noti dei piani di taglio
    %   info: struct:
    %     - x: soluzione ottima
    %     - B: base ottima
    %     - N: indici delle variabili non di base
    %     - R: indici delle variabili di base frazionarie
    %     - ATilde: matrice ATilde = A(:,B)^-1*A(:,N)
    %     - equazioni: piani di taglio
    %   ESEMPIO
    %   c = [4; 5; 2];
    %   A = [0 0.6 0.8;
    %       -1 2 0;
    %       1 0 -1;
    %       -1 0 0;
    %       0 -1 0;
    %       0 0 -1];
    %   b = [500; 0; 0; 0; 0; 0];
    %   [ATagli, bTagli, info] = gomoryDuale(b,A',c);
    %   info.equazioni
    
    if ~exist("c", "var") || ~exist("A", "var") || ~exist("b", "var")
        error("Errore. c, A e b devono essere tutti e tre specificati.");
    end
    controlloPrimale(b,A',c);
    res = simplessoPrimale(b,A',c);
    [~,nVariabili] = size(A);
    x = res.y; % soluzione ottima duale
    B = res.base; % base ottima duale
    N = setdiff(1:nVariabili, B);
    R = find(x ~= round(x))'; % trovo indici di x frazionari
    RIndice = x(B) ~= round(x(B)); % trovo indici di x(B) frazionari
    ATilde = A(:,B)^-1*A(:,N); % calcolo ATilde
    ATagli = sym(zeros(length(R),nVariabili));
    ATagli(:,N) = -parteFrazionaria(ATilde(RIndice,:)); % calcolo ATagli
    bTagli = -parteFrazionaria(x(R)); % calcolo bTagli
    if ~isempty(ATagli) && ~isempty(bTagli)
        [ATagli,bTagli] = semplificaVincoli(ATagli,bTagli); % semplificazioni
    end
    % salvo informazioni
    info.x = x;
    info.B = B;
    info.N = N;
    info.R = R;
    info.ATilde = ATilde;
    xSym = sym('x', [nVariabili 1]);
    info.equazioni = ATagli*xSym <= bTagli;
end

