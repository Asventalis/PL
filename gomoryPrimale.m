function [ATagli,bTagli,info] = gomoryPrimale(c,A,b)
    %GOMORYPRIMALE Trova i tagli di Gomory per un problema in forma
    %primale standard. Funziona solo se il problema ha TUTTE LE VARIABILI
    %DEFINITE POSITIVE.
    %   L'algoritmo trasforma prima il primale in un duale usando
    %   dualeConMosse. Applica poi gomoryDuale al problema trasformato, e
    %   poi trasforma i piani di taglio in modo che siano espressi solo in
    %   funzione delle variabili del problema primale originario.
    %   Quest'ultima operazione è eseguita ipotizzando che tutte le
    %   variabili del problema in forma primale standard siano definite
    %   positive. L'algoritmo non funzionerà se ciò non è vero, e darà
    %   errore. Non esiste un modo per passare da un primale con variabili
    %   non tutte definite positive a un duale utilizzando le mosse
    %   standard, poi aggiungere i piani di taglio, e poi tornare dal duale
    %   nuovo al primale, perché aggiungendo dei nuovi vincoli arbitrari si
    %   perde la biunivocità della trasformazione con mosse. Questa
    %   affermazione si basa sul fatto che i piani di taglio siano
    %   arbitrari rispetto alla trasformazione con mosse. Ciò non è
    %   necessariamente vero.
    %   PARAMETRI
    %   c: vettore colonna della funzione obiettivo
    %   A: matrice dei vincoli
    %   b: vettore colonna dei termini noti
    %   OUTPUT
    %   ATagli: matrice dei vincoli dei piani di taglio nel primale
    %   bTagli: vettore dei termini noti dei piani di taglio nel primale
    %   info: struct:
    %     - cDuale: funzione obiettivo del problema in forma duale
    %     - ADuale: matrice dei vincoli del problema in forma duale
    %     - bDuale: termini noti del problema in forma duale
    %     - x: soluzione ottima del problema in forma duale
    %     - B: indici delle variabili della base ottima
    %     - N: indici delle variabili non di base ottima
    %     - R: indici delle variabili di base ottima non intere
    %     - ATilde: matrice ATilde = A(:,B)^-1*A(:,N)
    %     - ATagliDuale: matrice dei vincoli dei tagli del duale
    %     - bTagliDuale: vettore dei termini noti dei tagli del duale
    %     - equazioniDuale: piani di taglio del problema in forma duale
    %     - equazioni: piani di taglio del problema in forma primale
    %   ESEMPIO
    %   c = [4; 5; 2];
    %   A = [0 0.6 0.8;
    %       -1 2 0;
    %       1 0 -1;
    %       -1 0 0;
    %       0 -1 0;
    %       0 0 -1];
    %   b = [500; 0; 0; 0; 0; 0];
    %   [ATagli, bTagli, info] = gomoryPrimale(c,A,b);
    %   info.equazioniDuale
    %   info.equazioni
    
    if ~exist("c", "var") || ~exist("A", "var") || ~exist("b", "var")
        error("Errore. c, A e b devono essere tutti presenti come parametri.");
    end
    controlloPrimale(c,A,b);
    [lb,~] = limitiPrimale(A,b);
    if any(lb<0)
        error("Errore. Tutte le variabili devono essere definite positive.")
    end
    [cd,Ad,bd] = dualeConMosse(c,A,b); % conversione del problema da primale a duale con mosse
    [ATagliDuale,bTagliDuale,infoDuale] = gomoryDuale(cd,Ad,bd); % calcolo equazioni del duale
    [~,nVariabili] = size(A);
    [~,nVariabiliDuale] = size(Ad);
    ATagli = ATagliDuale(:,1:nVariabili) - ATagliDuale(:,nVariabili+1:nVariabiliDuale)*Ad(:,1:nVariabili); % sostituzioni per ATagli
    bTagli = bTagliDuale - ATagliDuale(:,nVariabili+1:nVariabiliDuale)*bd; % sostituzioni per bTagli
    if ~isempty(ATagli) && ~isempty(bTagli)
        [ATagli,bTagli] = semplificaVincoli(ATagli,bTagli); % semplificazioni
    end
    % salvo informazioni
    info.cDuale = cd;
    info.ADuale = Ad;
    info.bDuale = bd;
    info.x = infoDuale.x;
    info.B = infoDuale.B;
    info.N = infoDuale.N;
    info.R = infoDuale.R;
    info.ATilde = infoDuale.ATilde;
    info.ATagliDuale = ATagliDuale;
    info.bTagliDuale = bTagliDuale;
    info.equazioniDuale = infoDuale.equazioni;
    xSym = sym("x",[nVariabili 1]);
    info.equazioni = ATagli*xSym <= bTagli;
end

