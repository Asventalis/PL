function [dc,dA,db,dB] = dualeConMosse(c,A,b,B)
    %DUALECONMOSSE calcola la forma duale e la corrispettiva base di un
    %problema in forma primale con una base nel primale.
    %   Parte da un problema primale della forma max c'*x e con vincoli
    %   della forma A*x <= b e usa delle mosse per portarlo in un problema
    %   in forma duale, ovvero nella forma min dc'*x con vincoli della forma
    %   dA*x = db. Inoltre, avendo una base del problema primale,
    %   restituisce una nuova base, che sarebbe la base corrispondente alla
    %   vecchia ma nel problema duale. La base può essere omessa per
    %   evitare questo passaggio. Tale base è definita in modo che le prime
    %   componenti della soluzione di base del nuovo problema usando la
    %   nuova base siano uguali alla soluzione di base del vecchio problema
    %   usando la vecchia base.
    %   PARAMETRI
    %   c: vettore colonna della funzione obiettivo
    %   A: matrice dei vincoli
    %   b: vettore colonna dei termini noti
    %   B: una base del problema primale
    %   OUTPUT
    %   dc: funzione obiettivo del duale
    %   dA: matrice dei vincoli del duale
    %   db: vettore colonna dei termini noti del duale
    %   dB: base del duale corrispondente alla base del primale
    %   ESEMPIO
    %   c = [4; 5; 2];
    %   A = [0 0.6 0.8;
    %       -1 2 0;
    %       1 0 -1;
    %       -1 0 0;
    %       0 -1 0;
    %       0 0 -1];
    %   b = [500; 0; 0; 0; 0; 0];
    %   B = [1 2 3];
    %   [dc,dA,db,dB] = dualeConMosse(c,A,b,B)

    if ~exist("c", "var") || ~exist("A", "var") || ~exist("b", "var")
        error("Errore. c, A e b devono essere tutti e tre specificati.");
    end
    controlloPrimale(c,A,b);
    c = sym(c);
    A = sym(A);
    b = sym(b);
    [nVincoli,nVariabili] = size(A);
    if exist("B", "var")
        controlloBase(B, A);
    end

    % cambiamo di segno la funzione obiettivo
    dc = -c;

    % troviamo indici di variabili definite positive
    [lb,~] = limitiPrimale(A,b);
    posVar = find(lb>=0);

    % troviamo indici di variabili non definite positive
    undVar = setdiff(1:nVariabili, posVar);

    % per tutte le variabili non definite positive xi, applichiamo la
    % sostituzione con due variabili fittizie positive tali per cui
    % xi = xa - xb. Questo è sempre vero per ogni xi reale e xa e
    % xb positivi
    dA = [A -A(:, undVar)];
    dc = [dc; -dc(undVar)];

    % rimuoviamo vincoli della forma -xi <= 0 perché equivalgono a
    % dire che la variabile è definita positiva, ma è un vincolo
    % ridondante nel problema duale, in cui le variabili sono tutte
    % definite positive implicitamente
    vincoliRidondanti = ismember(A,-eye(nVariabili),"rows") & b==0;
    dA(vincoliRidondanti,:) = [];
    db = b;
    db(vincoliRidondanti) = [];
    nNuoviVincoli = size(dA,1);
    
    % trasforma tutte le disuguaglianze in uguaglianze aggiungendo
    % a ogni disequazione una variabile di scarto positiva.
    dA = [dA eye(nNuoviVincoli)];
    dc = [dc; zeros(nNuoviVincoli,1)];

    if ~exist("B","var")
        return;
    end
    % Scegliere la nuova base è un processo di eliminazione.
    % Abbiamo: 
    % variabili originali + 
    % variabili fittizie per numeri non definiti positivi + 
    % variabili di scarto per ogni vincolo non ridondante
    % La base del duale deve essere lunga quanto i vincoli. I vincoli del
    % duale sono pari al numero di vincoli non ridondanti del primale.
    % Quindi dobbiamo sottrarre, da tutte le variabili del duale,
    % esattamente un numero di variabili pari a:
    % variabili originali + 
    % variabili fittizie per numeri non definiti positivi
    % Le variabili le rimuoveremo in base a due criteri:
    % 1. Per ogni variabile non definita positiva, abbiamo due variabili
    %    positive. Sia per esempio x = a - b una di esse. Allora si 
    %    dimostra che possiamo sempre rimuovere almeno una di esse.
    %    Infatti, immaginiamo che x sia il valore assunto nella soluzione
    %    di base del primale della variabile non definita positiva. Allora:
    %    a. x >= 0, allora ipotizzeremo che x = a e b = 0 e scarteremo b. 
    %    b. x < 0, allora ipotizzeremo che x = -b e a = 0, e quindi 
    %       scarteremo a. 
    % 2. Per ogni vincolo di base nel primale, possiamo rimuovere una
    %    variabile nel duale. Distinguiamo tra vincoli ridondanti e vincoli
    %    non ridondanti:
    %    a. Se il vincolo è ridondante ed è in base, vuol dire che è
    %       rispettato con l'uguale e che corrisponde con il fatto che una
    %       variabile originale sia nulla e quindi possiamo levarla dalla
    %       base duale.
    %    b. Se il vincolo non è ridondante, ed è in base, allora è
    %       rispettato con l'uguale. Ciò implica che, nel corrispettivo 
    %       duale, la variabile di scarto deve essere nulla.
    % Intanto dobbiamo scegliere per ogni variabile non definita positiva
    % se prendere la colonna della parte positiva o della parte
    % negativa, e questo dipende dal valore della variabile nella
    % soluzione di base: se sappiamo che viene positiva, allora
    % dovremo prendere la colonna positiva, ovvero i, mentre se
    % sappiamo che viene negativa, allora dovremo prendere la
    % colonna negativa, ovvero n_var+i.

    %1.
    dB = 1:nVariabili+length(undVar)+nNuoviVincoli;
    xOld = A(B,:)^-1*b(B);
    daRimuovere = zeros(1,length(undVar));
    for i=1:length(undVar)
        if xOld(i)>=0 % 1.a.
            daRimuovere(i) = nVariabili + i; % rimuoviamo la colonna negativa
        else % 1.b.
            daRimuovere(i) = i; % rimuoviamo la colonna positiva
        end
    end
    dB(daRimuovere) = [];
    % 2.a
    % adesso invece dobbiamo controllare se nella base vecchia ci
    % sono delle righe ridondanti. Ciò equivale a dire che la variabile
    % corrispondente è uguale a 0 nella soluzione di base, e quindi nel
    % duale ciò vuol dire che la variabile è fuori dalla base.
    % Ciò implica che dobbiamo rimuovere la variabile corrispondente al
    % vincolo ridondante dalla base nel duale.
    % Per fare ciò, calcolo le variabili corrispondenti ai vincoli
    % ridondanti come la posizione che ogni vincolo ridondante ha
    % all'interno della matrice identità.
    [~,varCorrisp] = ismember(A(vincoliRidondanti,:),-eye(nVariabili),"rows"); 
    % creo un vettore di coppie (a,b) con a indice del vincolo ridondante e
    % b indice della variabile corrispondente
    vinRidVar = [find(vincoliRidondanti) varCorrisp];
    % trovo i vincoli ridondanti che sono in base vecchia
    daRimuovere = ismember(vinRidVar(:,1),B);
    % rimuovo le variabili corrispondenti dalla base nuova
    dB = setdiff(dB, vinRidVar(daRimuovere,2));
    
    % 2.b
    % L'ultima cosa che ci resta da fare è rimuovere le variabili di scarto
    % corrispondenti ai vincoli rimanenti nella base vecchia, escludendo
    % però i vincoli ridondanti che abbiamo già considerato. Come prima
    % cosa rimuovo quindi i vincoli ridondanti dalla base vecchia
    B_non_rid = setdiff(B, vinRidVar(:,1));
    % dopodiché mi calcolo la nuova posizione dei vecchi vincoli della base, che è
    % cambiata a seguito della rimozione dei vincoli ridondanti.
    vecchiVincoli = 1:nVincoli;
    nuoviVincoli = vecchiVincoli;
    nuoviVincoli(vincoliRidondanti) = []; % base vecchia a cui rimuoviamo i vincoli ridondanti
    % sappiamo che nuova_pos non contiene valori nulli perché abbiamo già
    % rimosso i vincoli ridondanti sia dalla base che da nuoviVincoli.
    [~,nuova_pos] = ismember(B_non_rid,nuoviVincoli); % nuova_pos = posizione di B_i in nuoviVincoli
    % le variabili di scarto vengono in ordine, però vengono dopo le
    % variabili già esistenti e le variabili fittizie per garantire la
    % positività. Quindi, dato un vincolo i, la sua variabile di scarto
    % sarà nella posizione i+n_var+n_und_var
    nuova_pos = nuova_pos + (length(undVar)+nVariabili)*ones(1,length(nuova_pos));
    % rimuoviamo queste variabili dalla base e abbiamo finito.
    dB = setdiff(dB, nuova_pos);
end

