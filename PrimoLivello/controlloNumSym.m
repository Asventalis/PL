function [] = controlloNumSym(x)
% CONTROLLONUMSYM Controlla se l'input è numerico o simbolico.
%   Questa funzione controlla se l'input è una variabile numerica o una
%   variabile simbolica di tipo numerico intero o razionale.
%   PARAMETRI
%   x: variabile da controllare
    % Verifica che l'input sia presente
    if ~exist("x", "var")
        % Cattura le informazioni sul chiamante e sul nome della variabile
        st = dbstack;
        if numel(st) > 1
            caller = st(2).name;
        else
            caller = 'Base Workspace';
        end
        vname = inputname(1);
        if isempty(vname)
            vname = 'unknown';
        end
        error("Errore in %s: serve un parametro alla funzione. Input variable: %s", caller, vname);
    end

    % Verifica che l'input sia numerico o simbolico come specificato
    if ~(isnumeric(x) || (class(x)=="sym" && all(isSymType(x,"real"),"all")))
        % Cattura le informazioni sul chiamante e sul nome della variabile
        st = dbstack;
        if numel(st) > 1
            caller = st(2).name;
        else
            caller = 'Base Workspace';
        end
        vname = inputname(1);
        if isempty(vname)
            vname = 'unknown';
        end
        error("Errore in %s: L'input (%s) deve essere un tipo numerico o un tipo simbolico rappresentante un numero razionale", caller, vname);
    end
end
