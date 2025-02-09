% Esame 29/01/2024

c = [4; 5; 2];
A = [0 0.6 0.8;
    -1 2 0;
    1 0 -1;
    -1 0 0;
    0 -1 0;
    0 0 -1];
b = [500; 0; 0; 0; 0; 0];
B = [1 4 5];
% fare un passo del simplesso partendo da B
resP = simplessoPrimale(c,A(1:3,:),b(1:3,:));
stampaStruct(resP.steps);
% calcolare il primo taglio di Gomory
% [ATagli,bTagli,info] = gomoryPrimale(c,A,b);
% info.x'
% info.B
% info.ATilde
% ATagli
% bTagli
% info.equazioniDuale
% info.equazioni
% % utilizzare il primo taglio per vedere se siamo arrivati all'ottimo PLI
% ANuovo = [A; ATagli(1,:)];
% bNuovo = [b; bTagli(1)];
% resP = simplessoPrimale(c,ANuovo,bNuovo);
% stampaStruct(resP.steps);

% PROVE
[lb, ub] = limitiPrimale(A(1:3,:),b(1:3,:))
[cd,Ad,bd,Bd] = dualeConMosse(c,A(1:3,:),b(1:3,:), resP.base) % conversione del problema da primale a duale con mosse
% risolvi con il simplesso
resD = simplessoPrimale(bd,Ad',cd)
stampaStruct(resD.steps)
[ATagliDuale,bTagliDuale,infoDuale] = gomoryDuale(cd,Ad,bd) % calcolo equazioni del duale
visualizzaPoliedro3D(c,A(1:4,:),b(1:4,:))