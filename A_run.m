function A_run(file,times,zeta,type,theta,nu,vcon,out)
% A_AnIso(file,times,zeta,type,theta,nu,vcon,out)
% file - starting mesh
% times - how often
% zeta - Zulaessigkeitsbed
% type - art des Tests
% theta - alle verfeinern oder nur wichtige? 
% nu - wie verfeinern iso oder aniso
% vcon - soll Vorkonditioniert werden? 0 1
% out - Dateiname der Ausgabe
%
% P. Schaefer

%Parameter testen
% assert(length(zeta)==max(type)-1||length(zeta)==1,...
%   'Pleas set right type and zeta parameters');

%Netz laden
A_loadMeshF(file)

%globale Speicherung initialisieren
global G_T;
global G_D;

%Schrittweise Verfeinern
for i = 1:times
    
 %benötigte Zeit Schätzen
 if(size(G_T,1)>2)
  nextF = G_T(end,1)*4;
  nextS = interp1(1:size(G_T,1),G_T(:,1)',size(G_T,1)+1,'slpine');
  cald = size(G_T,1)+(-1:0);
  calc = size(G_T,1)+(-2:0);
  nextTime = [nextS interp1(G_T(calc,1)',G_T(calc,2)',nextF,'spline') ...
    interp1(G_T(cald,1)',G_T(cald,3)',nextF,'spline') ...
    interp1(G_T(calc,1)',G_T(calc,4)',nextS,'spline')];
 end
 
 % Ein kompletter Verfeinerungschritt
 A_step(zeta,type,theta,nu,vcon,out);
 
 %Zeit Speichern
 usedTime = G_T(end,:);
 
 %Ergebnisse Abspeichern
 data = G_D(end,:)
 typeN = int2str(type);
 A_saveMesh([out typeN(typeN~=' ') '_' int2str(size(G_T,1))]);
 
end

