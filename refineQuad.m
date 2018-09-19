function [coo,ele,nei,f2s,sit,err] = refineQuad(coordinates,elements,neigh,sites,typ)
%
%  [coo,ele,nei,f2s,sit,err] = refineQuad(coordinates,elements,type)
%
% Verfeinert die markierten Elemente mit dem entsprechenden TYP und gibt
% auch die F2S Beziehungen zurueck. type muss von der Laenge der Anzahl der
% Elemente entsprechen oder genau 1 und die Eintraege koennen 1,2,3,4,5 sein.
%
% der Typ zu jedem Element entspricht dabei:
% 1 - keine Verfeinerung
% 2 - 4 neue Elemente
% 3 - 2 neue Elemente, uebereinander
% 4 - 2 neue Elemente, nebeneinander
% 5 - 4 neue Elemente, wird erst (3) und dann beide Elemente (4) geteilt
%
% P. Schaefer

err = [];

%Type: wenn nur ein Wert: aufblaehen
if([1 1] == size(typ))
    typ = repmat(typ, size(elements,1),1);
end

assert(size(elements,1)==size(neigh,1)&&size(elements,1)==length(typ),...
  'Dimensionen passen nicht');


%Globale Variabelen aufbauen
global G_ref_E;
global G_ref_C;
global G_ref_N;
global G_ref_S;
global G_ref_f2s; %finale Beziehung (VaterSohn)
global G_ref_t;  %wie soll verfeinert werden
global G_ref_tD; %wie wurde bereits verfeinert (in diesem Durchlauf)
global G_ref_f2sT; %temporaere Beziehung

%INTERNE globale Variablen zuweisen
G_ref_E = elements;
G_ref_C = coordinates;
G_ref_N = neigh;
G_ref_S = sites;
G_ref_t = typ;
G_ref_f2s = repmat((1:size(elements,1))',1,4);
G_ref_tD = ones(size(elements,1),1);

%Parameter freigeben (Speicher...)
clear elements coordinates neigh typ

% figure(11)
% plotShape(G_ref_C,G_ref_E,'s',G_ref_t);
% view(2)
% colorbar
% title('zuVerfeinern')

ref_old = [];
ref_old2 = [];

%jedes Element verfeinern
while(1==1)
  %jeden vierten 2er durch schrittweise Verfeinerung (5) ersetzen
  t_ref=find(G_ref_t==2);
  G_ref_t(t_ref(2:4:end)) = 5;
  
  %welche Elemente muessen bearbeitet werden
  ref = find(G_ref_t>1);
  ref = reshape(ref,1,length(ref));
  
  %muss noch weiter Verfeinert werden?
  if(isequal(ref,ref_old))
      assert(~isequal(ref_old2,G_ref_t,'Markierte sind verschieden'));
    break;
  end
  ref_old = ref;
  ref_old2  = G_ref_t;
  if(isempty(ref))
    break;
  end
  
  % Elementeweise bearbeiten
  for ele = ref % ref(randperm(length(ref)))

    % HangingNode Check
    Nt = find(G_ref_N(ele,5:8)==0);
    N = G_ref_N(ele,Nt);
    N2t = find(N~=0);
    N2 = N(N2t); %Nachbarn der Kanten mit nur einem Nachbar
    N2tt = Nt(N2t);  %Kante mit Nachbar ^
    
    % Nur Nachbarn verfeinern an Kanten die geteilt werden
    if(G_ref_t(ele)==3 || G_ref_t(ele)==5)
       N2 = N2(mod(N2tt,2)==0);
    elseif(G_ref_t(ele)==4)
       N2 = N2(mod(N2tt,2)==1);
    end

    %hat noch zu teilende Nachbarn?
    if(~isempty(N2))
      N3t = mod(find((G_ref_N(N2',:)==ele)')-1,4)+1; %Nachbarseiten    
      N4t = find(diag(G_ref_N(N2',(N3t + 4)'))~=0)';  %Nachbarn mit 2Nachbarn
        if(~isempty(N4t))
          for i = N4t  %Elemente die noch verfeinert werden muessen
            assert(G_ref_tD(N2(i))~=2,'Element wurde schon voll verfeinert')
          
            
            %Element ist schon markiert
            if(G_ref_t(N2(i)) == 2 || G_ref_t(N2(i)) ==5)
                continue;
            end
            
            %Wie muss das Element verfeinert werden.
            if(G_ref_tD(N2(i)) ~= 1)
              assert(G_ref_tD(N2(i))~=mod(N3t(i),2)+3,...
                  'Element wurde in der Richtung schon verfeinert')
              G_ref_t(N2(i)) = mod(N3t(i),2)+3;
            elseif(G_ref_t(N2(i)) == 1)
              G_ref_t(N2(i))=mod(N3t(i),2)+3;
            elseif(G_ref_t(N2(i)) == mod(N3t(i)+1,2)+3)
              G_ref_t(N2(i)) = 2;     
            end
            
          end
          
          % Da Nachbarn noch verfeinert werden muessen, erst mal weiter
          continue;
        end
    end

    % Wenn alle Ueberpruefungen durchgelaufen sind
     assert(G_ref_tD(ele)~=2,'Element ist schon verfeinert')
     assert(G_ref_t(ele)>1,'Element ist nicht Markiert')
     G_ref_f2sT = ones(1,4)*ele;
     refineE(ele); %Element teilen
     updateN(ele); %Nachbarn des Elements aktualisieren
     updateF2S(ele); %VaterSohn Relation setzen
  end
end

%Rueckgabe zuweisen
coo = G_ref_C;
ele = G_ref_E;
nei = G_ref_N;
f2s = G_ref_f2s;
sit = G_ref_S;

%doppelte Koordinaten loeschen
[coo , ~, pos] = unique(coo,'rows');
pos = pos';
ele = pos(ele);

 %INTERNE globale Variablen freigeben
clear G_ref_E G_ref_C G_ref_N G_ref_f2s G_ref_t G_ref_tD G_ref_s
end

%% Element verfeinern ! sollte nur ausgefuehrt werden wenn wirklich moeglich
function refineE(ele)
% Element wird gnadenlos Verfeinert

global G_ref_E;
global G_ref_C;
global G_ref_S;
global G_ref_f2sT;
global G_ref_t;

    c_ele = size(G_ref_E,1);
    c_coo = size(G_ref_C,1);
    
    
    el = G_ref_E(ele,:);
    if(G_ref_t(ele)==2)
        G_ref_C(c_coo+1,:) = (G_ref_C(el(1),:)+G_ref_C(el(2),:))/2;
        G_ref_C(c_coo+2,:) = (G_ref_C(el(2),:)+G_ref_C(el(3),:))/2;
        G_ref_C(c_coo+3,:) = (G_ref_C(el(3),:)+G_ref_C(el(4),:))/2;
        G_ref_C(c_coo+4,:) = (G_ref_C(el(1),:)+G_ref_C(el(4),:))/2;
        G_ref_C(c_coo+5,:) = (G_ref_C(el(1),:)+G_ref_C(el(3),:))/2;

        G_ref_E(ele,:) = [c_coo+4,c_coo+5,c_coo+3,el(4)];
        G_ref_E(c_ele+1,:) = [c_coo+5,c_coo+2,el(3),c_coo+3];
        G_ref_E(c_ele+2,:) = [c_coo+1,el(2),c_coo+2,c_coo+5];
        G_ref_E(c_ele+3,:) = [el(1),c_coo+1,c_coo+5,c_coo+4];
        
        G_ref_S([ele (c_ele+1):(c_ele+3)]',:) = repmat(G_ref_S(ele,:)+[1 1],4,1);
        
        G_ref_f2sT(1:3)=c_ele+3:-1:c_ele+1;
    elseif(G_ref_t(ele)==3||G_ref_t(ele)==5)
        G_ref_C(c_coo+1,:) = (G_ref_C(el(1),:)+G_ref_C(el(4),:))/2;
        G_ref_C(c_coo+2,:) = (G_ref_C(el(2),:)+G_ref_C(el(3),:))/2;
        
        G_ref_E(ele,1) = c_coo+1;
        G_ref_E(ele,2) = c_coo+2;        
        G_ref_E(c_ele+1,:) = [el(1),el(2),c_coo+2,c_coo+1];
        
        G_ref_S([ele (c_ele+1)]',:) = repmat(G_ref_S(ele,:)+[0 1],2,1);
        
        G_ref_f2sT([1 2])=c_ele+1;
     elseif(G_ref_t(ele)==4)
        G_ref_C(c_coo+1,:) = (G_ref_C(el(1),:)+G_ref_C(el(2),:))/2;
        G_ref_C(c_coo+2,:) = (G_ref_C(el(4),:)+G_ref_C(el(3),:))/2;
        
        G_ref_E(ele,2) = c_coo+1;
        G_ref_E(ele,3) = c_coo+2;       
        G_ref_E(c_ele+1,:) = [c_coo+1,el(2),el(3),c_coo+2];
        
        G_ref_S([ele (c_ele+1)]',:) = repmat(G_ref_S(ele,:)+[1 0],2,1);
        
        G_ref_f2sT([2 3])=c_ele+1;
    end
end


%% Aktualisieren der Nachbarn ! sollte nur ausgefuehrt werden wenn wirklich moeglich
function updateN(ele)
% Nachbarschaften werden neu gesetzt (nach N und f2s)

global G_ref_E;
global G_ref_N;
global G_ref_f2sT;
global G_ref_t;

this = G_ref_N(ele,:);
%an welchen Kanten habe ich Nachbarn
S = find(mod((this(1:4)~=0).*(this(5:8)==0),2))'; %einen Nachbar (Single)
D = find(this(5:8)~=0)';    %zwei Nachbarn (Double)

%an welchen Kanten bin ich Nachbar
MSt = mod(find((G_ref_N(this(S),:)==ele)')-1,8)+1;
MS = mod(MSt-1,4)+1;  % (Single)    %OHNE MOD????
MD = mod(find((G_ref_N(this([D D+4]),:)==ele)')-1,4)+1;
MD = reshape(MD,length(MD)/2,2); % (Double)

G_ref_N(G_ref_f2sT,1:8) = 0;

if(G_ref_t(ele)==3||G_ref_t(ele)==5)
  G_ref_N(G_ref_f2sT([1 3])',1:4) = [ 0 0 G_ref_f2sT(3) 0;G_ref_f2sT(1) 0 0 0];    %Innere Beziehung
  
      % Beziehungen fuer Kanten mit einem Nachbar
    for i = 1:length(S)
      if(mod(S(i),2)==0)
        G_ref_N(this(S(i)),[MS(i) MS(i)+4]) = [G_ref_f2sT(S(i)) G_ref_f2sT(mod(S(i),4)+1)];
        G_ref_N([G_ref_f2sT(S(i)) G_ref_f2sT(mod(S(i),4)+1)]',S(i))=this(S(i));
      else
        G_ref_N(this(S(i)),MSt(i)) = G_ref_f2sT(S(i));
        G_ref_N(G_ref_f2sT(S(i)),S(i)) = this(S(i));
      end
    end
  
    % Beziehungen fuer Kanten mit zwei Nachbarn
    for i = 1:length(D)
      if(mod(D(i),2)==0)
        if(length(unique([G_ref_E(this(D(i)),:) G_ref_E(G_ref_f2sT(D(i)),:)]))==7)
          G_ref_N(this(D(i)),MD(i,1)) = G_ref_f2sT(D(i));
          G_ref_N(this(D(i)+4),MD(i,2)) = G_ref_f2sT(mod(D(i),4)+1);
          G_ref_N(G_ref_f2sT(D(i)),D(i)) = this(D(i));
          G_ref_N(G_ref_f2sT(mod(D(i),4)+1),D(i)) = this(D(i)+4);
        else
          G_ref_N(this(D(i)),MD(i,1)) = G_ref_f2sT(mod(D(i),4)+1);
          G_ref_N(this(D(i)+4),MD(i,2)) = G_ref_f2sT(D(i));
          G_ref_N(G_ref_f2sT(D(i)),D(i)) = this(D(i)+4);
          G_ref_N(G_ref_f2sT(mod(D(i),4)+1),D(i)) = this(D(i));      
        end
      else
        G_ref_N(this(D(i)),MD(i,1)) = G_ref_f2sT(D(i));
        G_ref_N(this(D(i)+4),MD(i,2)) = G_ref_f2sT(D(i));
        G_ref_N(G_ref_f2sT(D(i)),[D(i) D(i)+4]) = [this(D(i)) this(D(i)+4)];
      end
    end
  
elseif(G_ref_t(ele)==4)
  G_ref_N(G_ref_f2sT([1 2])',1:4) = [ 0 G_ref_f2sT(2) 0 0;0 0 0 G_ref_f2sT(1)];
    % Beziehungen fuer Kanten mit einem Nachbar
    for i = 1:length(S)
      if(mod(S(i),2)==1)
        G_ref_N(this(S(i)),[MS(i) MS(i)+4]) = [G_ref_f2sT(S(i)) G_ref_f2sT(mod(S(i),4)+1)];
        G_ref_N([G_ref_f2sT(S(i)) G_ref_f2sT(mod(S(i),4)+1)]',S(i))=this(S(i));
      else
        G_ref_N(this(S(i)),MSt(i)) = G_ref_f2sT(S(i));
        G_ref_N(G_ref_f2sT(S(i)),S(i)) = this(S(i));
      end
    end
  
    % Beziehungen fuer Kanten mit zwei Nachbarn
    for i = 1:length(D)
      if(mod(D(i),2)==1)
        if(length(unique([G_ref_E(this(D(i)),:) G_ref_E(G_ref_f2sT(D(i)),:)]))==7)
          G_ref_N(this(D(i)),MD(i,1)) = G_ref_f2sT(D(i));
          G_ref_N(this(D(i)+4),MD(i,2)) = G_ref_f2sT(mod(D(i),4)+1);
          G_ref_N(G_ref_f2sT(D(i)),D(i)) = this(D(i));
          G_ref_N(G_ref_f2sT(mod(D(i),4)+1),D(i)) = this(D(i)+4);
        else
          G_ref_N(this(D(i)),MD(i,1)) = G_ref_f2sT(mod(D(i),4)+1);
          G_ref_N(this(D(i)+4),MD(i,2)) = G_ref_f2sT(D(i));
          G_ref_N(G_ref_f2sT(D(i)),D(i)) = this(D(i)+4);
          G_ref_N(G_ref_f2sT(mod(D(i),4)+1),D(i)) = this(D(i));      
        end
      else
        G_ref_N(this(D(i)),MD(i,1)) = G_ref_f2sT(D(i));
        G_ref_N(this(D(i)+4),MD(i,2)) = G_ref_f2sT(D(i));
        G_ref_N(G_ref_f2sT(D(i)),[D(i) D(i)+4]) = [this(D(i)) this(D(i)+4)];
      end
    end
  
elseif(G_ref_t(ele)==2)
  G_ref_N(G_ref_f2sT',1:4) = [0 G_ref_f2sT(2) G_ref_f2sT(4) 0; 0 0 G_ref_f2sT(3) G_ref_f2sT(1);...
                         G_ref_f2sT(2) 0 0 G_ref_f2sT(4); G_ref_f2sT(1) G_ref_f2sT(3) 0 0];
  
    % Beziehungen fuer Kanten mit einem Nachbar
    for i = 1:length(S)
      G_ref_N(this(S(i)),[MS(i) MS(i)+4]) = [G_ref_f2sT(S(i)) G_ref_f2sT(mod(S(i),4)+1)];
      G_ref_N([G_ref_f2sT(S(i)) G_ref_f2sT(mod(S(i),4)+1)]',S(i))=this(S(i));
    end
  
    % Beziehungen fuer Kanten mit zwei Nachbarn
    for i = 1:length(D)
      if(length(unique([G_ref_E(this(D(i)),:) G_ref_E(G_ref_f2sT(D(i)),:)]))==7)
        G_ref_N(this(D(i)),MD(i,1)) = G_ref_f2sT(D(i));
        G_ref_N(this(D(i)+4),MD(i,2)) = G_ref_f2sT(mod(D(i),4)+1);
        G_ref_N(G_ref_f2sT(D(i)),D(i)) = this(D(i));
        G_ref_N(G_ref_f2sT(mod(D(i),4)+1),D(i)) = this(D(i)+4);
      else
        G_ref_N(this(D(i)),MD(i,1)) = G_ref_f2sT(mod(D(i),4)+1);
        G_ref_N(this(D(i)+4),MD(i,2)) = G_ref_f2sT(D(i));
        G_ref_N(G_ref_f2sT(D(i)),D(i)) = this(D(i)+4);
        G_ref_N(G_ref_f2sT(mod(D(i),4)+1),D(i)) = this(D(i));      
      end
    end
end


end

%% aktualisieren der VaterSohn Beziehung ! sollte nur ausgefuehrt werden wenn wirklich moeglich
function updateF2S(ele)
%VaterSohn Beziehungen richtig Setzen
    global G_ref_f2s;
    global G_ref_f2sT;
    global G_ref_t;
    global G_ref_tD;
   
    if(G_ref_tD(ele)==1) %wenn Element zum ersten Mal verfeinert wird
        G_ref_f2s(ele,:) = G_ref_f2sT;
        if(G_ref_t(ele)<5)
            G_ref_tD(G_ref_f2sT) = G_ref_t(ele);
            G_ref_t(G_ref_f2sT) = 0;
        else %Element muss noch mal geteilt werden
            G_ref_tD(G_ref_f2sT) = 3;
            G_ref_t(G_ref_f2sT) = 4;
        end
    else %wenn Element zum zweiten Mal verfeinert wird
        G_ref_tD(G_ref_f2sT) = 2;
        org=floor((find(G_ref_f2s'==ele,1)-1)/4)+1;
        pos=find(G_ref_f2s(org,:)==ele,1);
        
        if(G_ref_t(ele)==3)
            if(pos==1)
                G_ref_f2s(org,[1 4]) = G_ref_f2sT([1 3]);
            else
                G_ref_f2s(org,[2 3]) = G_ref_f2sT([3 1]);
            end
        else
            if(pos==1)
                G_ref_f2s(org,[1 2]) = G_ref_f2sT([1 2]);
            else
                G_ref_f2s(org,[3 4]) = G_ref_f2sT([2 1]);
            end
        end
        G_ref_t(G_ref_f2sT) = 0;
    end
end
