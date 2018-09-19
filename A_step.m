function [data time er] = A_step(zeta,type,theta,nu,vcon,varargin)
% [data time er] = A_step(zeta,type,theta,nu,vcon [, out])
% Führt einen Verfeinerungsschritt aus
%
% zeta & type - Bestimmen die Art der Matrix Berechnung
% theta - adaptiv?
% nu - isotrop?
% vcon - Vorkonditionierung der Matrix? 1 oder 0
% out - (optional) Dateizusatz um Netz, Lösung & co zu speichern
%
% P. Schaefer

%globale Speicherung initialisieren
global G_E;
global G_C;
global G_N;
global G_T;
global G_D;
global G_S;

%Netz darf nicht leer sein
assert(~isempty(G_E) && ~isempty(G_C) && ~isempty(G_N),...
  'Please use A_loadMesh first')

%Steuerungs Parameter prüfen
% assert(length(zeta)==max(type)-1||length(zeta)==1,...
%   'Pleas set right type and zeta parameters')

%Steuerungsparameter korrigieren
% if(length(zeta)==1)
%   zeta = repmat(zeta,max(type)-1,1);
% end

%Vektor für Zeitmessung anlegen
time = zeros(1,3);
  tic
  %uniformIsotrop Verfeinern
  [coo_fine,ele_fine,neigh_fine,f2s,sit_fine]=refineQuad(G_C,G_E,G_N,G_S,2);
  time(1) = toc;

  %Flaecheninhalte Berechnen (rhs)
  b_fine = areaQuad(sit_fine);b_fine = areaQuad(sit_fine);
  
  b = areaQuad(G_S);
  hmin = 2.^-max(G_S,[],2);
  hmax = 2.^-min(G_S,[],2);

  tic
  %data -> ErgebnisMatrix
  data = size(G_E,1);
  %save_* -> ZwischenSpeicherung
%  save_A = {};
%  save_x = {};
%  save_A_fine = {};
%  save_x_fine = {};
  
  %Alle MatrixBrechenungsArten mit dem selben Netz berechnen
  for i = 1:length(type)
    %Matrix aufbauen -> MEX
    V_fine = mex_build_V(coo_fine,ele_fine,zeta,type(i));
    
    %Testet auf Fehlerhafte Einträge (NaN +/-Inf)
    [r c] = find(isnan(V_fine)~=isinf(V_fine));
    if(~isempty(r))
        figure(9)
        plotShape(coo_fine,ele_fine(unique([r c]),:),'')
        plotMark([r';c'],coo_fine,ele_fine)
        title('Fehlerhafte Elemente')
    end
 
    if(~vcon)
    %Lösung Berechnen
        x_fine = V_fine\b_fine;
        con = cond(V_fine);
    else
    %Vorkonditionierte Lösung!
        D = diag(diag(V_fine.^(-1/2)));

        A = D * V_fine * D;
        c = D*b_fine;
        y = A\c;
        x_fine = D*y;
        con = cond(A);
    end
    
    % \tilde \mu ( \Pi h -h + L_2 )
    tmu = hmin.* b .* sum((x_fine(f2s)'-repmat(sum(x_fine(f2s)',1)/4,4,1)).^2)' /4;
    
    %Fehlerschätzer 2 aufbauen
    V = mex_build_V(G_C,G_E,zeta,type(i));
    
    if(~vcon)
        x = V\b;
    else
        D = diag(diag(V.^(-1/2)));

        A = D * V * D;
        c = D*b;
        y = A\c;
        x = D*y;
    end
    
    xo_fine(f2s) = repmat(x,1,4);
    xd_fine = xo_fine'-x_fine;
    
    % |||h/2 -h|||
%     eta = xd_fine'*A_fine*xd_fine;
    
    % \tilde \mu ( h/2 -h + L_2 )
    mu = hmin.*b.*sum((x_fine(f2s)'-repmat(x',4,1)).^2)'/4;


    %Energienorm^2 Berechnen |||h||| & |||h/2|||
%     xe_fine = x_fine'*A_fine*x_fine;
    xe_fine = b_fine'*x_fine;
%     xe = x'*A*x;
    xe = b'*x;
    
    %\tilde \mu 2 = ( ||\Pi h|| - ||h||)
    tmu2 = hmin.* b.* (...
    sum((x_fine(f2s)').^2)'-sum(repmat(sum(x_fine(f2s)',1)/4,4,1).^2)'...
    ) /4;
    
    eta = abs(xe_fine-xe);
    
%    save_A_fine{i} = V_fine;
%    save_x_fine{i} = x_fine;
%    
%    save_A{i} = V;
%    save_x{i} = x;
   
    
    
    
    data = [data ...
        type(i) ...
        sqrt(sum(tmu))...
        sqrt(eta) ...
        xe ...
        sqrt(sum(mu))...
        min(hmin)/max(hmax)...
        min(hmax)/max(hmax)...
        min(hmin./hmax) con...
        sqrt(sum(tmu2))...
        xe_fine...
        ];
  end
  time(2) = toc;
  
  % nur RandElemente Verfeinern
%   marked = ones(1,size(G_E,1));
%   marked(find(sum((G_N(:,1:4)==0),2))) = 2;

  % Markieren mit gewählten Parametern
  marked = mark(x_fine(f2s)',tmu,theta,nu);
  
  % Bunt Plotten!
%   figure(1)
%   plotShape(G_C,G_E,'s',tmu);
%   title('Elemente mit Fehlerschaetzer')
%   colorbar
%   view(2)
%   plotMark(find(marked>1),G_C,G_E);
  
  assert(size(G_E,1)==length(marked),'MarkierungsVektor ist fehlerhaft')
  
  tic
  er = [];

  
  old_C = G_C;
  old_E = G_E;
  old_S = G_S;
  
  %Netz Verfeinern, wie durch marked bestimmt
  [G_C, G_E, G_N, f, G_S, er] = refineQuad(G_C,G_E,G_N,G_S,marked);
  
  %Vater Sohn test
  assert(sum(areaQuad(old_S))==sum(areaQuad(G_S)),'Gesamtinhalt Fehlerhaft')
  
  
  time(3) = toc;
%   if(~isempty(er))
%     figure(10)
%     plotShape(G_C,G_E,'');
%     plotMark(er,G_C,G_E);
%   end
 
  %ErgebnisWerte und Zeiten Speichern
  G_T(size(G_T,1)+1,1:4) = [size(G_E,1) time];
  G_D(size(G_D,1)+1,1:length(data)) = data;

  
  %Alle Relevanten zwischenInformationen Speichern
%   out = '_';
%   if(length(varargin)~=0)
%      out = varargin{1};
%    end
%   typeN = int2str(type);
%  save (['meshSave/fine_' out typeN(typeN~=' ') '_' int2str(size(G_T,1))],...
%      'coo_fine', 'ele_fine','neigh_fine','f2s','data',...
%      'save_A','save_x','save_A_fine','save_x_fine','b','b_fine')

%   clear 'coo_fine' 'ele_fine' 'neigh_fine' 'f2s' 
%   plotShape(G_C,G_E,'');

end

