function [data, er, fileo] = compute(file,times,zeta,typ,theta,nu,vcon)
%***************************************************************************
%* [data er fileo] = compute(file,times,zeta,typ,theta,nu,vcon)            *
%*                                                                         *
%*   Fuehrt times Verfeinerungsschritte mit gewaehlten Parametern aus.     *
%*                                                                         *
%*   file - StartNetz                                                      *
%*   times - wie oft Verfeinert werden soll                                *
%*   zeta & typ - Bestimmen die Art der Matrix Berechnung                  *
%*   theta - adaptiv?                                                      *
%*   nu - isotrop?                                                         *
%*   vcon - Vorkonditionierung der Matrix? 1 oder 0                        *
%*                                                                         *
%***************************************************************************
%* Author: Peter Schaefer                             schaeferpm@gmail.com *
%* Version: 1.0  (2012)                                                    *
%***************************************************************************

% Datei laden
load(file)

% Container fuer Messungen initialisieren
if(~exist('data','var'))
    data=[];
end

kap3 = 0;

tic

%times mal verfeinern oder bis Elementanzahl fuer times > 40
for j = 1:times
    
  disp(['number of elements'  size(elements,1)])
  %beende Schleife wenn Elementanzahl erreicht
  if(times>40 && size(elements,1) > times)
    break;
  end
  
  fprintf("test 1 \n")
  
  
  %altes Netz mit Loesung merken
  if(exist('coo_fine','var'))
    old_C_fine = coordinates;

    old_E_fine = elements;

    
  else
    old_C_fine = coordinates;
    old_E_fine = elements;
  end

  %uniformIsotrop verfeinern
  [coo_fine,ele_fine,neigh_fine,f2s,sit_fine]...
    =refineQuad(coordinates,elements,neigh,sites,2);

 fprintf("refine ok \n")

  %Flaecheninhalte berechnen (rhs)
  b_fine = 2.^-sum(sit_fine,2);
  
  b = 2.^-sum(sites,2);
  hmin = 2.^-max(sites,[],2);
  hmax = 2.^-min(sites,[],2);


  %data -> ErgebnisVektor
  dataS = size(elements,1);
  
  %zeta vorbereiten
  if(~iscell(zeta))
    zeta_tmp = zeta;
    clear zeta;
    for i = 1:length(typ)
      zeta{i} = zeta_tmp;
    end
  end
  
  fprintf("prep ok compute matrix \n")
  %alle MatrixBrechenungsArten mit dem selben Netz berechnen
  
  for i = 1:length(typ)
    start_time = toc;
    disp([ '[' num2str(j) ',' num2str(i) '] ' ...
      num2str(size(elements,1)) ' : ' t2str(toc) ' ->' num2str(typ(i))])
    %Matrix aufbauen -> MEX
    fprintf("jump into the mex file \n")
    
    V_fine = mex_build_V(coo_fine,ele_fine,zeta{i},typ(i));
    build_time = toc;
    %testet auf fehlerhafte Eintraege (NaN +/-Inf)
    fprintf("finished max \n")
    [r c] = find(isnan(V_fine)~=isinf(V_fine));
    if(~isempty(r))
        figure(9)
        plotShape(coo_fine,ele_fine(unique([r c]),:),'')
        plotMark([r';c'],coo_fine,ele_fine)
        title('Fehlerhafte Elemente')
    end
 
    
     fprintf("precon \n")
    if(~vcon)
    %Loesung berechnen
        x_fine = V_fine\b_fine;
        
    else
    %Vorkonditionierte Loesung!
        D = diag(V_fine).^(-1/2);
        for k = 1:length(V_fine)
          for l = 1:length(V_fine)
            V_fine(k,l) = V_fine(k,l)*D(k)*D(l);
          end
        end
        c = D.*b_fine;
        y = V_fine\c;
        x_fine = D.*y;
    end
    
    solve_time = toc;
    
    %Konditionszahl aufstellen
    con = cond(V_fine);
    
    clear V_fine
    
    % \tilde \mu ( \Rho h -h + L_2 )
    tmu = hmin.* b .* ...
      sum((x_fine(f2s)'-repmat(sum(x_fine(f2s)',1)/4,4,1)).^2)' /4;
    
    %Fehlerschaetzer 2 aufbauen
    V = mex_build_V(coordinates,elements,zeta{i},typ(i));
    
    if(~vcon)
        x = V\b;
    else
        D = diag(V).^(-1/2);
        for k = 1:length(V)
          for l = 1:length(V)
            V(k,l) = V(k,l)*D(k)*D(l);
          end
        end
        c = D.*b;
        y = V\c;
        x = D.*y;
    end
    
    clear V
    
    xo_fine(f2s) = repmat(x,1,4);
    xd_fine = abs(xo_fine'-x_fine);
    
    % \mu ( h/2 -h + L_2 )
    mu = hmin.*b.*sum((x_fine(f2s)'-repmat(x',4,1)).^2)'/4;


    %Energienorm^2 berechnen |||h||| & |||h/2|||
%     xe_fine = x_fine'*A_fine*x_fine;
    xe_fine = b_fine'*x_fine;
%     xe = x'*A*x;
    xe = b'*x;
    
    
    %Enorm^2 Elementweise vergleichen
    if(exist('old_F_fine','var'))
      for k = 1:size(ele_fine,1)
        e_f = find(sum(k==f2s,2));
        e_f_p = find(f2s(e_f,:)==k);
        assert(size(e_f_p,1) ==1,'Error: e_f_p wurde mehrfach gefunden')
        e_ff = find(sum(e_f==f,2));
        e_ff_p = find(f(e_ff,:)==e_f);
        assert(size(e_ff_p,1)==1,'Error: e_ff_p wurde mehrfach gefunden')
        e_ff_s = length(e_ff_p);
        if e_ff_s==4
           e_p = e_f_p;
        elseif e_ff_s == 1
           e_p = e_ff_p;  
        elseif e_ff_s == 2      
           if(sum(e_ff_p) == 5)
             e_p = e_ff_p(floor((e_f_p-1)/2)+1);
           else
             inh = [1 2 4 3];
             e_p = inh(e_ff_p(mod(floor((e_f_p)/2),2)+1));
           end
        end        
        xf_fine(k,1) = abs(x_fine(k) - old_x_fine(old_F_fine(e_ff,e_p)));
      end

      kap3 = b_fine'*(xf_fine);
    end

    
    %\tilde \mu 2 = ( ||\Pi h|| - ||h||)
    tmu2 = hmin.* b.* (sum((x_fine(f2s)').^2)'...
      -sum(repmat(sum(x_fine(f2s)',1)/4,4,1).^2)') /4;
  
    % |||h/2 -h|||
%     eta = xd_fine'*A_fine*xd_fine;    
    eta = abs(xe_fine-xe);
     
    end_time = toc - start_time;
%     disp(['Relative Zeit ' num2str(100*(build_time - start_time)/end_time)...
%       '% absolute Zeit '  t2str(build_time - start_time)]);
    time_build = build_time - start_time;
    time_solve = solve_time - build_time;
    
    dataS = [dataS ...
        typ(i) ... berechnet mit Typ
        sqrt(sum(tmu))... tilde mu
        sqrt(eta) ... eta
        xe ... error (kappa 2)
        sqrt(sum(mu))... mu
        min(hmin)/max(hmax)...
        min(hmax)/max(hmax)...
        min(hmin./hmax)...
        con... Kondition
        sqrt(sum(tmu2))... tilde mu 2
        xe_fine... (kappa)
        kap3 ... kappa 3
        time_build ... benoetigte Zeit (Aufbauen)
        time_solve ... benoetigte Zeit (Berechnen)
        ];
  end

  % markieren mit gewaehlten Parametern
  marked = mark(x_fine(f2s)',tmu,theta,nu);
  
  % Netz bunt plotten!
   figure(1)
   plotShape(old_C_fine,old_E_fine,'s',tmu);
   title('Elemente mit Fehlerschaetzer')
%   colorbar
%   view(2)
%   plotMark(find(marked>1),G_C,G_E);
  
  assert(size(elements,1)==length(marked),...
    'MarkierungsVektor ist fehlerhaft')
  
  old_C = coordinates;
  old_E = elements;
  old_S = sites;
  
  %Netz verfeinern, wie durch marked bestimmt
  [coordinates, elements, neigh, f, sites, er]...
    = refineQuad(coordinates,elements,neigh,sites,marked);
  
  %Vater Sohn test
  assert(sum(2.^-sum(old_S,2))==sum(2.^-sum(sites,2)),...
    'Gesamtinhalt Fehlerhaft')
  
  e = 0;
  for k = 1:size(data,1)
    e = e +  abs(2.^-sum(old_S(k,:))...
      - sum(2.^-sum(sites(unique(f(k,:)'),:),2)));
  end
  assert(e==0, 'Einzelne Inhalte Fehlerhaft')
  
  % Fehlerhafte Elemente Plotten
%   if(~isempty(er))
%     figure(10)
%     plotShape(G_C,G_E,'');
%     plotMark(er,G_C,G_E);
%   end
   
   p1 = find(file=='/',1,'last')+1;
   if(isempty(p1))
     p1=1;
   end

   p2 = find(file(p1:end)=='_',2)+p1;
   if(length(p2)<2)
      p2(2) = length(file);
   else
     p2(2) = p2(2) -2;
   end
  
  %ErgebnisWerte speichern
  data(size(data,1)+1,1:length(dataS)) = dataS;
  typeN = int2str(typ);
%   fileo = [typeN(typeN~=' ')...
%     't' regexprep(num2str(theta,2),'\.','')...
%     'n' regexprep(num2str(nu,2),'\.','') '_'...
%     file(p2(1):p2(2)) '_'];
fileo = 'test'
save (['meshSave/' fileo int2str(size(data,1))]...
    , 'coordinates', 'elements','neigh','data', 'sites')
  
end
end


function str = t2str(time)
  type = 's';
  
  if(time/60>1)
    time =  time/60;
    type = 'm';
    
    if(time/60>1)
      time =  time/60;
      type = 'h';  
    end 
  end

str = [num2str(round(time*1000)/1000) type];
end
