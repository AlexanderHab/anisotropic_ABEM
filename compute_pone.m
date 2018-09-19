function [data, er, fileo] = compute_pone(file,times,zeta,typ,theta,nu,vcon)
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
% /*************************************************************************/
% /* Code is based on the bachelor thesis of Peter Schäfer                 */
% /* adaptiert von Alexx 2018                                              */
% /*************************************************************************/

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
  % do some step of uniform refinement
  num_steps = 1; 
  
  coo_fine = coordinates;
  ele_fine = elements; 
  neigh_fine = neigh; 
  sit_fine = sites; 
  
  for ell=1:num_steps 
      [coo_fine,ele_fine,neigh_fine,f2s,sit_fine]...
        =refineQuad(coo_fine,ele_fine,neigh_fine,sit_fine,2);
  end


 fprintf("refine ok \n")

  %Flaecheninhalte berechnen (rhs)
  b_fine = 2.^-sum(sit_fine,2);
  b = 2.^-sum(sites,2);
  
  
  el_size = 1/4;
  el_size_fine = 1/4^(num_steps+1);
  
  b_pone = repmat([el_size;el_size*1/4;el_size*1/4;el_size*1/4*1/4],24,1);
  %b_pone = repmat([1/4*0.25;1/4*1/4;1/4*1/4;1/4*1/4*1/4],24,1);
  b_pone_fine = repmat([el_size_fine;el_size_fine*1/4;el_size_fine*1/4;el_size_fine*1/4*1/4],24*4^num_steps,1);
  
  
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
    
    fprintf("jump into the mex file \n");
    
    
    V_pone = mex_build_V_pone(coordinates,elements,zeta{i},typ(i));
    build_time = toc;
    %testet auf fehlerhafte Eintraege (NaN +/-Inf)
    [r c] = find(isnan(V_pone)~=isinf(V_pone));
    
    
    V = mex_build_V(coordinates,elements,zeta{i},typ(i));
    build_time = toc;
    %testet auf fehlerhafte Eintraege (NaN +/-Inf)
    [r c] = find(isnan(V)~=isinf(V));
    
    
    
    V_fine_pone = mex_build_V_pone(coo_fine,ele_fine,zeta{i},typ(i));
    build_time = toc;
    [r c] = find(isnan(V_fine_pone)~=isinf(V_fine_pone));
    
    
    V_fine = mex_build_V(coo_fine,ele_fine,zeta{i},typ(i));
    build_time = toc;
    %testet auf fehlerhafte Eintraege (NaN +/-Inf)
    [r c] = find(isnan(V_fine)~=isinf(V_fine));

 
    
     fprintf("precon \n")
    if(~vcon)
    %Loesung berechnen
       x_fine = V_fine\b_fine;
       x_fine_pone = V_fine_pone\b_pone_fine;
       
       x = V\b;
       x_pone = V_pone\b_pone;
        
    solve_time = toc;
    
    %Konditionszahl aufstellen
    con = cond(V);
    con_pone = cond(V_pone);
   
    disp("condition number P1");
    disp(con_pone);
    disp("difference in condition number P0");
    disp(con);
    
    
  
%     clear V_fine
    
    % \tilde \mu ( \Rho h -h + L_2 )
%     tmu = hmin.* b .* ...
%       sum((x_fine(f2s)'-repmat(sum(x_fine(f2s)',1)/4,4,1)).^2)' /4;
    
%     %Fehlerschaetzer 2 aufbauen
%     V = mex_build_V(coordinates,elements,zeta{i},typ(i));
%     
%     if(~vcon)
%         x = V\b;
%     else
%         D = diag(V).^(-1/2);
%         for k = 1:length(V)
%           for l = 1:length(V)
%             V(k,l) = V(k,l)*D(k)*D(l);
%           end
%         end
%         c = D.*b;
%         y = V\c;
%         x = D.*y;
%     end
%     
% %     clear V
%     
%     xo_fine(f2s) = repmat(x,1,4);
%     xd_fine = abs(xo_fine'-x_fine);
    
    % \mu ( h/2 -h + L_2 )
%     mu = hmin.*b.*sum((x_fine(f2s)'-repmat(x',4,1)).^2)'/4;


    %Energienorm^2 berechnen |||h||| & |||h/2|||
    xe_fine = b_fine'*x_fine;
    xe_fine_pone = b_pone_fine'*x_fine_pone;
    
    xe = b'*x;
    xe_pone = b_pone'*x_pone;
    
    
    disp("difference in energy norm coarse ");
    disp(abs(xe_pone-xe));
    disp("energynorm P0 coarse :");
    disp(xe);
    disp("energynorm P1 coarse :");
    disp(xe_pone);
    
    
%     
    disp("difference in energy norm fine ");
    disp(abs(xe_fine_pone-xe_fine));
    disp("energynorm P0 fine :");
    disp(xe_fine);
    disp("energynorm P1 fine :");
    disp(xe_fine_pone);
    
    
    
    disp("end");
    
%     xe = x'*A*x;
%    xe = b'*x;
    
    
%     %Enorm^2 Elementweise vergleichen
%     if(exist('old_F_fine','var'))
%       for k = 1:size(ele_fine,1)
%         e_f = find(sum(k==f2s,2));
%         e_f_p = find(f2s(e_f,:)==k);
%         assert(size(e_f_p,1) ==1,'Error: e_f_p wurde mehrfach gefunden')
%         e_ff = find(sum(e_f==f,2));
%         e_ff_p = find(f(e_ff,:)==e_f);
%         assert(size(e_ff_p,1)==1,'Error: e_ff_p wurde mehrfach gefunden')
%         e_ff_s = length(e_ff_p);
%         if e_ff_s==4
%            e_p = e_f_p;
%         elseif e_ff_s == 1
%            e_p = e_ff_p;  
%         elseif e_ff_s == 2      
%            if(sum(e_ff_p) == 5)
%              e_p = e_ff_p(floor((e_f_p-1)/2)+1);
%            else
%              inh = [1 2 4 3];
%              e_p = inh(e_ff_p(mod(floor((e_f_p)/2),2)+1));
%            end
%         end        
%         xf_fine(k,1) = abs(x_fine(k) - old_x_fine(old_F_fine(e_ff,e_p)));
%       end
% 
%       kap3 = b_fine'*(xf_fine);
%     end

% %     
% %     %\tilde \mu 2 = ( ||\Pi h|| - ||h||)
% %     tmu2 = hmin.* b.* (sum((x_fine(f2s)').^2)'...
% %       -sum(repmat(sum(x_fine(f2s)',1)/4,4,1).^2)') /4;
% %   
% %     % |||h/2 -h|||
% % %     eta = xd_fine'*A_fine*xd_fine;    
% %     eta = abs(xe_fine-xe);
% %      
% %     end_time = toc - start_time;
% % %     disp(['Relative Zeit ' num2str(100*(build_time - start_time)/end_time)...
% % %       '% absolute Zeit '  t2str(build_time - start_time)]);
% %     time_build = build_time - start_time;
% %     time_solve = solve_time - build_time;
% %     


  
end
  end
end



