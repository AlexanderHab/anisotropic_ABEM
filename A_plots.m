function A_plots(files,printt)

global type2sy
figure(4)
% figure(5)
figure(6)
figure(7)
figure(8)

% type2str = {'Analytisch  '   'Quad Element'   'Quad Achse  '  'Quad Seite  '};
type2str = {'A  ' 'QA'  'QEQA'   'QE  '  '?   '};
type2sy = {'o-' 'v-' 'x-' '.-' '^-','+-','p-','s-','d-'};
type2color = [[0 .3 .3]; [0 .6 .6]; [.6 0 .6 ]; [0 0 .9]];

set(0,'DefaultAxesColorOrder',type2color,...
      'DefaultAxesLineStyleOrder','o-|v-|x-')
    
    set(4,'DefaultAxesColorOrder',type2color(1,:)...
...       ,'DefaultAxesLineStyleOrder','.-|-o|-x|-*|-.|--|:'...
    );
    set(6,'DefaultAxesColorOrder',type2color(1,:)...
...       ,'DefaultAxesLineStyleOrder','.-|-o|-x|-*|-.|--|:'...
    );

G_D=[];
X = [];
leg0 = {};
leg1 = {};
leg2 = {};
leg3 = {};
leg4 = {};
sym = {};

rows = 14;


for i = 1:length(files)
   
   load(files{i});
%    data
      
   [m n] = size(data);
   step = round(n/rows);
   if step<1
    disp (['Error: No Data to show. : ' i])
    continue;
   end
   
   if(size(X,1)>m&&size(X,1)~=0)
       data(size(X,1),end) = 0;
   end
   if(size(X,1)<m&&size(X,1)~=0)
       G_D(m,1) = 0;
       X(m,1) = 0;
   end
   
   G_D = [G_D data(:,2:end)];
   X = [X repmat(data(:,1),1,step)];
   
   p1 = find(files{i}=='/',1,'last')+1;
   if(isempty(p1))
     p1=1;
   end

   p2 = find(files{i}(p1:end)=='_',2)+p1-1;
   if(length(p2)<2)
      p2(2) = length(files{i});
   else
     p2(2) = p2(2) -1;
   end
   
   t0 = ['(' files{i}(p2(1)+1:p2(2)) ')'];
   l0 = [files{i}(p1:p2(1)-1) ' '];
   l1 = {type2str{data(1,[2+(0:step-1)*rows])}};
%    l00 = l0;
   for i = 1:step
%       l0 = [l00 num2str(i-1) ' '];
      leg0 = {leg0{:}...
          ['tmu ' l0 l1{i}]...
...          ['tilde \mu_2 ' l0 l1{i}]...
          ['eta ' l0 l1{i}]...
          ['fehler ' l0 l1{i}]...
...          ['\mu ' l0 l1{i}]...
...          ['\kappa ' l0 l1{i}]...
...          ['\kappa2 ' l0 l1{i}]...
...          ['\kappa3 ' l0 l1{i}]...
          }';
      leg1 = {leg1{:}...
          [ l0 l1{i}]...
          }';
      leg2 = {leg2{:}...
          ['min hmin/max hmax ' l0 l1{i}]...
          ['min hmax/max hmax ' l0 l1{i}]...
          ['min hmin/hmax ' l0 l1{i}]...
          }';
      
      leg3 = {leg3{:}...
          ['cond ' l0 l1{i}]...
          }';
      leg4 = {leg4{:}...
          ['Zeit aufbauen ' l0 l1{i}]...
          ['Zeit loesen ' l0 l1{i}]...
          }';
      sym = {sym{:} type2sym(data(1,[2+(i-1)*rows]))}';
   end
   
end
  
%   [u ia ic] = unique(leg0)

   [m n] = size(G_D);
   step = round(n/rows);

if step<1
    disp ('Error: No Data to show.')
else

%     data((end-8):end,[1 [3 4 5]])
    
% sol = interp1(1./X((round(1)):(end),1)',G_D((round(1)):(end),4)',0,'spline')

% G_D(end,4)

  sol = 8.30235;
   
   
if(strcmp(t0, '(3DFichCube)') || strcmp(t0, '(3DFichCube2)'))
  sol = 16.2265;  % Ferraz-Leite Paper
%   sol = 50;
elseif(strcmp(t0, '(2DQuad)'))
  sol = 4.609193;
elseif(strcmp(t0, '(2DLShape)'))
%   sol = 8.28457; % Ferraz-Leite Dipl.
  sol = 8.28466; % Ferraz-Leite Paper
elseif(strcmp(t0, '(3DCube)') || strcmp(t0, '(3DCube2)'))
  sol = 8.30235;
elseif(strcmp(t0, '(3DCube3)'))
%   sol = 1.0379;
  sol = 16.604703 % Ferraz-Leite Dipl.
end

% G_D(:,2+2)
% 
% abs(log10((sqrt(sol-G_D(:,2+2))-G_D(:,2+1))))




% G_D

%% Plotte Fehlerschätzer
figure(4)
i=0;
% [shift k] = min(G_D(1:end-5,2+1+rows*i)*G_D(1,2)/G_D(1,3)-G_D(1:end-5,2+0+rows*i));
% shift = shift+shift/10;

% eta = G_D(:,2+1+rows*i)*(G_D(k,2)-shift)/G_D(k,3)
% 
% error = sqrt(abs(sol - G_D(:,2+2+rows*i)))
% [shift2 l] = min(error*eta(1)/error(1) - eta);
% shift2 = shift2+shift2/10;
% error*(eta(l)-shift2)/error(l)

 first = 0 +find(([ (G_D(2:end,2+rows*i)-G_D(1:end-1,2+rows*i))./G_D(2:end,2+rows*i)])>=0,1)
 Ferr = (X(first,i+1)+X(first-1,i+1))/2;

loglog(X(:,i+1),G_D(:,2+rows*i),type2sym(i*3+1), ...
...     G_D(:,2+8+rows*i)...
    X(:,i+1),G_D(:,2+1+rows*i),type2sym(i*3+2),...*(G_D(k,2)-shift)/G_D(k,3)...*G_D(1,2)/G_D(1,2+1+rows*i) ...
    X(:,i+1),sqrt(abs(sol - G_D(:,2+2+rows*i))),type2sym(i*3+3),...*(G_D(k,2)-shift)/G_D(k,3)...
...    G_D(:,2+3+rows*i)...*G_D(1,2)/G_D(1,2+3+rows*i) ...
...    [ 0; sqrt(G_D(2:end,2+9+rows*i)-G_D(1:end-1,2+9+rows*i))]...
...    [ 0; sqrt(G_D(2:end,2+2+rows*i)-G_D(1:end-1,2+2+rows*i))]...
...    G_D(:,2+10+rows*i)...
    'color', type2color(i+1,:));
hold on

for i = 1:step-1
loglog(X(:,i+1),G_D(:,2+rows*i),type2sym(i*3+1), ...
...     G_D(:,2+8+rows*i)...
    X(:,i+1),G_D(:,2+1+rows*i),type2sym(i*3+2),...*(G_D(k,2)-shift)/G_D(k,3)...*G_D(1,2)/G_D(1,2+1+rows*i) ...
    X(:,i+1),sqrt(abs(sol - G_D(:,2+2+rows*i))),type2sym(i*3+3),...*(G_D(k,2)-shift)/G_D(k,3)...
...    G_D(:,2+3+rows*i)...*G_D(1,2)/G_D(1,2+3+rows*i) ...
...    [ 0; sqrt(G_D(2:end,2+9+rows*i)-G_D(1:end-1,2+9+rows*i))]...
...    [ 0; sqrt(G_D(2:end,2+2+rows*i)-G_D(1:end-1,2+2+rows*i))]...
...    G_D(:,2+10+rows*i)...
    'color', type2color(i+1,:));
    
end

loglog(X(:,1),9*X(:,1).^(-1/2),'-.',...
  X(:,1),3*X(:,1).^(-1/4),':',...
  X(:,1),2*X(:,1).^(-3/4),'--','color', [.9 .6 0])

if(~isempty(Ferr))
loglog([Ferr; Ferr],...
  [min(G_D(:,3+rows*[0:step-1]))*0.6; max(G_D(:,2+rows*[0:step-1]))*1.3],...
  '--','color',[.9 0 0]);
end

hold off

% title(['Fehler ' t0])
xlabel('Elementanzahl');
ylabel('Schaetzer');
legend({leg0{:} ...
      'N-12'  'N-14'  'N-34'...
      'erste Fehler' ...
     } ,'location','southwest','box', 'off');

%print('-r600','-depsc',[printt '_error.eps'])
 
% %% Plotte eNorm
% figure(5)
% loglog(repmat(X(:,1),1,1),G_D(:,2+2),type2sym{1});
% hold on
% for i = 1:step-1
%     loglog(repmat(X(:,i+1),1,1),G_D(:,2+2+i*rows),'color', type2color(i+1,:));
% end
% loglog(X(:,1),repmat(sol,size(X,1),1),'-.', 'color', [.9 0 0])
% if(~isempty(Ferr))
% loglog([Ferr; Ferr],...
%   [min(G_D(:,4+rows*[0:step-1]))*0.6; max(G_D(:,4+rows*[0:step-1]))*1.1],...
%   '--','color',[.9 0 0]);
% end
% hold off
% yl = ylim();
% ylim([yl(1) 1.005*yl(2)])
% % title(['Energie Norm ' t0])
% xlabel('Elementanzahl');
% ylabel('eNorm^2');
% legend({leg1{:} , 'extrapoliert', 'erste Fehler' },'location','SouthEast');
% 
% 
% print('-r600','-depsc',[printt '_norm.eps'])

%% Plotte HMIN HMAX
figure(6)
i=0;
loglog(...
    X(:,i+1),G_D(:,2+4+rows*i),type2sym(i*3+1),...
    X(:,i+1),G_D(:,2+5+rows*i),type2sym(i*3+2),...
    X(:,i+1),G_D(:,2+6+rows*i),type2sym(i*3+3),...
    'color', type2color(i+1,:));
hold on
for i = 1:step-1
loglog(...
    X(:,i+1),G_D(:,2+4+rows*i),type2sym(i*3+1),...
    X(:,i+1),G_D(:,2+5+rows*i),type2sym(i*3+2),...
    X(:,i+1),G_D(:,2+6+rows*i),type2sym(i*3+3),...
    'color', type2color(i+1,:));
end
% loglog(X(:,1),[7*X(:,1).^(-1/2),3*X(:,1).^(-1/4),2*X(:,1).^(-3/4)],'-.')

if(~isempty(Ferr))
loglog([Ferr; Ferr],...
  [min(G_D(:,6+rows*[0:step-1]))*0.6; max(G_D(:,6+rows*[0:step-1]))*1.3],...
  '--','color',[.9 0 0]);
end
yl = ylim();
ylim([yl(1) 2])


hold off

% title(['hmin hmax ' t0])
xlabel('Elementanzahl');
ylabel('Verhaeltnis');
legend({leg2{:} ...
      'erste Fehler' ...
%       'N^{-1/2}'  'N^{-1/4}'  'N^{-3/4}'...
     } ,'location','southwest','box', 'off');

%print('-r600','-depsc',[printt '_hminmax.eps'])


%% Plotte Kond
figure(7)
i=0;
loglog(repmat(X(:,i+1),1,1),...
    G_D(:,2+7+rows*i)...
    ,type2sym(i+1),'color', type2color(i+1,:));
hold on
for i = 1:step-1
loglog(repmat(X(:,i+1),1,1),...
    G_D(:,2+7+rows*i)...
    ,type2sym(i+1),'color', type2color(i+1,:));
end
% loglog(X(:,1),[7*X(:,1).^(-1/2),3*X(:,1).^(-1/4),2*X(:,1).^(-3/4)],'-.')

if(~isempty(Ferr))
loglog([Ferr; Ferr],...
  [min(G_D(:,9+rows*[0:step-1]))*0.6; max(G_D(:,9+rows*[0:step-1]))*1.3],...
  '--','color',[.9 0 0]);
end

hold off

% title(['KonditionierungsZahlen von V ' t0])
xlabel('Elementanzahl');
ylabel('Kondition');
legend({leg3{:} ...
      'erste Fehler' ...
%       'N^{-1/2}'  'N^{-1/4}'  'N^{-3/4}'...
     } ,'location','northwest','box', 'off');

%print('-r600','-depsc',[printt '_cond.eps'])

%% Plotte Zeit
figure(8)
i=0;
loglog(...
    X(:,i+1),G_D(:,2+11+rows*i),type2sym(i*3+1),...  
    X(:,i+1),G_D(:,2+12+rows*i),type2sym(i*3+2),...
      'color', type2color(i+1,:));
hold on
for i = 1:step-1
loglog(...
    X(:,i+1),G_D(:,2+11+rows*i),type2sym(i*3+1),...  
    X(:,i+1),G_D(:,2+12+rows*i),type2sym(i*3+2),...
      'color', type2color(i+1,:));
end
% loglog(X(:,1),[7*X(:,1).^(-1/2),3*X(:,1).^(-1/4),2*X(:,1).^(-3/4)],'-.')
hold off

% title(['Benoetigte Zeit ' t0])
xlabel('Elementanzahl');
ylabel('Sekunden');
legend({leg4{:} ...
%       'N^{-1/2}'  'N^{-1/4}'  'N^{-3/4}'...
     } ,'location','northwest','box', 'off');

%print('-r600','-depsc',[printt '_time.eps'])

end

end


function str = type2sym(num)
global type2sy
  str = type2sy{mod(num-1,length(type2sy))+1};
end