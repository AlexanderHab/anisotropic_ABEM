function in = surfDoubleQuad(f,a,b,c,d,r,t,u,v,w)
%
% in = surfDoubleQuad(f,a,b,c,d,r,t,u,v,w)
%
% Berechnet das Vierfachintegral über die Funktion f(x1,x2,y1,y2) mit den 
% Grenzen [a b]x[c d]x[r t]x[u v]. Dazu wird eine Gaussquadratur über w 
% Punkte verwendet.
%
% P. Schaefer

    [x1n x1w] = gauss(w,a,b);
    [x2n x2w] = gauss(w,c,d);
    [y1n y1w] = gauss(w,r,t);
    [y2n y2w] = gauss(w,u,v);

    
    in = 0;
    for i=1:length(y2n)
        for j=1:length(y1n)
            parfor k=1:length(x2n)
                for l = 1:length(x1n)
                    
                    in = in + y2w(i) * y1w(j) * x2w(k) * x1w(l) *...
                      f(x1n(l),x2n(k),y1n(j),y2n(i));
                    
                end
            end
        end
    end

end