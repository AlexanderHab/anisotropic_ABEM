function A = build_A2(coordinates,elements)
    N = size(elements,1);

     A1 = zeros(N);
     A1 = mex_build_V(coordinates,elements,1,1);

    % untere schranke s t obere schranke k l
    intF = @(f,a,b,c,d,r,t,u,v)...
        f(b,d,t,v)-f(b,d,t,u)-f(b,d,r,v)+f(b,d,r,u)...
        -f(b,c,t,v)+f(b,c,t,u)+f(b,c,r,v)-f(b,c,r,u)...
        -f(a,d,t,v)+f(a,d,t,u)+f(a,d,r,v)-f(a,d,r,u)...
        +f(a,c,t,v)-f(a,c,t,u)-f(a,c,r,v)+f(a,c,r,u);
    
%     matlabpool
    
    for j = 1:N
        for k = j+1:N
            ej = coordinates(elements(j,[1,2,4])',:);
            ek = coordinates(elements(k,[1,2,4])',:);
            
            d = ej(1,:) - ek(1,:)
            
            ej = ej - repmat(ej(1,:),3,1)
            ek = ek - repmat(ek(1,:),3,1)

%             d = zeros(1,3);
%             if(j~=k)
            A1(j,k) = surfDoubleQuad(@(x1,x2,y1,y2) 1/sqrt((x1-y1-d(1)).^2+(x2-y2-d(2)).^2+d(3).^2)...
                ,ej(1,1),ej(2,1),ej(1,2), ej(3,2),ek(1,1), ek(2,1),ek(1,2), ek(3,2),8);
            A1(k,j) = A1(j,k);
%             else
                
%             end

%             A1(j,k) = intF(@(x1,x2,y1,y2) mex_Fpar(x1,x2,y1,y2,d(1),d(2),d(3))...
%                 ,ej(1,1),ej(2,1),ej(1,2), ej(3,2),ek(1,1), ek(2,1),ek(1,2), ek(3,2));
        end
    end
    
%     matlabpool close
    
    A = A1/(4*pi);

end