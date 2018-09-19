function REF = mark(xF2S,ind,theta,nu)
% function REF = mark(xF2S,ind,theta,nu)
% xF2S - Father son relation
% ind - error estimator
% theta - refine element? (0..1, 1 = All)
% nu - refine how? (0...1, 0 = Isotrop)
% REF - vector with entries [1 : 4] 

if(size(xF2S,1)==1)
  xF2S = xF2S';
end

T4 = [1 1 1 1; 1 -1 1 -1; 1 1 -1 -1;1 -1 -1 1]/4;

REF=ones(1,size(xF2S,2));
t1=zeros(1,size(xF2S,2));
t3 =0; t4 = 0;

Ct = T4*xF2S;

%% muss ueberhaupt verfeinert werden (welche sollen nicht verfeinert werden)
if(theta <1)
    [s_ind idx] = sort(ind,'descend');

    sum_ind = cumsum(s_ind,1);

    ell = find(sum_ind  >= sum_ind(end) * theta,1);
    
    %Symmetrisieren
    ell = ell + find(abs(( sum_ind(ell)-sum_ind(ell:end)))/sum_ind(ell)>10^-2,1);

    t1(idx(ell+1:end)) = 1;   % Nicht verfeinern
end


%% wie muss verfeinert werden
if(nu > 0) % Horizontal oder Vertikal
    t3 = (nu*abs(Ct(3,:)) >= sqrt(Ct(2,:).^2+Ct(4,:).^2));
    REF(t3-t1==1) = 3;   
    t4 = (nu*abs(Ct(4,:)) >= sqrt(Ct(2,:).^2+Ct(3,:).^2));
    REF(t4-t1==1) = 4;
end
REF(~(t4+t3+t1)) = 2;   % Rest wird horizontal UND vertikal geteilt


end
