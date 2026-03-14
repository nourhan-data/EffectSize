%_________________________________________________________________________________
% Jaya
%
% Modified by Mahamed Omran to work with IEEE CEC 2011, 8-Sep-2020.
%____________________________________________________________________________________

function [best_f, ccurve, dcurve]=HBA_CEC2011(pop_size,maxFEs,lb,ub,dim, fun,opt_f, FuncAddr)

beta = 6;     % the ability of HB to get the food  Eq.(4)
C = 2;     %constant in Eq. (3)
vec_flag=[1,-1];

x = zeros(pop_size,dim);
x_new = x;


% Initialization
x = lb + rand(pop_size, dim).*(ub - lb);

for i=1:pop_size
    f(i)= fEval(x(i,:), fun, opt_f, FuncAddr);
end
f_new = f;

[best_f, gbest] = min(f);
Xprey = x(gbest,:);
    
nfe = pop_size;

ccurve = zeros(1,maxFEs);
ccurve(1:nfe) = min(f);
%diversity
dcurve = zeros(1,maxFEs);
% Determine the initial diversity of the archive    
prum = mean(x);
prum_mat=repmat(prum,pop_size,1);
diam2 = (sum( sum( (x-prum_mat).*(x-prum_mat) ) ) / pop_size);
diversity = sqrt(diam2);
dcurve(1:nfe) = diversity;

while nfe <= maxFEs

    alpha=C*exp(-nfe/maxFEs);   %density factor in Eq. (3)
    I=Intensity(pop_size,Xprey,x); %intensity in Eq. (2)
    
    old_nfe = nfe;
     
    for i=1:pop_size

        r =rand();
        F=vec_flag(floor(2*rand()+1));
        
        for j=1:1:dim
            di=((Xprey(j)-x(i,j)));
            if r<.5
                r3=rand;                r4=rand;                r5=rand;
                
                x_new(i,j)=Xprey(j) +F*beta*I(i)* Xprey(j)+F*r3*alpha*(di)*abs(cos(2*pi*r4)*(1-cos(2*pi*r5)));
            else
                r7=rand;
                x_new(i,j)=Xprey(j)+F*r7*alpha*di;
            end
        end

     %   FU=x_new(i,:)>ub;FL=x_new(i,:)<lb;x_new(i,:)=(x_new(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        x_new(i,:) = max(min(x_new(i,:),ub),lb);
       
        f_new(i)= fEval(x_new(i,:), fun, opt_f, FuncAddr);
    end
    
    nfe = nfe + pop_size;

    for i=1:pop_size    
        if f_new(i) < f(i)
            x(i,:) = x_new(i,:);
            f(i) = f_new(i);
        end
    end
    
    ccurve(old_nfe+1:nfe) = min(f);
    % Determine the diversity of the archive    
    prum = mean(x);
    prum_mat=repmat(prum,pop_size,1);
    diam2 = (sum( sum( (x-prum_mat).*(x-prum_mat) ) ) / pop_size);
    diversity = sqrt(diam2);
    dcurve(old_nfe+1:nfe) = diversity;

    [best_f, gbest] = min(f);
    Xprey = x(gbest,:);

end % while

ccurve(nfe:maxFEs) = best_f;
dcurve(nfe:maxFEs) = diversity;

end % HBA

function I=Intensity(N,Xprey,X)
    for i=1:N-1
        di(i) =( norm((X(i,:)-Xprey+eps))).^2;
        S(i)=( norm((X(i,:)-X(i+1,:)+eps))).^2;
    end
    di(N)=( norm((X(N,:)-Xprey+eps))).^2;
    S(N)=( norm((X(N,:)-X(1,:)+eps))).^2;
    for i=1:N
        r2=rand;
        I(i)=r2*S(i)/(4*pi*di(i));
    end
end

%========================================================================= fEval
  function fx=fEval(x, fun, opt_f,FuncAddr)

      fx = FuncAddr(x, fun, opt_f);
    
    % Unfortunately, there are  sometimes bugs in function evaluations
    % and the position can not be evaluated
    %  (see for example CEC 2011 competition, functions 5 and 6)
    if isnan(fx)
      fx=Inf;
    end
  end
