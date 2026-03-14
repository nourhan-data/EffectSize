%_________________________________________________________________________________
% Jaya
%
% Modified by Mahamed Omran to work with IEEE CEC 2011, 8-Sep-2020.
%____________________________________________________________________________________

function [best_f, ccurve, dcurve]=EnSO_CEC2011(pop_size,maxFEs,lb,ub,dim, fun,opt_f, FuncAddr, b)

if nargin < 9
    b = 2;
end

x = zeros(pop_size,dim);
x_new = x;

f = zeros(1,pop_size);
f_new = f;

% Initialization
x = lb + rand(pop_size, dim).*(ub - lb);

for i=1:pop_size
    f(i)= fEval(x(i,:), fun, opt_f, FuncAddr);
end

f_new = f;
    
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

C = 1; % as suggested in the original paper of BFPA

hpop_size = fix(pop_size/3);

W = 5;

while nfe <= maxFEs
    
    old_nfe = nfe;
    
   % Global knowlege base
   pop1 = x(1:hpop_size,:);
   fpop1 = f(1:hpop_size);
   [~, indices] = sort(fpop1);
   pop1 = pop1(indices,:);
   
   pop2 = x(hpop_size+1:2*hpop_size,:);
   fpop2 = f(hpop_size+1:2*hpop_size);
   [~, indices] = sort(fpop2);
   pop2 = pop2(indices,:);
   
   pop3 = x(2*hpop_size+1:end,:);
   fpop3 = f(2*hpop_size+1:end);
   [~, indices] = sort(fpop3);
   pop3 = pop3(indices,:);
   
   elite = [pop1(1:b,:); pop2(1:b,:); pop3(1:b,:)];
   best1 = elite(randi(3*b),:);
   best2 = elite(randi(3*b),:);
   best3 = elite(randi(3*b),:);
   worst = pop3(end,:);
   
   RimeFactor = (rand - 0.5)*2*cos((pi*nfe/(maxFEs/10)))*(1 - round(nfe*W/maxFEs)/W);
   E = sqrt(nfe/maxFEs);
   normalized_rime_rates = normr(fpop2);
      
    for i=1:pop_size
        
        if i <= hpop_size % BFPA
            r = 2*C*rand - C;
            M = 1:hpop_size;
            
            M(i) = [];  % remove i from {1,..,pop_size}
            p = M(randi(hpop_size-1));
            q = M(randi(hpop_size-1));
            x_new(i,:) = x(i,:) + r*(best2 - x(q,:)); %r*(x(p,:) - x(q,:));
        elseif i <= 2*hpop_size  % RIME optimization algorithm
            
            for j=1:dim
                r1 = rand();
                if r1 < E
                    x_new(i,j) = best3(j) + RimeFactor*((ub(j) - lb(j))*rand + lb(j));
                end
                
                r2 = rand();
                if r2 < normalized_rime_rates(i-hpop_size)
                    x_new(i,j) = best3(j);
                end
            end

        else
             x_new(i,:) = x(i,:) + rand(1,dim).*(best1 - x(i,:)) - rand(1,dim).*(worst - x(i,:));
        end
        
        x_new(i,:) = max(min(x_new(i,:),ub),lb);
        f_new(i)= fEval(x_new(i,:), fun, opt_f, FuncAddr);
        
        if f_new(i) < f(i)
            x(i,:) = x_new(i,:);
            f(i) = f_new(i);
        end
    end
  
    nfe = nfe + pop_size;
    
       
    ccurve(old_nfe+1:nfe) = min(f);
    % Determine the diversity of the archive    
    prum = mean(x);
    prum_mat=repmat(prum,pop_size,1);
    diam2 = (sum( sum( (x-prum_mat).*(x-prum_mat) ) ) / pop_size);
    diversity = sqrt(diam2);
    dcurve(old_nfe+1:nfe) = diversity;

end % while

best_f = min(f);

ccurve(nfe:maxFEs) = best_f;
dcurve(nfe:maxFEs) = diversity;

end % Jaya

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
  
 function [best] = choose(i, pop_size, x, f)
    
    left = i - 1;
    if (left<1) 
        left = pop_size; 
    end
    
    right = i + 1;
    if (right>pop_size) 
        right = 1; 
    end
    
    xx = [x(left,:);x(right,:)];
    ff = [f(left); f(right)];
    
    [~,best_i] = min(ff);
     
    best = xx(best_i,:);
 
end
