function [Population1] = SWO(Problem,Population1,Pop,popm)
%-----------------------------------------------------------------------------------------------------------
%This code is based on the MATLAB implementation of the single-objective spider wasp optimization algorithm. 
%The corresponding reference for the single-objective version is as follows:
%Abdel-Basset M, Mohamed R, Jameel M, et al. Spider wasp optimizer: 
%A novel meta-heuristic optimization algorithm[J]. Artificial Intelligence Review, 2023, 56(10): 11675-11738.
%------------------------------------------------------------------------------------------------------------
TR=0.4;
Cr=0.6;
[Pop, F]=NonDominatedSorting(Pop);
F1 = Pop(F{1});
t=Problem.FE/Problem.N;
MaxIt=Problem.maxFE/Problem.N;

 for i=1:Problem.N
    a=2-2*(t/MaxIt); 
    a2=-1+-1*(t/MaxIt);
     k=(1-t/MaxIt); 
    JK=randperm(Problem.N); 
    Dim=Problem.D;
     DirectVector=zeros(Problem.N,Dim);
  
     for j = 1:Dim
            DirectVector(i, j) = randi([0 1]);
     end
     
      if rand<TR
            r1=rand(); 
            r2=rand(); 
            r3=rand();
            p = rand();  
            C=a*(2*r1-1); 
            l=(a2-1)*rand+1;  
            L=Levy(1);
            vc = unifrnd(-k,k,1,Dim);
           
            if i<k*Problem.N
               if p<(1-t/MaxIt) 
                    if r1<r2
                       popm(i).Position=F1(randi([1 numel(F1)])).Position+randn*DirectVector(i,:).*(F1(randi([1 numel(F1)])).Position-Pop(i).Position);
                    else          
                           popm(i).Position=(Problem.lower+rand*(Problem.upper-Problem.lower));  % Eq. (6)                
                    end   
                else 
                    if r1<r2
                        popm(i).Position=Pop(i).Position+C*abs(2*rand*Pop(JK(3)).Position-Pop(i).Position); %% Eq. (10)
                    else
                        popm(i).Position=Pop(i).Position.*vc;
                    end 
                end
            else  
                if r1<r2
                    popm(i).Position=F1(randi([1 numel(F1)])).Position+cos(2*l*pi)*(F1(randi([1 numel(F1)])).Position - Pop(i).Position);
                else
                    popm(i).Position=Pop(JK(1)).Position+r3*abs(L)*(Pop(JK(1)).Position-Pop(i).Position)+(1-r3)*(rand>rand)*(Pop(JK(3)).Position-Pop(JK(2)).Position);
                end 
            end  
            for j=1:Dim
                    popm(i).Position(j)= max(popm(i).Position(j), Problem.lower(j));
                    popm(i).Position(j)= min(popm(i).Position(j), Problem.upper(j));
            end
             popm(i).Cost=Problem.Evaluation( popm(i).Position);
           
            
        else  
            l=(a2-1)*rand+1;   
            if Dominates(Pop(JK(1)).Cost,Pop(JK(i)).Cost)
                v1=Pop(JK(1)).Position-Pop(JK(i)).Position;
            else
                v1=Pop(JK(i)).Position-Pop(JK(1)).Position;
            end
            if Dominates(Pop(JK(2)).Cost,Pop(JK(3)).Cost) 
                v2=Pop(JK(2)).Position-Pop(JK(3)).Position;
            else
                v2=Pop(JK(3)).Position-Pop(JK(2)).Position;
            end
            rn1=randn; 
            rn2=randn; 
            SW_m= Pop(i).Position+(exp(l))*abs(rn1)*v1+(1-exp(l))*abs(rn2)*v2;      
            if(rand<Cr) 
                popm(i).Position=SW_m;             
            else
                popm(i).Position=F1(randi([1 numel(F1)])).Position+cos(2*l*pi)*(F1(randi([1 numel(F1)])).Position - Pop(i).Position);
            end  
            for j=1:Dim
                    popm(i).Position(j)= max(popm(i).Position(j), Problem.lower(j));
                    popm(i).Position(j)= min(popm(i).Position(j), Problem.upper(j));
            end  
             popm(i).Cost=Problem.Evaluation( popm(i).Position);
      end  
 end

 
 for i=1:Problem.N
 Population1(i)=Problem.Evaluation( popm(i).Position);
 end
 
end
function [pop, F]=NonDominatedSorting(pop)

    nPop=numel(pop);

    for i=1:nPop
        pop(i).DominationSet=[];
        pop(i).DominatedCount=0;
    end
    
    F{1}=[];
    
    for i=1:nPop
        for j=i+1:nPop
            p=pop(i);
            q=pop(j);
            
            if Dominates(p.Cost,q.Cost)
                p.DominationSet=[p.DominationSet j];
                q.DominatedCount=q.DominatedCount+1;
            end
            
            if Dominates(q.Cost,p.Cost)
                q.DominationSet=[q.DominationSet i];
                p.DominatedCount=p.DominatedCount+1;
            end
            
            pop(i)=p;
            pop(j)=q;
        end
        
        if pop(i).DominatedCount==0  
            F{1}=[F{1} i];
            pop(i).Rank=1;
        end
    end
    
    k=1;
    
    while true
        Q=[];
        for i=F{k}
            p=pop(i);
            for j=p.DominationSet
                q=pop(j);             
                q.DominatedCount=q.DominatedCount-1;
                if q.DominatedCount==0
                    Q=[Q j]; 
                    q.Rank=k+1;
                end
                pop(j)=q;
            end
        end
        if isempty(Q)
            break;
        end
        F{k+1}=Q; 
        k=k+1;    
    end
end
function pop=DetermineDomination(pop)

 nPop=numel(pop);

    for i=1:nPop
        pop(i).Dominated=false;
    end
    
    for i=1:nPop
        for j=i+1:nPop
            if Dominates(pop(i).Cost,pop(j).Cost)
              
                pop(j).Dominated=true;
                
            elseif Dominates(pop(j).Cost,pop(i).Cost)
                pop(i).Dominated=true;
                break;
            else
                if    all (pop(i).Cost==pop(j).Cost) %remove the same individual
                      pop(i).Dominated=true;
                end
                
            end
        end
    end

end
function dom=Dominates(x,y)
    
    dom=all(x<=y) && any(x<y);

end
