classdef SSCMOSWO < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Multi-objective spider wasp optimization algorithm based on sparse species conservation

%------------------------------- Reference --------------------------------
% A multi-objective spider wasp algorithm based on sparse species conservation for the 
%application in flexible charging of household electric vehicles
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);
Individual.Position=[];
Individual.Cost=[];
Individual.Dominated=false;
Individual.CrowdingDistance=0;
Individual.Rank=[];
Individual.DominationSet=[];
Individual.DominatedCount=[];
Pop=repmat(Individual,Problem.N,1);
for i=1:Problem.N 
Pop(i).Position=Population(i).dec;
Pop(i).Cost=Population(i).obj;
end 
popm=Pop;
            %Optimization
            while Algorithm.NotTerminated(Population)
                Population1 = SWO(Problem,Population,Pop,popm);
                Zmin       = min([Zmin;Population1(all(Population1.cons<=0,2)).objs],[],1);
                Population = NDSSSC([Population,Population1],Problem.N,Z,Zmin);
            end
        end
    end
end