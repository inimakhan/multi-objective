

%__________________________________________________________________ %
%                          Multi-Objective                          %
%        Crystal Structure Algorithm (CryStAl) (MOCryStAl)          %
%                                                                   %
%                                                                   %
%                  Developed in MATLAB R2021a (MacOs)               %
%                                                                   %
%                      Author and programmer                        %
%                ---------------------------------                  %
%                      Nima Khodadadi (ʘ‿ʘ)                         %
%                       Siamak Talatahari                           %
%                         Mahdi Azizi                               %
%                         Pooya Sareh                               %
%                                                                   %
%                             e-Mail                                %
%                ---------------------------------                  %
%                         inimakhan@me.com                          % 
%                                                                   %
%                            Homepage                               %
%                ---------------------------------                  %
%                    https://nimakhodadadi.com                      %
%                                                                   %
%                                                                   %
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------------------------------------------------------- %

function [Archive_costs]=MOCRY(MaxIteation,Archive_size,Cr_Number,Var_Number,method,m)

disp('MOCRY is running')

%% Simple Problem Definition


if method==3
    
    TestProblem=sprintf('P%d',m);
    
    ObjFuncName= Ptest(TestProblem);
    
    xrange  = xboundaryP(TestProblem);
    Var_Number=max(size(xrange));
    % Lower bound and upper bound
    LB=xrange(:,1)';
    UB=xrange(:,2)';
    
end


%% Updating the Size of ProblemParameters

alpha=0.1;  % Grid Inflation Parameter
nGrid=30;   % Number of Grids per each Dimension
beta=4; %=4;    % Leader Selection Pressure Parameter
gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure


if length(LB)==1
    LB=repmat(LB,1,Var_Number);
end
if length(UB)==1
    UB=repmat(UB,1,Var_Number);
end
%% Initialization


Crystal=CreateEmptyParticle(Cr_Number);

% Initializing the Position of first probs
for i=1:Cr_Number
        
    Crystal(i).Velocity=0;
    Crystal(i).Position=zeros(1,Var_Number);
    for j=1:Var_Number
        Crystal(i).Position(1,j)=unifrnd(LB(j),UB(j),1);
    end
    Crystal(i).Cost=ObjFuncName(Crystal(i).Position')';
    
    Fun_eval(i,:)=norm(Crystal(i).Cost);
    Crystal(i).Best.Position=Crystal(i).Position;
    Crystal(i).Best.Cost=Crystal(i).Cost;
end

[~,idbest]=min(Fun_eval);
Crb=Crystal(idbest,:);
Crystal=DetermineDominations(Crystal);
Archive=GetNonDominatedParticles(Crystal);

Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end
% The best Crystal

%% Search Process
Iter=1;
while Iter<MaxIteation
    for i=1:Cr_Number
        Leader=SelectLeader(Archive,beta);
        %% Generate New Crystals
        % Main Crystal
        Crmain=Crystal(randperm(Cr_Number,1),:);
        % Random-selected Crystals
        RandNumber=randperm(Cr_Number,1);
        RandSelectCrystal=randperm(Cr_Number,RandNumber);
        CrystalP=vertcat( Crystal.Position );
        % Mean of randomly-selected Crystals
        Fc=mean(CrystalP(RandSelectCrystal,:)).*(length(RandSelectCrystal)~=1)...
            +CrystalP(RandSelectCrystal(1,1),:).*(length(RandSelectCrystal)==1);
        % Random numbers (-1,1)
        r=2*rand-1;       r1=2*rand-1;
        r2=2*rand-1;     r3=2*rand-1;
        % New Crystals
        Crystal(1+Cr_Number,:).Position=Leader.Position+r*Crmain.Position;
        Crystal(2+Cr_Number,:).Position=Leader.Position+r1*Crmain.Position+r2*Crb.Position;
        Crystal(3+Cr_Number,:).Position=Leader.Position+r1*Crmain.Position+r2*Fc;
        Crystal(4+Cr_Number,:).Position=Leader.Position+r1*Crmain.Position+r2*Crb.Position+r3*Fc;
        
        for i2=1:4
            % Checking/Updating the boundary limits for Crystals
            Crystal(i2+Cr_Number,:).Position=min(max(Crystal(i2+Cr_Number).Position,LB),UB);
            
            % Evaluating New Crystals
            Crystal(i2+Cr_Number,:).Cost=ObjFuncName(Crystal(i2+Cr_Number).Position')';
            Crystal(i2+Cr_Number).Best.Position=Crystal(i2+Cr_Number).Position;
            Crystal(i2+Cr_Number).Best.Cost=Crystal(i2+Cr_Number).Cost;
            % Updating the Crystals
            
        end
    end % End of One Iteration
    
    Crystal=DetermineDominations(Crystal);
    non_dominated_Crystal=GetNonDominatedParticles(Crystal);
    
    Archive=[Archive
        non_dominated_Crystal];
    
    Archive=DetermineDominations(Archive);
    Archive=GetNonDominatedParticles(Archive);
    
    for i=1:numel(Archive)
        [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end
    
    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);
        
        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,alpha);
        
    end
    disp(['In iteration ' num2str(Iter) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
    save results
    
    
    Archive_costs=GetCosts(Archive);
    
    
    Iter=Iter+1;
    % The best Crystal
    
end % End of Main Looping
end

%% Boundary Handling
function x=bound(x,UB,LB)
x(x>UB)=UB(x>UB); x(x<LB)=LB(x<LB);
end
