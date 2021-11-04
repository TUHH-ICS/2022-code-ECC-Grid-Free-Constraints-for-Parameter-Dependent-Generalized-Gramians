
function gramian = calcGramianBoundFBSP(sysUbLft,blkStructPlant,uncertainty,gramianFlag,gramWeight,maxRate,onlyFeasibility)
%CALCGRAMIANBOUNDFBSP Computes affine generalized Gramian in gridless manner
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    % System specs
    n = order(sysUbLft);
    numParams = numel(blkStructPlant);

    numBases = size(gramWeight,1) / n;
    gramOrder = (numBases - 1) / numParams;


    % Kind of Gramian
    if(strcmp(gramianFlag,'c'))

        calcCtrb = true;
    elseif(strcmp(gramianFlag,'o'))

        calcCtrb = false;

    else
        error('Flag unspecified')
    end


    % Matrix dimensions
    nDPlant = 0;
    for i = 1:numParams

        nDPlant = sum([nDPlant, blkStructPlant(i).Occurrences]);
    end

    nU = size(sysUbLft.B,2) - nDPlant;
    nY = size(sysUbLft.C,1) - nDPlant;


    % Set components of plant LFT 
    A0 = sysUbLft.A;
    B1 = sysUbLft.B(:,1:nDPlant);
    C1 = sysUbLft.C(1:nDPlant,:);
    D11 = sysUbLft.D(1:nDPlant,1:nDPlant);

    if(calcCtrb)

        B0 = sysUbLft.B(:,nDPlant+1:end);
        D10 = sysUbLft.D(1:nDPlant,nDPlant+1:end);

        G{1,1} = D11';
        G{1,2} = B1';
        G{2,1} = [zeros(n,nDPlant);...
                  C1';...
                  D10'];
        G{2,2} = [eye(n);...
                  A0';...
                  B0'];

        identity = eye(nU);
    else

        C0 = sysUbLft.C(nDPlant+1:end,:);
        D01 = sysUbLft.D(nDPlant+1:end,1:nDPlant);

        G{1,1} = D11;
        G{1,2} = C1;
        G{2,1} = [zeros(n,nDPlant);...
                  B1;...
                  D01];
        G{2,2} = [eye(n);...
                  A0;...
                  C0];

        identity = eye(nY);
    end


    % Determin the grammian
    % Init yalmip
    yalmip('clear');
    sdpOptions = sdpsettings('solver', 'mosek','verbose',0);


    % Define decision variables

    if(gramOrder == 0)    

        % Gramian
        X = sdpvar(n,n,'symmetric');        % >0

        QRepeated = blkdiag(sdpvar(nDPlant,nDPlant,'symmetric'));       % < 0
        SRepeated = blkdiag(sdpvar(nDPlant,nDPlant,'skew'));
        UpsilonTilde = ones(nDPlant,1);

        listParams = createParameterList(blkStructPlant,nDPlant);

        
        for i = 1:size(QRepeated,1)
            
            % Eliminate couplings between different parameters
            for j = 1:size(QRepeated,2)

                if(~strcmp(listParams(i),listParams(j)))

                    QRepeated(i,j) = 0;
                    SRepeated(i,j) = 0;
                end
            end
            
            UpsilonTilde(i) = 1; % Assume a scaled LFT model
        end
        
        UpsilonTilde = diag(UpsilonTilde);

    else
        % Extend list of uncertainties by change rates
        uncertaintyExt = struct2cell(uncertainty);
        fieldNamesExt = fieldnames(uncertainty);
        
        for i = 1:numel(uncertaintyExt)
            
            fieldNamesExt{end+1,1} = [uncertaintyExt{i}.Name 'd'];
            uncertaintyExt{end+1,1} = ureal(fieldNamesExt{end},0,'Range',[-maxRate(i), maxRate(i)]);
        end
        uncertaintyExt = cell2struct(uncertaintyExt,fieldNamesExt);
        

        %Gramian Bases
        X = sdpvar(n*numBases,n*numBases,'symmetric');
        X(1:end-n,1:end-n) = zeros(size(X,1)-n);

        %Gramian LFT
        fieldNames = fieldnames(uncertainty);

        Trho = eye(n);

        for i = 1:numel(fieldNames)        

            Trho = [uncertainty.(fieldNames{i}) * eye(n); Trho];
        end

        [M,~,blkStructGramian] = lftdata(Trho);     %ssbal()?

        T{1,1} = M(1:end-numBases*n,1:end-n);
        T{1,2} = M(1:end-numBases*n,end-n+1:end);
        T{2,1} = M(end-numBases*n+1:end,1:end-n);
        T{2,2} = M(end-numBases*n+1:end,end-n+1:end);


        nDGramian = size(T{1,1},1);
        
        THat{1,1} = kron([0 0 1;1 1 0;0 0 1],T{1,1});

        THat{1,2} = kron([1 0;0 1;1 0],T{1,2});
        THat{1,2} = [THat{1,2}, zeros(size(THat{1,2},1),size(identity,1))];

        if(calcCtrb)

            THat{2,1} = kron([0 0 1;-1 1 0],T{2,1});
        else

            THat{2,1} = kron([0 0 1;1 1 0],T{2,1});
        end
        THat{2,1} = [THat{2,1}; zeros(size(identity,1),size(THat{2,1},2))];

        THat{2,2} = blkdiag(T{2,2},T{2,2},identity);


        % Combine plant and gramian LFT
        GHat{1,1} = [THat{1,1},                     THat{1,2}*G{2,1};...
                     zeros(nDPlant,size(THat{1,1},2)),   G{1,1}];
        GHat{1,2} = [THat{1,2}*G{2,2};...
                     G{1,2}];
        GHat{2,1} = [THat{2,1}, THat{2,2}*G{2,1}];
        GHat{2,2} = THat{2,2}*G{2,2};

        %Overwrite
        G = GHat;


        %Multiplier1
        QRepeated = blkdiag(sdpvar(3*nDGramian+nDPlant,3*nDGramian+nDPlant,'symmetric'));
        SRepeated = blkdiag(sdpvar(3*nDGramian+nDPlant,3*nDGramian+nDPlant,'skew'));

        listParams = [createParameterList(blkStructGramian,nDGramian,'d');
                      createParameterList(blkStructGramian,nDGramian);
                      createParameterList(blkStructGramian,nDGramian);
                      createParameterList(blkStructPlant,nDPlant)];

        for i = 1:size(QRepeated,1)
            
            % Eliminate couplings between different parameters
            for j = 1:size(QRepeated,2)

                if(~strcmp(listParams(i),listParams(j)))

                    QRepeated(i,j) = 0;
                    SRepeated(i,j) = 0;
                end
            end
            
            if i<nDGramian+1
                
                UpsilonTilde(i) = abs(uncertaintyExt.(listParams{i}).Range(1));
            else
                
                UpsilonTilde(i) = 1; % Assume a scaled LFT model
            end
        end
        
        UpsilonTilde = diag(UpsilonTilde);


        QBlock = [];
        SBlock = [];
        UpsilonBlock = [];

        for i = 1:numParams

            QBlock = blkdiag(QBlock,sdpvar(blkStructGramian(i).Occurrences,blkStructGramian(i).Occurrences,'symmetric'));
            SBlock = blkdiag(SBlock,sdpvar(blkStructGramian(i).Occurrences,blkStructGramian(i).Occurrences,'skew'));
            UpsilonBlock = [UpsilonBlock; ones(blkStructGramian(i).Occurrences,1)];
        end
        
        UpsilonBlock = diag(UpsilonBlock);
    end


    % Define LMIs

    % Gramian must solve LMI
    PiRepeated = [UpsilonTilde*QRepeated*UpsilonTilde,     SRepeated;
                  SRepeated',                             -QRepeated];

    M = blkdiag([zeros(size(X,1)), X; X, zeros(size(X,1))],identity);

    calSCon = [G{1,1},          G{1,2};
               eye(size(PiRepeated,1)/2), zeros(size(PiRepeated,1)/2,n);
               G{2,1},          G{2,2}];

    nominalCon = calSCon' * blkdiag(PiRepeated, M) * calSCon;

    LMIConstraints = [nominalCon <= 0, QRepeated >= 0];

    % Gramian must be positive definite

    if(gramOrder == 0)

        LMIConstraints = [LMIConstraints, X >= 0];
    else
        Pi = [-UpsilonBlock*QBlock*UpsilonBlock,   -SBlock;
              -SBlock',                             QBlock];

        calS = [T{1,1},  T{1,2};
                eye(nDGramian), zeros(nDGramian,n);
                T{2,1},  T{2,2}];

        nominal = calS' * blkdiag(Pi,X) * calS;

        LMIConstraints = [LMIConstraints, nominal >= 0, QBlock >= 0];
    end

    % Solve LMIs

    if(~exist('onlyFeasibility','var'))

        onlyFeasibility = false;
    end

    if(onlyFeasibility)

        diagnostic = optimize(LMIConstraints,[],sdpOptions);    % Check feasibility
    else

        diagnostic = optimize(LMIConstraints,costFcn(X,gramWeight,uncertainty),sdpOptions);     % Find the tightest bound
    end


    if(isempty(strfind(diagnostic.info,'Successfully solved (')))

        gramian = NaN;
    else

        gramian = value(X);
    end
end


%%

function list = createParameterList(blkStruct,nD,appendix)

    if(~exist('appendix','var'))
        
        appendix = [];
    end

    numParams = numel(blkStruct);

    list = cell(nD,1);
    
    index = 1;
    
    for i=1:numParams
       
        for j=1:blkStruct(i).Occurrences
        
            list{index,1} = [blkStruct(i).Name appendix];
            index = index+1;
        end
    end
end