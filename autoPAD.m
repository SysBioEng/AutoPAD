function model=autoPAD(model,pHset,pHzero, pKa,modelName,reportFormat,directionBool)
%autoPAD (AUTOmatic pH ADjustment) algorithm adjusts the charge and chemical
%formula of the metabolites present in a genome-scale metabolic model to the
%formula and charge of the protonation state that is more abundant at the pH
%set fixed by the user. Then, the algorithms balances the reactions included
%in the model that are imbalanced due to these mass and charge
%modifications.
%autoPAD avoids: exchange rxns, and reactions that include metabolites that
%don't have its chemical formula in the model
%
%USAGE:
%      model=autoPAD(model,pHset,pHstart, pKa)  
%
%INPUTS:
%model          COBRA model structure. Requires vectors metFormula and
%               metCharges. vector comps is desired
%pHset          list of pH values to set for each compartment [cx1]
%pHzero        list of the current pH values for each compartment [cx1] 
%pKa            numeric array. Contains row-wise all the pKa related to the
%               metabolites in the model. Each row correspond to a metabolite
%
%OPTIONAL INPUT:
%modelName      name of the COBRA model structure set to the pHset values
%               selected by the user (default: adjustedModel)
%reportFormat   determines the format of the printed document. Options: 0
%               for no printed archive, 1 for excel archive; 2 for text file
%               archive. default: 1
%directionBool  boolean array that indicates which transport reactions are
%               directed on the opposite direction as defined in the model,
%               valid for reversible reactions (default: none is modified)
%
%OUTPUTS:
%
%model          A COBRA model structure that has each compartment balanced
%               at the pHset defined by the user
%
%
%Authors:
%- Magdalena Ribbeck 1/18

%%
%INPUT VALIDATIONS
[model,metComps]=assignCompartments(model);

if length(model.comps) ~= length(pHset)
    error('number of compartments is different from pHset values');
end

if length(model.comps) ~= length(pHzero)
    error('number of compartments is different from initial pH values');
end

if not(isnumeric(pHset)) || not(isnumeric(pHset))
    error('pHset and pHzero should be numeric arrays');
end

if not(isnumeric(pKa))
    error('pKa should be a numeric array');
end

if length(pKa(:,1))~=length(model.mets)
    error('pKa table size is different to the quantity of mets in the model');
end

if nargin<5
    modelName='adjustedModel';
end

if not(isstr(modelName))
    error('modelName should be a string');
end

if nargin<6
    reportFormat=1;
end

if not(isnumeric(reportFormat))
    error('reportFormat should be a numeric input');
elseif reportFormat~=0 && reportFormat~=1 && reportFormat~=2
    error('reportFormat should be a numeric input in the range 0, 1 or 2');
end

if nargin<7
    directionBool=zeros(length(model.rxns),1);
end

if not(isnumeric(directionBool))
    error('directionBool should be a numeric array');
end

if length(directionBool)~=length(model.rxns)
    error('length of directionBool should be the same as the reactions present in the model');
end

sum=length(find(directionBool==0))+ length(find(directionBool==1));
if sum~=length(directionBool);
    error('directionBool can only contain 0s and 1s');
end

pKa=sort(pKa,2,'descend');

%%
%PART 1: ADJUSTMENT OF THE CHEMICAL AND CHARGE FORMULAS OF EACH METABOLITE

%DIVISION OF EACH CHEMICAL FORMULAS PER ELEMENT
%Elements = {'H','C', 'O','P','S', 'N','Mg','X','Fe','Zn','Co','R','K','Cl','Cd','Na','Ni','Mn', 'Cu'};
Elements = {'H','C', 'O','P','S', 'N','Mg','X','Fe','Zn','Co','R','K','Cl','Cd','Na','Ni','Mn','Cu','Ca','Y','I','F','Ag','FULLR'};

MetabElements=zeros(length(model.metNames),length(Elements));

for i=1:length(model.metNames)
    for j=1:length(Elements)
        if not(isempty(model.metFormulas{i}))
            MetabElements(i,j)=numAtomsOfElementInFormula(model.metFormulas{i},Elements{j});
        end
    end
end


%%
%PH ADJUSTMENT OF EACH METABOLITE

posbasePKA=zeros(length(model.mets),1);
posmodifiedPKA=zeros(length(model.mets),1);

for i=1:length(model.mets)
    posbasePKA(i)=findpH(pKa(i,:),pHzero(metComps(i)));
    posmodifiedPKA(i)=findpH(pKa(i,:),pHset(metComps(i)));
end

%to each metabolite, protons are added or taken when it corresponds.
protonVariation=posmodifiedPKA-posbasePKA;
MetabElements(:,1)=MetabElements(:,1)+protonVariation;


%%
% REFORMULATION OF CHEMICAL FORMULAS AND CHARGES


model.metCharges=model.metCharges+protonVariation;

metAdjustedFormulas=cell(length(model.mets),1);

for i=1:length(model.mets)
    for j=1:length(Elements)
        if MetabElements(i,j)>0
            metAdjustedFormulas{i}=strcat(metAdjustedFormulas{i},Elements{j},num2str(MetabElements(i,j)));
        end
    end
end

%model.metFormulas=metAdjustedFormulas;

%adjust the new formulas; if there were not elements found, it means
%it is possible that the metabolite has a ionic formula that was not counted.
%in that case, keep the original formula.
for i=1:length(model.mets)
    if not(isempty(metAdjustedFormulas{i})) 
        model.metFormulas{i}=metAdjustedFormulas{i};
    end
end

%%
%PART 2: ADJUSTMENT OF CHEMICAL REACTIONS, BALANCE OF REACTIONS WHICH ARE
%IMBALANCES DUE TO PROTONS


[massImbalance,imBalancedMass,imBalancedCharge,imBalancedBool,Elements]=checkMassChargeBalance(model,0);
protonImbalance=massImbalance(:,1);

imbalancedBool=findComplexImbalancedRxns(massImbalance);

[model,protonPosition]=identifyProtons(model,metComps);
nonAnalyzedRxns=[];

%Proton modifications according to protonImbalance
for i=1:length(model.rxns)
    compartment=assignRxnToCompartment(model,metComps,imbalancedBool,directionBool,pHset,i);
    if compartment>0
        model.S(protonPosition(compartment),i)=model.S(protonPosition(compartment),i)-protonImbalance(i);
    elseif compartment==0
        a=strcat('Reaction',{' '}, model.rxns(i), ' was not analyzed as not all of the metabolites included in the reaction possess an assigned chemical formula');
        warning(a{1}); 
        nonAnalyzedRxns=[nonAnalyzedRxns;model.rxns(i)];
    end
end


save(modelName,'model');

%%
%LOG GENERATION

if reportFormat==0
    
elseif reportFormat==1
    c=date;
    archiveName=strcat('AUTOPAD log for',{' '},modelName, {' '}, datestr(c));
    archiveName=archiveName{1};
    t=strcat('adjusted model =',{' '},modelName);
    beginning={datestr(c);t{1};'';'PH VALUES'};
    xlswrite(archiveName,beginning,'sheet1');
    xlswrite(archiveName,{'compartments';'pH values'},'sheet1', 'A5');
    xlswrite(archiveName,model.comps','sheet1', 'B5');
    xlswrite(archiveName,pHset,'sheet1', 'B6');
    xlswrite(archiveName,{'REACTIONS'},'sheet1', 'A8');
    xlswrite(archiveName,strcat(num2str(length(model.rxns)-length(nonAnalyzedRxns)), {' '}, 'out of',{' '}, num2str(length(model.rxns)),{' '},'reactions were analyzed and balanced if needed'),'sheet1', 'A9');

    xlswrite(archiveName, strcat('the following', {' '}, num2str(length(nonAnalyzedRxns)),{' '},'reactions were not, due to lack of chemical formula of one of the metabolites that participate in the reaction,'),'sheet1', 'A10');
    xlswrite(archiveName,{'or due to complex mass imbalances:'},'sheet1', 'A11');
    xlswrite(archiveName,nonAnalyzedRxns,'sheet1', 'A12');

    disp('AUTOPAD log was saved in an excel file');

else
    c=date;
    archiveName=strcat('AUTOPAD log for',{' '},modelName, {' '}, datestr(c),'.txt');
    archiveName=archiveName{1};
    t=strcat('adjusted model =',{' '},modelName);
    beginning={datestr(c);t{1};'';'PH VALUES'};
    fileID=fopen(archiveName,'w');
    fprintf(fileID,'%4s\r\n',beginning{1});
    fprintf(fileID,'%4s\r\n',beginning{2});
    fprintf(fileID,'%4s\r\n',beginning{3});
    fprintf(fileID,'%4s\r\n',beginning{4});
    fprintf(fileID,'%4s\t','compartments');
    for i=1:length(model.comps)-1
        fprintf(fileID, '%4s\t',model.comps{i});
    end
    fprintf(fileID, '%4s\n',model.comps{length(model.comps)});
    
    fprintf(fileID,'%4s\t','pH values');
    for i=1:length(pHset)-1
        fprintf(fileID, '%4d\t',pHset(i));
    end
    fprintf(fileID, '%4d\n',pHset(length(pHset)));
    
    fprintf(fileID,'%4s\r\n','REACTIONS');
    temp=strcat(num2str(length(model.rxns)-length(nonAnalyzedRxns)), {' '}, 'out of',{' '}, num2str(length(model.rxns)),{' '},'reactions were analyzed and balanced if needed');
    fprintf(fileID,'%4s\r\n',temp{1});
    temp=strcat('the following', {' '}, num2str(length(nonAnalyzedRxns)),{' '},'reactions were not, due to lack of chemical formula of one of the metabolites that participate in the reaction, or due to complex mass imbalances:');
    fprintf(fileID,'%4s\r\n',temp{1});
    
    for i=1:length(nonAnalyzedRxns)
        fprintf(fileID, '%4s\n',nonAnalyzedRxns{i});
    end
    fclose(fileID);

    disp('AUTOPAD log was saved in a text file');
end


end

function pos=findpH(pKa,pH)
%for one metabolite

    if isnan(pKa(length(pKa))) %if there are not pKa for the metabolite
        pos=0;
    else
        for i=1:length(pKa)
            if not(isnan(pKa(i)))
                if pH>pKa(i)
                    pos=i-1;
                    break;   
                elseif i==length(pKa)
                    pos=i;
                    break;
                end
            end
        end
    end

end

function [model,metCompsNumeric]=assignCompartments(model)
%assigns a compartment for each metabolite. If comps does not exist,
%creates comps too

%
comps={};

metComps={};

%check whether metabolite namimg format is [c] or _c
parenthesisFormat=checkMetaboliteFormat(model.mets);

for i=1:length(model.mets)
    
    if parenthesisFormat
        %find position of last '[' and ']'
        pos1=strfind(model.mets(i),'[');
        pos1=pos1{end}(end);
        pos2=strfind(model.mets(i),']');
        pos2=pos2{end}(end);
        metComps{i}=substrings(model.mets{i},pos1+1,pos2-1);
    else
        %find position of last "_"
        pos=strfind(model.mets(i),'_');
        pos=pos{end}(end);
        metComps{i}=substrings(model.mets{i},pos+1,length(model.mets{i}));
    end
    
    
    isThere=strmatch(metComps{i},comps);
    if isempty(isThere)
        comps=[comps; metComps{i}];
    end
    
end

%check whether compartments identified in this step are the same compartments
%that those indicated by the user in the model 

if isfield(model,'comps')
    if length(comps)>length(model.comps)
        error('AutoPAD identified more compartments than the compartments indicated in model.comps');
    elseif length(comps)<length(model.comps)
        error('AutoPAD identified less compartments than the compartments indicated in model.comps');
    else
        for i=1:length(comps)
            if not(ismember(comps{i},model.comps))
                error('AutoPAD identified different compartments than the compartments indicated in model.comps');
            end
        end
    end 
else
    warning('model did not include model.comps. Vector comps was generated automatically');
    comps
    model.comps=comps;
end

metComps=metComps';
metCompsNumeric=zeros(length(metComps),1);
for i=1:length(metComps)
    for j=1:length(model.comps)
        if strcmp(metComps{i},model.comps{j})
            metCompsNumeric(i)=j;
            break
        end  
    end
    
end



end

function [model,protonPosition]=identifyProtons(model,metComps)
%generates compartment-proton associations, and checks whether each
%compartment possess a proton. If not, the metabolite is generated.
%Returns the metabolite number in the model that each proton ocupies in the
%same order as the compartments.

%check whether metabolite namimg format is [c] or _c
parenthesisFormat=checkMetaboliteFormat(model.mets);

protons=[];

%select proton candidates
if parenthesisFormat
    protons=startWith(model.mets,'h[');
    protons=[protons; startWith(model.mets,'H[')];
    protons=[protons; startWith(model.mets,'h [')];
    protons=[protons; startWith(model.mets,'H [')];
    protons=[protons; startWith(model.mets,'h+[')];
    protons=[protons; startWith(model.mets,'H+[')];
    protons=[protons; startWith(model.mets,'h(+)[')];
    protons=[protons; startWith(model.mets,'H(+)[')];
    protons=[protons; startWith(model.mets,'h+ [')];
    protons=[protons; startWith(model.mets,'H+ [')];
    protons=[protons; startWith(model.mets,'h(+) [')];
    protons=[protons; startWith(model.mets,'H(+) [')];
else
    protons=startWith(model.mets,'h_');
    protons=[protons; startWith(model.mets,'H_')];
    protons=[protons; startWith(model.mets,'h+_')];
    protons=[protons; startWith(model.mets,'H+_')];
    protons=[protons; startWith(model.mets,'h(+)_')];
    protons=[protons; startWith(model.mets,'H(+)_')];
end


%check whether all compartments have its own proton for balancing. if not,
%create new met.

protonPosition=[];
if length(protons)>length(model.comps)
    error('Failure at proton identification, too many protons were found. Check the correct name assignment in mets');
    
elseif length(protons)<length(model.comps)
    for i=1:length(model.comps)
        isThere=0;
        for j=1:length(protons)
            pos=strmatch(protons{j},model.mets);
            compartment=metComps(pos);
            
            if compartment==i
                isThere=1;
                protonPosition=[protonPosition; pos];
                break
            end
            
        end
        
        if not(isThere)
            if parenthesisFormat
                metName=strcat('h[',model.comps{i},']');
            else
                metName=strcat('h_',model.comps{i});
            end
            model=addMetabolite(model,metName,metName,'H','','','','',1,0);
            protons=[protons; metName];
            protonPosition=[protonPosition; length(model.mets)];
            metComps=[metComps;i];
            disp(strcat('Due to absence in the compartment, a proton was added to compartment',{' '},model.comps{i},' named as',{' '},metName));
        end        
    end
else %equal amount of compartments and protons
    for k=1:length(model.comps)
        isThere=0;
        for l=1:length(protons)
            pos=strmatch(protons{l},model.mets);
            compartment=metComps(pos);
            if compartment==k
                protonPosition=[protonPosition;pos];
                isThere=1;
                break;
            end
        end
        if not(isThere)
            error('Failure at proton identification, too many protons were found. Check the correct name assignment in mets');
        end
    end
end


end

function imbalancedBool=findComplexImbalancedRxns(massImbalance)
%determines which reactions present complex imbalances, i.e, imbalances
%that are not due to protons but caused by other elements. 

massImbalance=massImbalance(:,[2:end]);
imbalancedBool=zeros(length(massImbalance(:,1)),1);

for i=1:length(massImbalance(:,1))
   if not(isempty(find(massImbalance(i,:))))
      imbalancedBool(i)=1;
   end 
end

end

function compartment=assignRxnToCompartment(model,metComps,imbalancedBool,directionBool,pHset,rxn)
%assigns each reaction to a compartment, thus any imbalanced reactions will
%be rebalanced ussing the proton related to that specific compartment.
%returns the number of the compartment in which the reaction is assigned.
%if reaction is exchange, returns -1. if reaction is mass imbalanced in other
%elements and thus was not analyzed, returns 0.

subPositions=find(model.S(:,rxn)<0);
prodPositions=find(model.S(:,rxn)>0);
subPlacedIn=metComps(subPositions);
prodPlacedIn=metComps(prodPositions);

%if rxn is exchange, it is already balanced
if isempty(subPositions) || isempty(prodPositions)
    compartment=-1;
    return
end

%if reaction presents imbalances in other elements, it is not analyzed
if imbalancedBool(rxn)
    compartment=0;
    return
end

%the same if any metFormula is missing
for i=1:length(subPositions)
    if isempty(model.metFormulas{subPositions(i)})
         compartment=0;
         return
    end
end
for i=1:length(prodPositions)
    if isempty(model.metFormulas{prodPositions(i)})
         compartment=0;
         return
    end
end


%classification of the reactions:
%this classification change to the opposite if reversibility is indicated
%by the user

%A(c)+B(c)->C(c)+D(c) same compartment
if all(subPlacedIn == subPlacedIn(1)) && all(prodPlacedIn == prodPlacedIn(1)) && subPlacedIn(1) == prodPlacedIn(1)
    if not(directionBool)
        compartment=subPlacedIn(1);
    else
        compartment=prodPlacedIn(1);
    end
%A(c)+B(c)->A(e)+B(e) simport/diffusion
elseif all(subPlacedIn == subPlacedIn(1)) && all(prodPlacedIn == prodPlacedIn(1))
    if not(directionBool)
        compartment=prodPlacedIn(1);
    else
        compartment=subPlacedIn(1);
    end
    
%A(c)+B(e)->C(c)+B(c)+D(c) ABC
elseif all(prodPlacedIn == prodPlacedIn(1))
    if not(directionBool)
        compartment=prodPlacedIn(1);
    else
        compartment=subPlacedIn(1);
    end
%A(c)+B(c)->A(e)+B(c) membrane bound rxns    
elseif all(subPlacedIn == subPlacedIn(1))
    if not(directionBool)
        compartment=subPlacedIn(1);
    else
        compartment=prodPlacedIn(1);
    end
    
%A(c)+B(e)->A(e)+B(c) OR A(c)+H(e)->H(c)+A(e)  antiport
else
    %select compartment with the lower pH of the possible compartments, as
    %in this one elements associate.
    lower=min(pHset(subPlacedIn));
    
    if not(directionBool)
        compartment=find(pHset==lower);
    else
        compartment=find(pHset==max(pHset(subPlacedIn)));
    end
end

end

function parenthesisFormat=checkMetaboliteFormat(mets)
%checks whether metabolite naming format is [c] or _c

parenthesisFormat=0;
checkFormat=strtrim(mets{1});
if strcmp(checkFormat(end),']')
    parenthesisFormat=1;
end

end

function sub=substrings(string,start,ending)
length=ending-start+1;

sub='';
for i=start:ending
    sub=[sub, string(i)];  
end

end

function candidates=startWith(array,element)

candidates=[];
for i=1:length(array)
    if length(array{i})>=length(element)
        if strcmp(substrings(array{i},1,length(element)),element)
            candidates=[candidates;array(i)];
        end
    end
end

end

