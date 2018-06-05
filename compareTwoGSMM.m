function compareTwoGSMM(model1,model2,reportFormat,archiveName,slideName)
%Compares two models that contain the same reactions and checks whether there
%is a stoichiometric difference in these reactions, and generates an excel
%archive with the reactions that present a difference
%
%USAGE:
%               compareTwoGSMM(model1,model2,archiveName)
%
%INPUTS:
%model1         COBRA model structure that possess model.metKEGGID
%model2         COBRA model structure that possess model.metKEGGID
%
%OPTIONAL INPUTS:
%reportFormat   determines the format of the printed document. Options: 
%               1 for excel archive; 2 for text file archive. default: 1
%archiveName    Name of the excel archive (.xls) in which the archive is
%               printed, default: 'modified reactions [current date]'
%slideName      Name of the excel slide in which the archive is printed,
%               default: 'sheet1'
%
%
%Authors:
%- Magdalena Ribbeck 1/18

%%
%INPUT VALIDATION

if nargin<1
    error('two models are required as input');
end

if nargin<3
    reportFormat=1;
end
if not(isnumeric(reportFormat))
    error('reportFormat should be a numeric input');
elseif reportFormat~=1 && reportFormat~=2
    error('reportFormat should be a numeric input in the range 1 or 2');
end


if nargin<4
    c=date;
    archiveName=strcat('modified reactions',{' '},datestr(c));
    archiveName=archiveName{1};
end
if not(isstr(archiveName))
    error('archiveName should be a string');
end

if nargin<5
   slideName='sheet1';
end

if not(isstr(slideName))
    error('slideName should be a string');
end

%%
%ALGORITHM
pos=0;
name=cell(1,1);
rxn1=cell(1,1);
rxn2=cell(1,1);

for i=1:length(model1.rxns)
    
    posModel2=find(strcmp(model1.rxns(i),model2.rxns));
    
    %obtain the metabolites that participate in the rxn
    [mets1 t coefs1]=find(model1.S(:,i));
    [mets2 t coefs2]=find(model2.S(:,posModel2));

    %compare the metabolites that participate in the rxn
    notEqual=0;
    if length(mets1) == length(mets2)
        for k=1:length(mets1)
            if isempty(find(strcmp(model1.metNames(mets1(k)),model2.metNames(mets2))))
                notEqual=1;
                %if one of the mets is different, rxns are too
            end
        end
    else
        notEqual=1; %or if different quantity of mets
    end

    if notEqual==0 && isequal(coefs1,coefs2)
    else
        pos=[pos; i];
        name=[name; model1.rxns{i}];
        rxn1=[rxn1; printRxnFormula(model1,model1.rxns(i),0,1,1)];
        rxn2=[rxn2; printRxnFormula(model2,model2.rxns(posModel2),0,1,1)];
         
    end
end

%%
%LOG GENERATION

if reportFormat==1
    xlswrite(archiveName,[name rxn1 rxn2],slideName);
    xlswrite(archiveName,{'reactions','model1','model2'},slideName,'A1');
    disp('comparison ready and available in excel file');
else
    fileID=fopen(archiveName,'w');
    fprintf(fileID,'%4s\t','reactions');
    fprintf(fileID,'%4s\t','model1');
    fprintf(fileID,'%4s\r\n','model2');
    for i=2:length(name)
        fprintf(fileID,'%4s\t',name{i});
        fprintf(fileID,'%4s\t',rxn1{i});
        fprintf(fileID,'%4s\n',rxn2{i});
    end
    fclose(fileID);
    disp('comparison ready and available in text file');
end


end