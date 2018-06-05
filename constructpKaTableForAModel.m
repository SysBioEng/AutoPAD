function [pKaTable,report]=constructpKaTableForAModel(model,archiveName,archiveFormat)
%generates a model-specific pKa table based on the Kegg Database pKa table.
%
%USAGE:
%           [pKaTable,report]=constructpKaTableForAModel(model,archiveName)
%
%INPUTS:
%model          COBRA model structure that possess model.metKEGGID
%
%OPTIONAL INPUTS:
%archiveName    Name of the excel archive (.xls) in which the variables are
%               printed, default: 'pKa [current date]'
%archiveFormat  Type of format used for returning the information. Options:
%               1 for excel file, 2 for text file, default: 1
%
%OUTPUTS:
%pKaTable       Double matrix that possess the pKa values assigned to each
%               met present in the model, in the corresponding order.
%report         boolean array, if the metabolite was found in the database,
%               one is assigned to the metabolite
%
%Authors:
%- Magdalena Ribbeck 1/18

%%
%INPUT VALIDATION

if nargin<1
    error('model is required as input');
end

if nargin<2
    c=date;
    archiveName=strcat('pKa',{' '},datestr(c));
    archiveName=archiveName{1};
end
if not(isstr(archiveName))
    error('archiveName should be a string')
end

if nargin<3
    archiveFormat=1;
end
if not(isnumeric(archiveFormat))
    error('reportFormat should be a numeric input');
elseif archiveFormat~=1 && archiveFormat~=2
    error('reportFormat should be a numeric input in the range 1 or 2');
end

%%
%ALGORITHM

[a,b,c]=xlsread('Base de datos kegg','Metabolite list');
pKaDB=c([2:end],[12:19]);
KeggIDDB=b([2:end],6);
clear a b c

keySet1=KeggIDDB;
valueSet1=[1:length(KeggIDDB)]';
dictionary=containers.Map(keySet1,valueSet1);

pKaTable=NaN(length(model.mets),length(pKaDB(1,:)));
report=zeros(length(model.mets),1);
for i=1:length(model.mets)
     if not(isempty(model.metKEGGID{i}))
        separateKeggIDs=strfind(model.metKEGGID{i}, ';');
        if isempty(separateKeggIDs)
            if isKey(dictionary,model.metKEGGID{i})
                for j=1:length(pKaDB(1,:))
                    pKaTable(i,j)=pKaDB{dictionary(model.metKEGGID{i}),j};
                end
                report(i)=1;
            end
        else
            separateKeggIDs=[0 separateKeggIDs (length(model.metKEGGID{i})+1)]; 
            for j=1:length(separateKeggIDs)-1
                sub=substrings(model.metKEGGID{i},separateKeggIDs(j)+1,separateKeggIDs(j+1)-1);
                if isKey(dictionary,sub)
                    for k=1:length(pKaDB(1,:))
                        pKaTable(i,k)=pKaDB{dictionary(sub),k};
                    end
                    report(i)=1;
                    break
                end
            end
        end
    end
end

%%
%LOG GENERATION
a=strcat(num2str(length(find(report))),{' '}, 'out of',{' '},num2str(length(model.mets)),{' '},'metabolites were found in the pKa database');

if archiveFormat==1
    xlswrite(archiveName,a,'report','A1');
    xlswrite(archiveName,{'metabolite','was the compound present in DB'}, 'report','A2');
    xlswrite(archiveName,model.mets, 'report','A3');
    xlswrite(archiveName,report, 'report','B3');
    xlswrite(archiveName,pKaTable, 'pKa','A1');
    
else
    archiveNameReport=strcat(archiveName,'_report.txt');   

    fileID=fopen(archiveNameReport,'w');
    fprintf(fileID,'%4s\n',a{1});
    fprintf(fileID,'%4s\t','metabolite');
    fprintf(fileID,'%4s\n','compound present in DB');
    for i=1:length(model.mets)
        fprintf(fileID,'%4s\t',model.mets{i});
        fprintf(fileID,'%4d\n',report(i));
    end
    fclose(fileID);

    archiveNameMain=strcat(archiveName,'.txt');
    fileID=fopen(archiveNameMain,'w');
    for i=1:length(pKaTable(:,1))
        for j=1:length(pKaTable(1,:))-1
            fprintf(fileID,'%4d\t',pKaTable(i,j));
        end
        fprintf(fileID,'%4d\n',pKaTable(i,length(pKaTable(1,:))));
    end
    
    fclose(fileID);
end

end

function sub=substrings(string,start,ending)
length=ending-start+1;

sub='';
for i=start:ending
    sub=[sub, string(i)];  
end

end