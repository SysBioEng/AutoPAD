function model=obtainChebiKeggIDForABiGGModel(model,archiveName, modelName)
%for a model that possess BiGG ID in model.mets, returns the KEGG ID and the
%first ChEBI ID that are given for each metabolite according to BiGG DB.
%
%USAGE:
%           model=obtainChebiKeggIDForABiGGGModel(model,archiveName, modelName)
%
%INPUTS:
%model          COBRA model structure that possess BiGG IDs in model.mets
%
%OPTIONAL INPUTS:
%archiveName    name of the archive (.mat) in which the model will be saved
%               default: 'newModel [current date]'
%modelName      name of the model struct variable for the archive,
%               default: model
%
%OUTPUTS:
%model          COBRA model structure that includes metChEBIID and metKEGGID
%
%Authors:
%- Magdalena Ribbeck 1/18

%%
if nargin<1
    error('a model is required as input'); 
end
if nargin<2
    c=date;
    archiveName=strcat('newModel',{' '},datestr(c));
    archiveName=archiveName{1};
end

if nargin<3
    modelName='model';
end


pause on


kegg=cell(length(model.mets),1);
chebi=cell(length(model.mets),1);
charges=NaN(length(model.mets),1);
for i=1:length(model.mets)
%
    pos=strfind(model.mets(i),'_');
    pos=pos{1}(end);
    met=substring(model.mets{i},1,pos-1);
    query1=strcat('http://bigg.ucsd.edu/api/v2/universal/metabolites/', met);
    try data=webread(query1);
        pause(0.5);%this pause is to avoid oversaturation of Bigg servers
        a=data.database_links.KEGGCompound.id;
        b=data.database_links.CHEBI.id;
        kegg{i}=a;
        chebi{i}=b;
        charges(i)=data.charges;
    catch
    end
    a=strcat('current state:',{' '}, num2str(i), {' '}, 'out of', {' '},num2str(length(model.mets)),{' '},'metabolites analyzed');
    disp(a);
end

model.metChEBIID=chebi;
model.metKEGGID=kegg;

if not(isfield(model,'metCharges'))
    model.metCharges=charges;
end

eval([modelName, ' =  model']);
save(archiveName,modelName);
end

function sub=substring(string,start,ending)
%small function to obtain substrings
length=ending-start+1;

sub='';
for i=start:ending
    sub=[sub, string(i)];  
end

end