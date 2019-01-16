% % % 
% % % 
% % % Trabajo Fin de Grado de Juan Manuel García Ilarramendi
% % % Contacto: a904429@alumni.unav.es
% % % 

%% Obtener molecular signatures 
filepath=('D:\PFG\Recon 2\Mol files\Signatures'); %% Introducir directorio donde están los archivos sscan
cd(filepath)
dinfo=dir(filepath);

names_files={dinfo.name};
names_files=(names_files(1,3:end));
for j=1:length(names_files)
    fileID=fopen(names_files{1,j});
    a=erase(names_files{1,j},'.scan3');
    B=textscan(fileID,'%s %s');
    fclose(fileID)
    eval(sprintf('Signatures.%s= horzcat(B{1,:});',a))
end    
%% Obtener reaction signatures 

load Recon2.v04.mat %% Elegir base de datos donde se extraerán coeficientes estequiométricos de las reacciones 
load 'Mol Sig.mat' %% Archivo donde está guardado el struct con la smolecular signatures 

Rxn=full(modelR204.S);
Names=modelR204.mets;

[nf,nc]=size(Rxn);
cont=0;
reac_non={};
for i=1:nf
    Names{i,1}=eraseBetween(Names{i,1},'[',']');
    Names{i,1}=erase(Names{i,1},'[]');
end 
Names1=unique(Names,'stable');

for i=1:nc
    reac=Rxn(:,i);
    names=Names(reac~=0);
    coef=reac(reac~=0);
    
    
   if all(ismember(names,Signatures.ref(:,2)))
           ref=names;
          signa={};
          
         for i=1:length(names) 
            eval(sprintf('signa{j,1}=Signatures.Sig_%s;',ref{i,1})); 
            signa{i,2}=coef(i);
            [nf,nc]=size(signa{i,1});
            for k=1:nf
                signa{i,1}{k,1}=signa{i,1}{k,1}*signa{i,2};
            end   
            
         end
         signa(:,2)=[];
         signa=vertcat(signa{:});
        uniques=unique(signa(:,2),'stable');
        for i=1:length(uniques)
            cof=signa(strcmp(uniques{i},signa(:,2))==1,1);
            uniques{i,2}=sum(cell2mat(cof));
        end
        uniques=uniques(~cell2mat(uniques(:,2))==0,:);
         uniques(strcmp(uniques(:,1),''),:)=[];
         eval(sprintf('RxnSig.Rxn_%d=uniques;',i))
         
   else
       cont=cont+1;
       reac_non{cont}=i;
   end    
                
end   
fprintf('%d reacciones no incluidas\n',cont)
A=fieldnames(RxnSig);
for i=1:length(fieldnames(RxnSig))
    if eval(sprintf('isempty(RxnSig.%s)',A{i}))
        eval(sprintf('RxnSig=rmfield(RxnSig,''%s'');',A{i}))
    end    
end   

%%  Obtener todos las reacciones donde haya atomic signatures de poliaminas en los substratos 




load ('Mol Sig.mat')
load('Reaction Sig.mat')
load Recon2.v04.mat


Putrescine=Signatures.Sig_ptrc;
Spermidine=Signatures.Sig_spmd;
Spermine=Signatures.Sig_sprm;




cont=0;

SUBS=cell(0,1);

files=fieldnames(RxnSig);
poly={'Putrescine';'Spermidine';'Spermine'};

for k=1:length(poly)
    cont=0;
    eval(sprintf('polyamine=%s;',poly{k}))
    for j=1:length(files)
        eval(sprintf('rxn=RxnSig.%s(cellfun(@(x) x<=0,RxnSig.%s(:,2)),1);',files{j},files{j}))

            for i=1:length(polyamine)
                if ismember(polyamine{i,2},rxn)==1
                    cont=cont+1;
                    eval(sprintf('SUBS_%d{cont,1}=''%s'';',k,files{j}))
                end
            end
    end  
end   

PUTR_subs=unique(SUBS_1,'stable');
PUTR_subs(:,2)=cellfun(@(x) sum(ismember(SUBS_1,x)),PUTR_subs(:,1),'un',0);

SPER_subs=unique(SUBS_2,'stable');
SPER_subs(:,2)=cellfun(@(x) sum(ismember(SUBS_2,x)),SPER_subs(:,1),'un',0);

SPER2_subs=unique(SUBS_3,'stable');
SPER2_subs(:,2)=cellfun(@(x) sum(ismember(SUBS_3,x)),SPER2_subs(:,1),'un',0);

[~, order] = sort(cell2mat(PUTR_subs(:, 2)),'descend');
PUTR_subs= PUTR_subs(order, :);

[~, order] = sort(cell2mat(SPER_subs(:, 2)),'descend');
SPER_subs= SPER_subs(order, :);

[~, order] = sort(cell2mat(SPER2_subs(:, 2)),'descend');
SPER2_subs= SPER2_subs(order, :);

load Recon2.v04.mat

modelR204.mets = cellfun(@(S) S(1:end-3), modelR204.mets, 'Uniform', 0);

Rxn=full(modelR204.S);

files={'PUTR_subs';'SPER_subs';'SPER2_subs'};

for m=1:length(files)
    eval(sprintf('polyamine=%s;',files{m}))

    for j=1:length(polyamine)
        polyamine{j,3}=modelR204.mets(Rxn(:,str2double(erase(polyamine{j,1},'Rxn_')))<0);

        for i=1:length(polyamine{j,3})
            eval(sprintf('sig=Signatures.Sig_%s;',polyamine{j,3}{i,1}))
            eval(sprintf('kernel=intersect(%s(:,2),sig(:,2));',poly{m}))
            if isempty(kernel) 
                   polyamine{j,3}{i,1}={};
            else   
                for k=1:length(kernel)
                    eval(sprintf('kernel{k,2}=%s{strcmp(kernel{k,1},%s(:,2)),1}*sig{strcmp(kernel{k,1},sig(:,2)),1};',poly{m},poly{m}))
                end

                eval(sprintf('polyamine{j,3}{i,2}=sum(cell2mat(kernel(:,2)))/(norm(cell2mat(%s(:,1)))*norm(cell2mat(sig(:,1))));',poly{m}))
            end 

        end   
    end 
    eval(sprintf('%s=polyamine;',files{m}))
end    

% %% Polyamine as product
% 
% for i=1:length(files)
%     eval(sprintf('rxn=RxnSig.%s(cellfun(@(x) x>0,RxnSig.%s(:,2)),1);',files{i},files{i}))
%     
%         for i=1:length(Putrescine)
%             if ismember(Putrescine{i,2},rxn)==1
%                 cont=cont+1;
%                 eval(sprintf('PROD{cont,1}=''%s'';',files{i}))
%             end
%         end
% end       
% 
% PUTR_prod=unique(PROD,'stable');
% PUTR_prod(:,2)=cellfun(@(x) sum(ismember(PROD,x)),PUTR_prod(:,1),'un',0);
% 
% [~, order] = sort(cell2mat(PUTR_prod(:, 2)),'descend');
% PUTR_prod= PUTR_prod(order, :);
% 
% 
% 
% 
% for i=1:length(PUTR_prod)
%     PUTR_prod{i,3}=modelR204.mets(Rxn(:,str2double(erase(PUTR_prod{i,1},'Rxn_')))<0);
% 
%     for i=1:length(PUTR_prod{i,3})
%         eval(sprintf('sig=Signatures.Sig_%s;',PUTR_prod{i,3}{i,1}))
%         kernel=intersect(Putrescine(:,2),sig(:,2));
%         if isempty(kernel) 
%             PUTR_prod{i,3}{i,1}={};
%         else   
%             for k=1:length(kernel)
%                 kernel{k,2}=Putrescine{strcmp(kernel{k,1},Putrescine(:,2)),1}*sig{strcmp(kernel{k,1},sig(:,2)),1};
%             end
% 
%             PUTR_prod{i,3}{i,2}=sum(cell2mat(kernel(:,2)))/(norm(cell2mat(Putrescine(:,1)))*norm(cell2mat(sig(:,1))));
%         end 
%         if PUTR_prod{i,3}{i,2}==1
%             PUTR_prod{i,3}{i,1}={};
%             PUTR_prod{i,3}{i,2}={};
%         end    
%     end   
% end    


%% Comprobar parte que reacciona en substrato coincide con la de la poliamina

load 'Mol Sig.mat'
load 'Reaction Sig.mat'

Putrescine=Signatures.Sig_ptrc;
Spermidine=Signatures.Sig_spmd;
Spermine=Signatures.Sig_sprm;

files={'PUTR_subs','Putrescine','ptrc';'SPER_subs','Spermidine','spmd';'SPER2_subs','Spermine','sprm'};

for m=1:length(files)
    cont=0;
    eval(sprintf('polyamine=%s;',files{m,1}))
    polyamine_borrar={};
    for i=1:length(polyamine)
        Mirar=polyamine{i,3};
        if length(Mirar(cell2mat(Mirar(:,2))>0,1))>1
            for k=1:length(Mirar(cell2mat(Mirar(:,2))>0,1))
                Caso_esp=Mirar(cell2mat(Mirar(:,2))>0,1);
                eval(sprintf('subs=Signatures.Sig_%s;',Caso_esp{k}))
                eval(sprintf('rxn=RxnSig.%s;',polyamine{i,1}))
                rxn=rxn(cellfun(@(x) x<0,rxn(:,2)),1);
                int1=intersect(rxn,subs(:,2),'stable');
                eval(sprintf('int2=intersect(rxn,%s(:,2),''stable'');',files{m,2}))
                if all(ismember(int1(:),int2(:))) && length(int1)==length(int2) && eval(sprintf('any(strcmp(''%s'',%s{i,3}(:,1)))==0',files{m,3},files{m,1}))
            
                    cont=cont+1;
                    eval(sprintf('polyamine_borrar(cont,1:3)=%s(i,:);',files{m,1}))
                    eval(sprintf('polyamine_borrar{cont,end+1}=''%s replaced'';',Caso_esp{k}))
                end
            end
        else   
                
        eval(sprintf('subs=Signatures.Sig_%s;',Mirar{cell2mat(Mirar(:,2))>0,1}))
        eval(sprintf('rxn=RxnSig.%s;',polyamine{i,1}))
        rxn=rxn(cellfun(@(x) x<0,rxn(:,2)),1);

        int1=intersect(rxn,subs(:,2),'stable');
        eval(sprintf('int2=intersect(rxn,%s(:,2),''stable'');',files{m,2}))   

        if all(ismember(int1(:),int2(:))) && length(int1)==length(int2) && eval(sprintf('any(strcmp(''%s'',%s{i,3}(:,1)))==0',files{m,3},files{m,1}))
            
            cont=cont+1;
            eval(sprintf('polyamine_borrar(cont,1:3)=%s(i,:);',files{m,1}))
        end
        end
    
    end
    
  eval(sprintf('%s=polyamine_borrar;',files{m,1}))  
end   

%% Obtener productos

clc


load 'Mol Sig.mat'
load 'Reaction Sig.mat'
load Recon2.v04.mat

modelR204.mets = cellfun(@(S) S(1:end-3), modelR204.mets, 'Uniform', 0);
Rxn=full(modelR204.S);
Putrescine=Signatures.Sig_ptrc;
Spermidine=Signatures.Sig_spmd;
Spermine=Signatures.Sig_sprm;


files={'PUTR_subs';'SPER_subs';'SPER2_subs'};
poly={'Putrescine';'Spermidine';'Spermine'};


for m=1:length(files)
    eval(sprintf('polyamine=%s;',files{m}))
    [nf1,nc1]=size(polyamine);

    for k=1:nf1

        Ejemplo=polyamine{k,3};
        Ejemplo=Ejemplo(cell2mat(Ejemplo(:,2))==0,:);
        if isempty(Ejemplo) && min(size(polyamine{k,3}))==2
            polyamine{k,7}='CASO ESPECIAL';
            [~,ind]=min(cell2mat(polyamine{k,3}(:,2)));
            Ejemplo=polyamine{k,3}(ind,:);
        end    
        [nf,nc]=size(Ejemplo);
        eval(sprintf('rxnsig=RxnSig.%s;',polyamine{k,1}))
        subs={};
        cont=0;
        for i=1:nf
                cont=cont+1;
                eval(sprintf('subs{cont}=Signatures.Sig_%s;',Ejemplo{i,1}))
            
        end 


        eval(sprintf('subs{cont+1,1}=%s;',poly{m}))
        subs=vertcat(subs{:});
        [~,ind]=unique(subs(:,2),'stable');
        subs1=subs(ind,:);
        for i=1:length(subs1)
            subs1{i,1}=sum(cell2mat(subs(strcmp(subs1{i,2},subs(:,2)),1)));
        end    
        prod={subs1(:,2);rxnsig(:,1)};
        prod=vertcat(prod{:});
        prod=unique(prod,'stable');

        for i=1:length(prod)
            if ismember(prod(i),subs1(:,2))==1 && ismember(prod{i},rxnsig(:,1))==1
                cof=cell2mat(subs1(strcmp(prod{i},subs1(:,2)),1))+ cell2mat(rxnsig(strcmp(prod{i},rxnsig(:,1)),2));
            else
                if ismember(prod{i},subs1(:,2))==1
                    cof= cell2mat(subs1(strcmp(prod{i},subs1(:,2)),1));
                end
                if ismember(prod{i},rxnsig(:,1))==1
                    cof= cell2mat(rxnsig(strcmp(prod{i},rxnsig(:,1)),2));
                end
            end
            prod{i,2}=cof;
        end   
        prod=prod(~cellfun(@(x) x<=0, prod(:,2)),:);
        polyamine{k,4}=fliplr(prod);
    end   
    eval(sprintf('%s=polyamine;',files{m}))   
end    

%% 
 
clc


load Recon2.v04.mat
load 'Mol Sig.mat'


files={'PUTR_subs','Putrescine';'SPER_subs','Spermidine';'SPER2_subs','Spermine'};
Rxn=full(modelR204.S);
modelR204.mets = cellfun(@(S) S(1:end-3), modelR204.mets, 'Uniform', 0);

Putrescine=Signatures.Sig_ptrc;
Spermidine=Signatures.Sig_spmd;
Spermine=Signatures.Sig_sprm;

for m=1:length(files)
    eval(sprintf('polyamine=%s;',files{m,1}))
    [nf1,nc1]=size(polyamine);

        for j=1:nf1
            if isempty(polyamine{j,4})==0  

                prod=modelR204.mets(Rxn(:,str2double(erase(polyamine{j,1},'Rxn_')))>0);
                    if isempty(prod)
                        continue
                    else
                        
                        for i=1:length(prod)
                            eval(sprintf('polyamine{j,5}{i,1}=Signatures.Sig_%s;',prod{i}))
                        end

                    polyamine{j,5}=vertcat(polyamine{j,5}{:});
                    [~,ind]=unique(polyamine{j,5}(:,2),'stable');
                    polyamine{j,5}=polyamine{j,5}(ind,:);
                    cont1=0;
                    
                    for k=1:length(polyamine{j,4})
                        for l=1:length(polyamine{j,5})
                            if strcmp(polyamine{j,4}{k,2},polyamine{j,5}{l,2}) && polyamine{j,4}{k,1}~=polyamine{j,5}{l,1} 
                                cont1=cont1+1;
                                polyamine{j,6}{cont1,1}=polyamine{j,4}{k,1};
                                
                                polyamine{j,6}{cont1,2}=polyamine{j,4}{k,2}
      
                            end
                            
                         end
                        
                        if ismember(polyamine{j,4}{k,2},polyamine{j,5}(:,2))==0
                            cont1=cont1+1;
                            polyamine{j,6}(cont1,:)=polyamine{j,4}(k,:);
                        end      
                    end   
                            
                    
                    end
                    
            else
                continue
            end    
        end
        polyamine(cellfun('isempty',polyamine(:,5)),:)=[];
        
        eval(sprintf('%s=polyamine;',files{m,1}))
%             [nf,nc]=size(polyamine{j,3});
%                 for k=1:nf
%                     
%                     if polyamine{j,3}{k,2}~=0 && polyamine{j,3}{k,2}~=1
%                         eval(sprintf('analyze=Signatures.Sig_%s;',polyamine{j,3}{k,1}))
%                         eval(sprintf('res=setdiff(%s(:,2),analyze(:,2));',files{m,2}))
%                     end
%                 end    
end              

% HACER RESTA DE COLUMNAS 5 Y 4 
% Distinguir entre aquellas moleculas formadas por incluir poliamina 
clc
clear all


load Products.mat
load 'Mol Sig.mat'
load 'Reaction Sig.mat'
load 'Rec_3 Sig.mat'

files={'PUTR_subs','Putrescine';'SPER_subs','Spermidine';'SPER2_subs','Spermine'};


Putrescine=Signatures.Sig_ptrc;
Spermidine=Signatures.Sig_spmd;
Spermine=Signatures.Sig_sprm;

 

% mol={};
% for i=1:length(SPER2_subs{1,3})
%       eval(sprintf('mol{i,1}=Signatures.Sig_%s;',SPER2_subs{1,3}{i,1}))
% end    
% 
% mol=vertcat(mol{:});
% 
% [~,ind]=unique(mol(:,2),'stable');
%         subs1=mol(ind,:);
%         for i=1:length(subs1)
%             subs1{i,1}=sum(cell2mat(mol(strcmp(subs1{i,2},mol(:,2)),1)));
%             if any(strcmp(subs1{i,2},B(:,1)))
%                 subs1{i,1}=subs1{i,1}+B{strcmp(subs1{i,2},B(:,1)),2};
%             end
%         end    
%  
% subs1=vertcat(subs1,fliplr(B(cellfun(@(x) x>0, B(:,2)),:)));
% subs1(cellfun(@(x) x==0, subs1(:,1)),:)=[];   


a=fieldnames(SignaturesRecon3);
a=a(startsWith(a,'Sig_'));
cont=0;
buscar=Spermidine(1:end-1,2);

for i=1:length(a)
    eval(sprintf('mole=SignaturesRecon3.%s(1:end-1,2);',a{i}))
    
    if any(ismember(buscar,mole)) 
        cont=cont+1;
        eval(sprintf('posible{cont,1}=SignaturesRecon3.%s(1:end-1,:);',a{i}))
        posible{cont,2}=erase(a{i},'Sig_');
    end
end    



for m=1:length(files)
    eval(sprintf('polyamine=%s;',files{m,1}))
    
    [nf,nc]=size(polyamine);
    cont=0;
    for i=1:nf
        
        comp=polyamine{i,4};
        
        
        for j=1:length(posible)
            
            if all(ismember(posible{j,1}(:,2),comp(:,2))) 
                cont=cont+1;
                eval(sprintf('producto_%s{cont,1}=''%s'';',files{m,2},polyamine{i,1}))
                eval(sprintf('producto_%s{cont,2}=posible{j,1};',files{m,2}))
                eval(sprintf('producto_%s{cont,3}=posible{j,2};',files{m,2}))
            end
        end
    end
end    




