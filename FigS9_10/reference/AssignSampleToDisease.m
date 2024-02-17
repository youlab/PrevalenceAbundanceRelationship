clear;
clc;
load('KEGG_Abundance.mat');
T=readtable('disease_to_sample.txt');
SampleToDisease=cell2table(SampleID);

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='arthritis'&&string(T{i,5})=='none'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                Arthritis_None(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='arthritis'&&string(T{i,5})=='moderate'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                Arthritis_Moderate(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='arthritis'&&string(T{i,5})=='high'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                Arthritis_High(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='carcinoma'&&string(T{i,5})=='controls'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                ColorectalCancer_Controls(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='carcinoma'&&string(T{i,5})=='carcinoma'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                ColorectalCancer_Carcinoma(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='carcinoma'&&string(T{i,5})=='advanced adenoma'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                ColorectalCancer_AdvancedAdenoma(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='cirrhosis'&&string(T{i,5})=='0'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                LiverCirrhosis_None(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='cirrhosis'&&string(T{i,5})=='1'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                LiverCirrhosis(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='crohns'&&string(T{i,5})=='0'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                CrohnsDisease_None(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='crohns'&&string(T{i,5})=='1'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                CrohnsDisease(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='obesity'&&string(T{i,5})=='0'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                Obesity_None(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='obesity'&&string(T{i,5})=='1'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                Obesity(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='t2d'&&string(T{i,5})=='0'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                Type2Diadetes_None(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='t2d'&&string(T{i,5})=='1'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                Type2Diadetes(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='ulcerative_colitis'&&string(T{i,5})=='0'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                UlcerativeColitis_None(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end

pin=1;
for i=1:size(T,1)
    i
    if string(T{i,4})=='ulcerative_colitis'&&string(T{i,5})=='1'
        for j=1:length(SampleID)
            if string(SampleID(j))==string(T{i,2})
                UlcerativeColitis(:,pin)=GeneAbundanceTable(:,j);
                pin=pin+1;
                break;
            end
        end
    end
end


