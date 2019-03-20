

%1 load tiger data
gbk = getgenbank('NC_010642');


seq = gbk.Sequence;
seq_len = length(seq)
first_300_bases = seqdisp(seq(1:300))

%plot composition of bases and nucleotides

figure
ntdensity(seq);


%calculate the basecount

bases = basecount(seq);

%find orfs on both the original and the complement 

orf = seqshoworfs(seq, 'MINIMUMLENGTH', 1, 'FRAMES', 'all', 'GENETICCODE', 2);
%orf.Start

%statistician comes in:
%a = randseq(length(seq))
orf1 = []
orf1Sizes = []
maxes = []
for j=1:1000
    a = seq(randperm(length(seq)))
    
    orf1Sizes = []
    orf1 = []

    orf1 = seqshoworfs(a, 'MINIMUMLENGTH', 1, 'FRAMES', 'all', 'GENETICCODE', 2, 'nodisplay', 'true');
    
    %find the lengths of all orfs in random sequence. make sure not to count
    %the starts  that dont end
    for i=1:6
        if(length(orf1(i).Start) > length(orf1(i).Stop))
            orf1(i).Start = orf1(i).Start(:,1:end-1);
        end;
        orf1Sizes = [orf1Sizes, orf1(i).Stop - orf1(i).Start];
    end;
    
    %keep track of the max every time
    maxes = [maxes, max(orf1Sizes)];
end;

%find the lengths of all orfs in random sequence. make sure not to count
%the starts  that dont end

%max_threshold=max(orf1Sizes)

max_lengths = maxk(maxes,10)

%visualize max lengths
figure
buckets = 20
hist(maxes,buckets)
hold on;
xlabel('Max Length (Nucleotides)'); ylabel('Number of Sequences');
stem(1068,1,'k')
title('Significance of Selected ORF');
hold off;


%find 99th percentile of maxes
max_threshold = prctile(maxes,99)

%convert to aminoacids so we can pass it to seqshoworfs
max_nucleotides = max_threshold/3    

%now check and find an orf 
orf = seqshoworfs(seq, 'MINIMUMLENGTH', max_nucleotides, 'FRAMES', 'all', 'GENETICCODE', 2);



possible_orf = seq(6280:7821);
possible_orf_length = length(possible_orf)

possible_orf2 = seq(13011:14516)

possible_orf3 = seq(9563:10507)
possible_orf4 = seq(4817:5884)





%now let's calculate its p value

probT = bases.T/length(seq);
probA = bases.A/length(seq);
probG = bases.G/length(seq);
probC = bases.C/length(seq);

probTAA = probT * probA * probA
probTAG = probT * probA * probG
probTGA = probT * probA * probG

stopCodonSum = probTAA + probTAG + probTGA

%divide by 3 cause we are doing codons
p_val = (1 - stopCodonSum)^(possible_orf_length/3);

%translate to protein based on genetic code 2
protein = nt2aa(possible_orf, 'GENETICCODE', 2);
protein2 = nt2aa(possible_orf2, 'GENETICCODE', 2);
protein3 = nt2aa(possible_orf3, 'GENETICCODE', 2);
protein4 = nt2aa(possible_orf4, 'GENETICCODE', 2);


protein_first_50 = protein4(1:50)



%time to play with the human mitochondria


gbk_h = getgenbank('NC_012920');

seq_h = gbk_h.Sequence

orf_h = seqshoworfs(seq_h, 'MINIMUMLENGTH', max_nucleotides, 'FRAMES', 'all', 'GENETICCODE', 2);


%locally align COX1 with the human mitochondrial DNA

[a,b] = swalign(possible_orf, seq_h, 'ALPHABET', 'NT');

showalignment(b);

%test with randoms
for i = 1:100
    perm = randperm(length(seq_h));
    globalscores(i) = swalign(possible_orf,seq_h(perm), 'ALPHABET', 'NT');
end

%visualize histogram
buckets = 10
hist(globalscores,buckets)
hold on;
xlabel('Score'); ylabel('Number of Sequences');
hold off;

perm = randperm(length(seq_h));
[c, d] = swalign(possible_orf,seq_h(perm), 'ALPHABET', 'NT');

showalignment(d);

%now let's find global alignments between tiger and human

[scores, similar_orfs] = nwalign(seq,seq_h,'ALPHABET', 'NT');
showalignment(similar_orfs);


%find all proteins in human genome
humanproteins={};
tigerproteins={};
for i = 1:6
    %first get all proteins from human genome
    if(length(orf_h(i).Start) > length(orf_h(i).Stop))
            orf_h(i).Start = orf_h(i).Start(:,1:end-1);
    end;
      
    if(length(orf_h(i).Stop) > 0)  
        for j = 1:length(orf_h(i).Stop)
            humanproteins = [humanproteins, seq_h(orf_h(i).Start(j):(orf_h(i).Stop(j)-1))];
        end;
    end;
    
    %now let's get all proteins from tiger genome
    if(length(orf(i).Start) > length(orf(i).Stop))
            orf(i).Start = orf(i).Start(:,1:end-1);
    end;
      
    if(length(orf(i).Stop) > 0)  
        for j = 1:length(orf(i).Stop)
            tigerproteins = [tigerproteins, seq(orf(i).Start(j):(orf(i).Stop(j)-1))];
        end;
    end;
end;

humanproteins={};
tigerproteins={};
for i = 1:6
    %first get all proteins from human genome
    if(length(orf_h(i).Start) > length(orf_h(i).Stop))
            orf_h(i).Start = orf_h(i).Start(:,1:end-1);
    end;
      
    if(length(orf_h(i).Stop) > 0)  
        for j = 1:length(orf_h(i).Stop)
            humanproteins = [humanproteins, nt2aa(seq_h(orf_h(i).Start(j):(orf_h(i).Stop(j)-1)), 'GENETICCODE', 2)];
        end;
    end;
    
    %now let's get all proteins from tiger genome
    if(length(orf(i).Start) > length(orf(i).Stop))
            orf(i).Start = orf(i).Start(:,1:end-1);
    end;
      
    if(length(orf(i).Stop) > 0)  
        for j = 1:length(orf(i).Stop)
            tigerproteins = [tigerproteins, nt2aa(seq(orf(i).Start(j):(orf(i).Stop(j)-1)), 'GENETICCODE', 2)];
        end;
    end;
end;


%remember to search complement for last part
protein_scores = []
for i = 1:10
    for j = 1:10
        protein_scores = [protein_scores, nwalign(cell2mat(humanproteins(i)),cell2mat(tigerproteins(j)))];
    end;
end;

[a,b] = nwalign(cell2mat(humanproteins(2)),cell2mat(tigerproteins(10)), 'ALPHABET', 'AA')

perm = randperm(length(cell2mat(humanproteins(1))));

nwalign(cell2mat(humanproteins(1)),randseq(length(cell2mat(humanproteins(1))),'ALPHABET', 'AA'))



























