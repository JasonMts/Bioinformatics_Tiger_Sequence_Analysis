

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
a = seq(randperm(length(seq)))
len_a = length(a)

orf1 = seqshoworfs(a, 'MINIMUMLENGTH', 1, 'FRAMES', 'all', 'GENETICCODE', 2);


%find the lengths of all orfs in random sequence. make sure not to count
%the starts  that dont end
orf1Sizes =[];
for i=1:6
    if(length(orf1(i).Start) > length(orf1(i).Stop))
        orf1(i).Start = orf1(i).Start(:,1:end-1);
    end;
    orf1Sizes = [orf1Sizes, orf1(i).Stop - orf1(i).Start];

end;

%max_threshold=max(orf1Sizes)

max_threshold = prctile(orf1Sizes,99)

%convert to aminoacids so we can pass it to seqshoworfs
max_nucleotides = max_threshold/3    

%now check and find an orf 
orf = seqshoworfs(seq, 'MINIMUMLENGTH', max_nucleotides, 'FRAMES', 'all', 'GENETICCODE', 2);

possible_orf = seq(6280:7821);
possible_orf_length = length(possible_orf)


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

protein_first_50 = protein(1:50)



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






















