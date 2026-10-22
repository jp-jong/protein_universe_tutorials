# Profile-based homology searching and SLURM commands
#### Authors: Joenilo II Paduhilao and Ursula Oggenfuss

## Introduction
We use hidden Markov models, or profile HMMs, to gather sequences that are homologous to our proteins of interest. Since profile HMMs model evolutionary changes given an alignment, it is more sensitive than using one sequence alone in gathering homologues. However, databases, like PFAM, only include a subset of sequences (usually from protein databases like UniProt) to create their alignment for their profile HMMs. Thus, we are unsure whether it is sensitive to gather all candidate homologues against a database that is not well-represented in these protein databases, in particular Asgard archaea. It is to our best of interest to refine our profile HMMs in N iterations to make it more sensitive. In this way, we create our own profile HMMs tailored for own databases. 

The flow chart below is a typical workflow in gathering homologues using profile HMMs, validate your hits, and how to execute iterations to make sensitive profile HMMs.

<image>

In this tutorial, we will use the HMMER suite and other tools (mentioned below) to address the following learning objectives.

### Objectives:
1. To gather homologues given a profile HMM;
2. To validate homologue hits; and
3. To create a custom profile HMM.

### Software and versions
1. HMMER suite
2. Seqtk
3. Seqkit
4. MAFFT
5. Interproscan
6. HHpred  https://toolkit.tuebingen.mpg.de/tools/HHpred
7. FoldSeek https://search.foldseek.com/search


## Tutorial
1. One of the most important steps in homology searches is to decide on the starting material, in this case the profile HMM. There are many databases for which we can obtain HMM files (e.g., PFAM, PantherDB) but for this tutorial, we will obtain PF01918 from PFAM. PF01918 is the PFAM profile HMM for the Alba superfamily. Go to this link https://www.ebi.ac.uk/interpro/entry/pfam/PF01918/ and downoad the raw HMM (.hmm) for this family (Profile HMM > Download). Place this in an appropriate directory, or follow the codes below. 

```bash
mkdir tutorial_1_homology_slurm

mv <path_of_your_hmm> tutorial_1_homology_slurm
```
Obviously, profile HMMs are for machines to read, not for you. But there are many tools out there that parses .hmm files to make it more interpretable for us, such as https://skylign.org/ and the hmmemit function from the HMM suite. For brevity, it is not included in this part of the tutorial. 

Note: It is always recommended to read the overview page because it gives you seminal papers regarding that specific protein family. For example, it was Aravind et al., (2003) that unified the Alba superfamily, including its eukaryotic homologues. You can see that this reference was also cited in PFAM.

2. We will use the raw profile HMM to do our first homology search against a toy dataset. Asgard56 is a dataset composed of, obviously, 56 MAGs/genomes of Asgard archaea. They were selected based on the best genome quality scores per Asgard archaea family. A protein database (usually in .faa or .fa format) is a concatenation of ALL proteins these MAGs/genomes encode. Therefore, it will take so much time and effort to go one by one and blast them to obtain your protein homologues. Even genome annotations are not 100% reliable for this task. 

```bash
# declare your HMM variables and database path
HMM=<path_to_hmms>
ASGARD_DB=<path_to_asgard_db>

# create directory for your hmm hits
mkdir hmm_hits

# for demonstration, we compare outputs of hmmsearch with and without --max flag. This flag controls the search filters of the program and we will see the difference later on

hmmsearch --cpu 30 --tblout hmm_hits/PF01918_Asgard56.tblout ${HMM} ${ASGARD_DB}

hmmsearch --max --cpu 30 --tblout hmm_hits/PF01918_Asgard56.max.tblout ${HMM} ${ASGARD_DB} 

# make tsv files to make reading easier 
## you can do this one-by-one 
grep -v "^#" hmm_hits/PF01918_Asgard56.tblout > hmm_hits/PF01918_Asgard56.tsv
gawk -i inplace '{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' hmm_hits/PF01918_Asgard56.tsv
sed -i '1i target_name\taccession\tquery_name\taccession\tE-value\tscore\tbias\td_E-value\td_score\td_bias' hmm_hits/PF01918_Asgard56.tsv


## or using a for loop
for tblout in hmm_hits/*.tblout; do
name=$(basename $tblout .tblout)
grep -v "^#" $tblout > hmm_hits/${name}.tsv
gawk -i inplace '{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' hmm_hits/${name}.tsv
sed -i '1i target_name\taccession\tquery_name\taccession\tE-value\tscore\tbias\td_E-value\td_score\td_bias' hmm_hits/${name}.tsv
done
```

Explore the TSV files. Take note what each column represents. You can also see the homologue hits and their names. By now, you should see the name "hypothetical protein". Try to look up what they are. 

hmmsearch reports two kinds of E-value and scores. Sequence-level statistics usually refer to the statistical significance of an entire protein sequence matching the model while domain-level statistics refer to the statistical signficance of a single domain within that sequence. If your protein of interest is composed of multiple-domains, it is usually recommended to look at both metrics. However, if you have only one protein domain in question, in this case the Alba domain, deciding on the significance threshold based on the domain-level seem more appropriate. In essence, e-values that are < 1e-5 are tagged as confident homologue hits. We can also be more inclusive by increasing this value to 1e-3. Take note of the domain bit scores for which these e-value thresholds are reported. They are important in exploring the twilight zone hits. 

Note: On your free time, look up the difference between e-value and domain bit scores. This is equally important if you are also using BLAST.

3. The twilight zone hits are those sequences with 25-30% similarity with that of our protein queries and are usually hits on the e-value/bit-scores thresholds borderline. We do not know if they are true homologues or not. Thus, we need to validate them. Because structures are more converseved than sequences, we can rely more on them to see whether these hits are indeed candidate homologues. We validate structural homology using HHpred and FoldSeek webservers. Importantly, we can also use InterProScan and other annotation tools to detect the domains we are interested in.

```bash
# first, we extract all sequences that were hit by this profile hmm. Use the output with the --max flag
mkdir sequences

# extract the --max hits
seqtk subseq ${ASGARD_DB} <(awk -F "\t" '{print $1}' hmm_hits/PF01918_Asgard56.max.tsv | sort -u) > sequences/PF01918_Asgard56.max.faa

# count how many you got and record this
grep -c ">" sequences/PF01918_Asgard56.max.faa

# InterProScan
# annotate with interpro
mkdir annotations
mkdir -p annotations/temp

INTERPRO="/local/two/jong2x/softwares/interproscan-5.77/interproscan-5.77-108.0"

sed 's/*//g' sequences/PF01918_Asgard56.max.faa | ${INTERPRO}/./interproscan.sh -appl Pfam-38.1 -f tsv -i - -b annotations/PF01918_Asgard56.max.pfam --tempdir annotations/temp
```

Check whether our sequences hit PF01918 or Alba domains. Since we included all of the homologues, we also obviously annotated the false positive ones. Observe the annotations. Are they also annotated as PF01918? Another important thing to note here is to observe whether you hit other domains for those that we are confident with. Usually, they are protein sequences with multiple domains, coiled coils, etc. 

Go back to our TSV file and explore 5 hits up and down of the threshold 1e^-3. Query these proteins sequences (from .faa file) to HHpred https://toolkit.tuebingen.mpg.de/tools/HHpred and/or FoldSeek https://search.foldseek.com/search . For HHpred, use the following databases:
* PDB_mmCIF70_20_Feb - protein structure database
* Pfam-A_v38.2 - sequence-based domain database
* COG_KOG_v1.0 - sequence-based orthologues database
* NCBI_Conserved_Domains(CD)_v3.19 - sequence-based domain database

Rule of thumb:
* If HHpred predicts a protein sequence more than 95% probability with an e-value less than 1e-3, I would trust them as candidate homologues and include them in our list.
* If HHpred predicts a protein sequence with more than 95% probability with slightly high e-values (more than 1e-3), I would query this protein in FoldSeek (predict > Structure with ESMFold > search). If it hits the protein that we are interested in, with TM-scores > 0.5 and RMSD < 5, I would trust the protein sequence and include it in my list of homologues. 
* If none of the HHpred rules is satisfied above, use FoldSeek and check if we hit the protein that we are interested in. Usually, these sequences are highly diverged so they might be problematic in the alignment or in building phylogenetic trees. However, at this point, we want to be as inclusive as possible and exclude them in downstream analyses after careful decisions. 
* If InterProScan, HHpred, and FoldSeek does not give us any homology signal, we do not trust it as homologous to our protein of interest.

At this point, keep on validating even deeper into the twilight zone hits. Once you pick up proteins lacking homology signal, record the last domain bit score (not the e-value). This becomes your cut-off for the first iteration. 

Note: Record as extensive as you can. Use excel if that is more convenient for you. Use word document if you like screenshots. This will help the future you on why you made certain decisions or why those results actually show up. Results might not make sense now but sooner or later, it will.

4. We then filter our protein sequences by repeating the seqtk process in no.3. We do this by specifying the domain bit score threshold on the 9th column, not 8th (for e-value). Re-extract the sequences and count how many you got compared to the --max output.
```bash
seqtk subseq ${ASGARD_DB} <(awk -F "\t" '$9 > threshold?' hmm_hits/PF01918_Asgard56.max.tsv | sort -u) > sequences/PF01918_Asgard56.d<threshold>.faa

# count how many passed your cut-off threshold
grep -c ">" sequences/PF01918_Asgard56.d<threshold>.faa
```

It is at this point that you also reflect whether there were protein sequences tagged as "confident hit" but you do not fully trust it because it had ambiguous domain hits from InterProScan or other predicted structures. Also, the next step of this tutorial is to do the 2nd iteration. What should you do with the sequences with other domain hits? It is good to think about it at this point. 

5. Our goal is to make a sensitive profile HMM tailed to our research objectives. Now that we sort of added homologue hits in the twilight zone, we also captured additional proteins that probably have more diverged sequences. We can model evolutionary changes given this set of candidate homologues by aligning the sequences and building our own profile HMM. This time around, you know that it will be an Asgard-specific profile HMM for obvious reasons. Although it does not mean it is not sensitive to other organisms because homology is a transitive property.

```bash
# create directory for alignments
mkdir alignments

# we use mafft-linsi since it's faily accurate and used extensively in our field

mafft-linsi --thread 16 sequences/PF01918_Asgard56.d<threshold>.faa > alignments/PF01918_Asgard56.d<threshold>.faa

# trimming it is not necessary as hmmbuild will search up the alignment and detect informative regions 
hmmbuild --amino hmms/PF01918_Asgard56.hmm alignments/PF01918_Asgard56.d<threshold>.faa

# investigate the new profile hmm in skylign and compare this from PF01918.hmm 
# if you have more than 1000 sequences, it is best to sample out 1000 because skylign has its own number of sequences limitations 
# for this tutorial, that is not needed
hmmemit -N 1000 hmms/PF01918_Asgard56.hmm > hmms/PF01918_Asgard56.hmm.sampled.faa

# it is also helpful to compare the lengths of your profile hmm using hmmstats
hmmstat hmms/PF01918_Asgard56.hmm 
hmmstat hmms/PF01918.hmm
```

6. Now, we do our 2nd iteration of homology search. Repeat steps 2-5 and evaluate whether you gather new homologues or if they are the same homologues. The steps below summarize what has been done previously in the most straightfoward manner. 
``` bash
# declare your HMM variables and database path
HMM=<path_to_new_hmms>
ASGARD_DB=<path_to_asgard_db>
N=<iteration_number>

# do the hmmersearch with --max

hmmsearch --max --cpu 30 --tblout hmm_hits/homs_Asgard56_iter.${N}.max.tblout ${HMM} ${ASGARD_DB} 

# make tsv files to make reading easier 
## you can do this one-by-one  or use the for loop above
grep -v "^#" hmm_hits/homs_Asgard56_iter.${N}.max.tblout > hmm_hits/homs_Asgard56_iter.${N}.max.tsv
gawk -i inplace '{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' 
hmm_hits/homs_Asgard56_iter.${N}.max.tsv
sed -i '1i target_name\taccession\tquery_name\taccession\tE-value\tscore\tbias\td_E-value\td_score\td_bias' hmm_hits/homs_Asgard56_iter.${N}.max.tsv

# extract the --max hits
seqtk subseq ${ASGARD_DB} <(awk -F "\t" '{print $1}' hmm_hits/PF01918_Asgard56.max.tsv | sort -u) > sequences/homs_Asgard56_iter.${N}.max.faa

# count how many you got and record this
grep -c ">" sequences/homs_Asgard56_iter.${N}.max.faa

# annotate with interpro
INTERPRO="/local/two/jong2x/softwares/interproscan-5.77/interproscan-5.77-108.0"

sed 's/*//g' sequences/homs_Asgard56_iter.${N}.max.faa | ${INTERPRO}/./interproscan.sh -appl Pfam-38.1 -f tsv -i - -b annotations/homs_Asgard56_iter.${N}.max.pfam --tempdir annotations/temp
```
 
Explore the twilight zone up/down of the borderline. Again, include 5-10 genes up the threshold because you might included hits that are not really trustable. This is possible as you added more sequences into your profile hmm. Equally important, explore more genes down the threshold if you still gather putative homologues. Validate with HHpred and FoldSeek. Be patient with this and do not rush. Lastly, record extensively from the get go.

7. Once you decide on a new threshold, record that domain bit score and re-build a profile HMM. 
```bash
seqtk subseq ${ASGARD_DB} <(awk -F "\t" '$9 > threshold?' hmm_hits/homs_Asgard56_iter.${N}.max.tsv | sort -u) > sequences/homs_Asgard56_iter.{$N}.d<threshold?>.faa

# count how many passed your cut-off threshold
grep -c ">" sequences/homs_Asgard56_iter.{$N}.d<threshold?>.faa

# align
mafft-linsi --thread 16 sequences/homs_Asgard56_iter.{$N}.d<threshold?>.faa > alignments/homs_Asgard56_iter.{$N}.d<threshold?>.linsi

# trimming it is not necessary as hmmbuild will search up the alignment and detect informative regions 
hmmbuild --amino hmms/homs_Asgard56_iter.${N}.d<threshold?>.hmm alignments/homs_Asgard56_iter.{$N}.d<threshold?>.linsi 

# see the stats if you think that is necessary

# repeat homology search again with this new profile hmm
```

8. Sometimes, you do not need multiple iterations. Sometimes, you need multiple iterations. It really depends on the nature of the protein family and also the database that you have. But for everytime you do an iteration, check if the homologues you have gathered is the same as the previous iteration. If they are, it converged! You now have a sensitive enough custom profile HMM. If you have done everything as faithful and as thorough as possible, you are quite assured you have gatherd a highly confident protein homologue set. Congrats. 

Note: Even profile HMMs can miss highly divergent sequences beyond the twilight zone hits. See papers like Sahakyan et al., (2025) and the midnight zone. This zone can be safely explored using Structural homology searches. But it can be computationally expensive and thus, it is not really recommended if sequence-based or profile-based searches can already answer your research questions. 

Another note: This tutorial is more on the foundations. More advanced means to do iterations can be done, for instance, extracting the domain subsequences from hits with multiple domain hits or ambiguous hits. The goal is to cleanse and start with "clean" sequences, as much as possible. Below is an example of how I do it.

```bash
# to extract sequences based on region specified in bed file:
# bed file format: <sequence_name> <start> <end>
# coordinates can be based on PFAM or any annotations like from HHpred

seqtk subseq <.faa file> <bed file> > <name>.extract_seq.faa

# to remove sequences - create exclude list
# format: <sequence_name_to_exclude_only>
seqkit grep -v -f <remove_seq.txt> <.faa file> > <name>.removed.faa

# you can then merge these sequences because <name>.extract_seq.faa only includes sequences that requires extraction
cat <name>.extract_seq.faa <name>.removed.faa > <name>.cleaned.faa

# then proceed with alignment, hmmbuild, hmmsearch, etc.
```

Once you finish collecting the homologue set, you can proceed to clustering and characterizing them (e.g., motifs, biochemical parameters, etc.). But before then, gather a good protein homologue set first. 
