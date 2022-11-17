import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

alphabet = 'ACGT'

def empirical_frequencies(counts):
    '''
    convert a nucleotide counts (dictionary {'A':3423, 'C':34543, ...}) to frequency
    '''
    total = sum(list(counts.values()))
    return np.array([counts[n]/total for n in alphabet])


def spectrum_to_matrix(spec):
    '''
    convert dictionary of mutation counts to mutation matrix
    '''
    M = np.zeros((4,4))
    for i1,n1 in enumerate(alphabet):
        for i2,n2 in enumerate(alphabet):
            if n1!=n2:
                M[i2,i1] = spec[f"{n1}to{n2}"]
    # normalize off-diagonal rates (just for standardization, doesn't affect the results)
    M /= M.sum()
    # will the diagonal with 'outflow' term to guarantee conservation of probability
    d = M.sum(axis=0)
    np.fill_diagonal(M,-d)

    return M


def equilibrium_probabilities(M):
    evals, evecs = np.linalg.eig(M)
    # find zero eigenvalue
    ii = np.argmin(np.abs(evals))
    assert np.abs(evals[ii])<1e15
    print(evecs)
    # pull out corresponding eigenvector, return normalized to sum_i p_i = 1
    p = evecs[:,ii]
    return p/p.sum()

def get_ffsyn(ref):
    from Bio import SeqIO
    features  = [{'pos':i, 'nuc':n} for i,n in enumerate(ref.seq)]
    ## make a map of codon positions, ignore orf9b (overlaps N)
    for feat in ref.features:
        if feat.type=='CDS':
            gene_name = feat.qualifiers.get('gene', ['-'])[0]
            if gene_name!='ORF9b':
                for gpos, pos in enumerate(feat.location):
                    tmp = {'pos_in_codon':(gpos%3)+1, 'gene':gene_name, 'four_fold':False}
                    if (gpos%3)+1==3:
                        codon12 = "".join(ref.seq[pos-2:pos])
                        if codon12 in ['TC','CT', 'CC', 'CG', 'AC', 'GT', 'GC','GG']:
                            tmp['four_fold'] = True
                    features[pos].update(tmp)

    return pd.DataFrame(features)

def get_sequence_from_genbank(accession):
    from Bio import SeqIO
    from Bio import Entrez
    Entrez.email = "richard.neher@unibas.ch"
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    return record


if __name__=="__main__":

    rates_by_clade =  pd.read_csv("results/synonymous_mut_rates/rates_by_clade.csv")

    clades = rates_by_clade.clade.unique()

    equilibrium_probabilities_by_clade = {}
    for clade in clades:
        subset = rates_by_clade.loc[rates_by_clade.clade==clade]
        tmp_spectrum = {}
        for rate in subset.itertuples():
            tmp_spectrum[rate.mut_type] = rate.rate

        equilibrium_probabilities_by_clade[clade] = equilibrium_probabilities(spectrum_to_matrix(tmp_spectrum))


    accessions = ["KY352407.1", "MN908947","FJ425184", "FJ938066", "HM211098",
                  "NC_017083", "HM211101", "LC494161", "KJ713296", "MF083115", "MH002339",
                  "MH687973", "OL674078", "MG197719", "OK017819", "MW064263"]
    empirical_frequencies_references = {}
    for accession in accessions:
        ref = get_sequence_from_genbank(accession)
        ref_genome = get_ffsyn(ref)
        empirical_frequencies_references[accession] = empirical_frequencies({i:r for i,r in ref_genome.loc[ref_genome.four_fold & (~ref_genome.four_fold.isna())].nuc.value_counts().iteritems()})

    plt.figure()
    for c, p in equilibrium_probabilities_by_clade.items():
        plt.plot(p,label=c)
    for a, f in empirical_frequencies_references.items():
        plt.plot(f, label=a, ls='--')


    plt.xticks(np.arange(len(alphabet)), alphabet)

    plt.legend(ncol=3)
    plt.yscale('log')
    plt.ylabel('equilibrium probabilities')

