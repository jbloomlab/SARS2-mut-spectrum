import urllib, json, os
from Bio import SeqIO
from collections import defaultdict
from equilibrium_frequencies import equilibrium_probabilities, empirical_frequencies, spectrum_to_matrix

def assign_mutation_counts(n, mutation_counts):
    for mut in n['branch_attrs']['mutations'].get('nuc',[]):
        a, pos, d = mut[0], int(mut[1:-1]), mut[-1]
        if pos in ff_syn_pos and a in alphabet and d in alphabet:
            mutation_counts[f"{a}to{d}"] += 1
    if "children" in n:
        for c in n['children']:
            assign_mutation_counts(c, mutation_counts)

def get_nextstrain_tree(lineage, segment, resolution):
    url_str = f"https://nextstrain.org/charon/getDataset?prefix=/groups/neherlab/flu/seasonal/{lineage}/{segment}/{resolution}"
    with urllib.request.urlopen(url_str) as url:
        data = json.loads(url.read())
    return data['tree']

def get_reference(lineage, segment):
    from io import StringIO
    url_str = f"https://raw.githubusercontent.com/nextstrain/seasonal-flu/5beeb5e53e073ef4ef66df17bdb53e76958c6a0c/config/reference_{lineage}_{segment}.gb"
    with urllib.request.urlopen(url_str) as url:
        ref = SeqIO.read(StringIO(url.read().decode()), "genbank")

    return ref

if __name__ == "__main__":

    alphabet = 'ACGT'
    lineages = {"h3n2":"60y", "h1n1pdm":"12y", "vic":"60y", "yam":"60y"}

    results = {}

    for lineage in lineages:
        ff_syn_pos = set()
        ff_syn_nucs = defaultdict(int)
        mutation_counts = defaultdict(int)

        for segment in ['pb1', 'pb2', 'pa', 'ha', 'np', 'na']:
            T = get_nextstrain_tree(lineage, segment, lineages[lineage])
            ref = get_reference(lineage, segment)

            for feature in ref.features:
                if feature.type == "CDS":
                    for fpos, ref_pos in enumerate(range(feature.location.start, feature.location.end)):
                        if fpos%3==2:
                            if ref.seq[ref_pos-2:ref_pos] in ['TC','CT', 'CC', 'CG', 'AC', 'GT', 'GC','GG']:
                                ff_syn_pos.add(ref_pos+1)
                                ff_syn_nucs[ref.seq[ref_pos]] += 1

            assign_mutation_counts(T, mutation_counts)

        res = {}
        res['mutation_counts'] = mutation_counts
        res['mutation_spectrum'] = {k: v/ff_syn_nucs[k[0]] for k, v in mutation_counts.items()}
        res['ff_syn_states'] = ff_syn_nucs
        res['empirical_frequencies'] = {n:float(x) for n,x in zip(alphabet, empirical_frequencies(ff_syn_nucs))}
        res["equilibrium_frequencies"] = {n:float(x) for n,x in zip(alphabet,
                            equilibrium_probabilities(spectrum_to_matrix(res['mutation_spectrum'])))}
        results[lineage] = res

    os.makedirs("results/flu_spectra", exist_ok=True)
    with open("results/flu_spectra/flu_spectra.json", "w") as f:
        json.dump(results, f)


