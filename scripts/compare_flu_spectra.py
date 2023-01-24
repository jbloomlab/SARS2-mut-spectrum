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

def get_nextstrain_tree(dset):
    url_str = f"https://nextstrain.org/charon/getDataset?prefix={dset}"
    with urllib.request.urlopen(url_str) as url:
        data = json.loads(url.read())
    return data['tree']

def get_reference(url_str):
    from io import StringIO
    with urllib.request.urlopen(url_str) as url:
        ref = SeqIO.read(StringIO(url.read().decode()), "genbank")

    return ref

if __name__ == "__main__":

    alphabet = 'ACGT'
    flu_segments = ['pb1', 'pb2', 'pa', 'ha', 'np', 'na']
    builds = {
        "flu_h3n2":{"url":"/groups/neherlab/flu/seasonal/h3n2/{segment}/60y", "segments":flu_segments, "ref":"https://raw.githubusercontent.com/nextstrain/seasonal-flu/5beeb5e53e073ef4ef66df17bdb53e76958c6a0c/config/reference_h3n2_{segment}.gb"},
        "flu_h1n1pdm":{"url":"/groups/neherlab/flu/seasonal/h1n1pdm/{segment}/12y", "segments":flu_segments, "ref":"https://raw.githubusercontent.com/nextstrain/seasonal-flu/5beeb5e53e073ef4ef66df17bdb53e76958c6a0c/config/reference_h1n1pdm_{segment}.gb"},
        "flu_vic":{"url":"/groups/neherlab/flu/seasonal/vic/{segment}/60y", "segments":flu_segments, "ref":"https://raw.githubusercontent.com/nextstrain/seasonal-flu/5beeb5e53e073ef4ef66df17bdb53e76958c6a0c/config/reference_vic_{segment}.gb"},
        "flu_yam":{"url":"/groups/neherlab/flu/seasonal/yam/{segment}/60y", "segments":flu_segments, "ref":"https://raw.githubusercontent.com/nextstrain/seasonal-flu/5beeb5e53e073ef4ef66df17bdb53e76958c6a0c/config/reference_yam_{segment}.gb"},
        "rsv-a":{"url":"/rsv/a/{segment}", "segments":["genome"], 'ref':"https://raw.githubusercontent.com/nextstrain/rsv/master/config/areference.gbk"},
        "rsv-b":{"url":"/rsv/b/{segment}", "segments":["genome"], 'ref':"https://raw.githubusercontent.com/nextstrain/rsv/master/config/breference.gbk"},
        "evD68":{"url":"enterovirus/d68/{segment}", "segments":["genome"], "ref":"https://raw.githubusercontent.com/nextstrain/enterovirus_d68/master/genome/config/ev_d68_reference_genome.gb"},
        "denv1":{"url":"dengue/denv1{segment}", "segments":[""], "ref":"https://raw.githubusercontent.com/nextstrain/dengue/main/config/reference_dengue_denv1.gb"},
        "denv2":{"url":"dengue/denv2{segment}", "segments":[""], "ref":"https://raw.githubusercontent.com/nextstrain/dengue/main/config/reference_dengue_denv2.gb"},
        "denv3":{"url":"dengue/denv3{segment}", "segments":[""], "ref":"https://raw.githubusercontent.com/nextstrain/dengue/main/config/reference_dengue_denv3.gb"},
        "denv4":{"url":"dengue/denv4{segment}", "segments":[""], "ref":"https://raw.githubusercontent.com/nextstrain/dengue/main/config/reference_dengue_denv4.gb"},
        "WNV":{"url":"WNV/NA{segment}", "segments":[""], "ref":"https://raw.githubusercontent.com/grubaughlab/WNV-nextstrain/master/config/reference.gb"},
}

    results = {}

    for bname, build in builds.items():
        ff_syn_pos = set()
        ff_syn_nucs = defaultdict(int)
        mutation_counts = defaultdict(int)

        for segment in build["segments"]:
            T = get_nextstrain_tree(build["url"].format(segment=segment))
            ref = get_reference(build["ref"].format(segment=segment))

            for feature in ref.features:
                if feature.type == "CDS":
                    if 'rsv' in bname and feature.qualifiers['gene'][0]=='G':
                        continue
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
        results[bname] = res

    os.makedirs("results/flu_spectra", exist_ok=True)
    with open("results/flu_spectra/flu_spectra.json", "w") as f:
        json.dump(results, f)


