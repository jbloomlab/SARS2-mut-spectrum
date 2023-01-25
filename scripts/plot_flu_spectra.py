import json



if __name__=="__main__":
    with open("results/other_virus_spectra.json") as f:
        results = json.load(f)

    import matplotlib.pyplot as plt

    ms = 50
    fig = plt.figure(figsize=(15,8))
    for li, lineage in enumerate(results):
        res = results[lineage]
        for ki, (k, v) in enumerate(res['empirical_frequencies'].items()):
            plt.scatter(v, li-0.2, color=f'C{ki}', marker='o', label=k if li==0 else None, s=ms)
        plt.text(0.03, li-0.2, 'empirical', ha='left', va='center', fontsize=12)

        for ki, (k, v) in enumerate(res['equilibrium_frequencies'].items()):
            plt.scatter(v, li+0.2, color=f'C{ki}', marker='+', label=k if li==0 else None, s=ms)
        plt.text(0.03, li+0.2, 'predicted', ha='left', va='center', fontsize=12)

        for ki, (k, v) in enumerate(res['empirical_frequencies'].items()):
            plt.plot([v, res['equilibrium_frequencies'][k]], [li-0.2, li+0.2], color=f'C{ki}')

        plt.fill_between([0,0.5], li-0.5, li+0.5, color='grey', alpha=0.1 + 0.2*(li%2), edgecolor='none')

    plt.yticks(range(len(results)), list(results.keys()))
    plt.tick_params(axis='y', labelsize=12)
    plt.xlabel("Frequency")
    plt.savefig("results/other_virus_spectra.pdf", bbox_inches='tight')
