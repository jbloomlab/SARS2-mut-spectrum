import json



if __name__=="__main__":
    with open("results/flu_spectra/flu_spectra.json") as f:
        results = json.load(f)

    import matplotlib.pyplot as plt

    ms = 50
    fig = plt.figure(figsize=(6, 3.25))
    colors = ["tab:blue", "tab:orange", "tab:red", "tab:cyan"]
    xmax = 0.8  # harmonize with other plot
    for li, lineage in enumerate(results):
        res = results[lineage]
        for ki, (k, v) in enumerate(res['empirical_frequencies'].items()):
            plt.scatter(v, li-0.2, color=colors[ki], marker='o', label=k if li==0 else None, s=ms)
        plt.text(0.01, li-0.2, 'empirical', ha='left', va='center', fontsize=10)

        for ki, (k, v) in enumerate(res['equilibrium_frequencies'].items()):
            plt.scatter(v, li+0.2, color=colors[ki], marker='+', label=k if li==0 else None, s=ms)
        plt.text(0.01, li+0.2, 'predicted', ha='left', va='center', fontsize=10)

        for ki, (k, v) in enumerate(res['empirical_frequencies'].items()):
            plt.plot([v, res['equilibrium_frequencies'][k]], [li-0.2, li+0.2], color=colors[ki])

        plt.fill_between([0,xmax], li-0.5, li+0.5, color='grey', alpha=0.1 + 0.2*(li%2), edgecolor='none')

    label_names = {
        "yam": "B/Yamagata",
        "vic": "B/Victoria",
        "h1n1pdm": "A/H1N1pdm",
        "h3n2": "A/H3N2",
    }
    plt.yticks(range(len(results)), [label_names[name] for name in results])
    plt.xlim([0, xmax])
    plt.ylim([-0.5, 3.5])
    plt.tick_params(axis='y', labelsize=12)
    plt.xlabel("nucleotide frequency", fontsize=12, fontweight="bold")
    plt.savefig("results/flu_spectra/flu_spectra.pdf", bbox_inches='tight')
