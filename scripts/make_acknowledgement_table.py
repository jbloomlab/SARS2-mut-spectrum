from compare_flu_spectra import get_nextstrain_tree

def get_labs(node, isolates):
    if "children" in node:
        for c in node["children"]:
            get_labs(c, isolates)
    else:
        if "accession" in node["node_attrs"] and "originating_lab" in node["node_attrs"]:
            isolates.append({"epi_isl": node["node_attrs"]["accession"],
                         "submitting_lab": node["node_attrs"].get("submitting_lab", {'value':'NA'})['value'],
                         "originating_lab": node["node_attrs"]["originating_lab"]['value']})



if __name__=="__main__":

    lineages = {"h3n2":"60y", "h1n1pdm":"12y", "vic":"60y", "yam":"60y"}

    isolates = []

    for lineage in lineages:
        T = get_nextstrain_tree(lineage, 'ha', lineages[lineage])
        get_labs(T,isolates)

    import pandas as pd

    pd.DataFrame(isolates).to_csv('flu_acknowledgement.tsv', sep='\t')
