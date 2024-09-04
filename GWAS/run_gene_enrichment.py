import pandas as pd # DataFrame
import re # Regular expressions
import os # Files and file paths
import argparse # Importing arguments from command-line call
from scipy.stats import fisher_exact # Enrichment function (Fisher's Exact Test)
from statsmodels.stats.multitest import fdrcorrection # For FDR correction
from tqdm import tqdm # Progress bar

def extract_fbgn_from_annot(fgenes):
    # List storing our results
    fbgn_values = []
    
    # Define the regular expression pattern to match FBgn values and the gene values
    pattern = r'(FBgn\d+)\|([^|]+)'
    
    # For all rows of DataFrame
    for i in range(0, len(fgenes)):
        annot = fgenes[i]
        matches = re.findall(pattern, annot)
        fbgn_values.append([match[0] for match in matches])
    
    fbgn_unique_values = set([item for sublist in fbgn_values for item in sublist])
    
    return fbgn_unique_values

def enrichment(phenotype_of_interest, fphenotyping, fbackground, fgwas):
    # Filter the DataFrame
    fbgn_phenotype = set(fphenotyping[fphenotyping['feature_id'].isin(phenotype_of_interest)]['GeneID'])
    #print("There are", len(fbgn_phenotype), " genes associated with this phenotype:", phenotype_of_interest)
    fbgn_phenotype = fbgn_phenotype.intersection(fbackground)
    #print("There are", len(fbgn_phenotype), " genes (overlapping with background) associated with this phenotype:", phenotype_of_interest)
    
    overlap = fgwas.intersection(set(fbgn_phenotype))
    #print("There are", len(overlap), " genes overlapping")
    
    # Create the 2x2 contingency table
    contingency_table = pd.DataFrame([[len(overlap), len(fbgn_phenotype) - len(overlap)], [len(fgwas) - len(overlap), len(fbackground) - len(fgwas) - len(fbgn_phenotype) + len(overlap)]])
    #print(contingency_table)
    
    # Perform Fisher's exact test
    results_fisher = fisher_exact(contingency_table)
    return results_fisher.statistic, results_fisher.pvalue, len(overlap), len(fbgn_phenotype), len(fgwas)

def run_enrichment(fdatagwas, fgeneset, fbackground):
    fbgn_gwas = extract_fbgn_from_annot(fdatagwas.gene_annotation)
    #print("There are", len(fbgn_gwas), "significant genes for this GWAS results")
    fbgn_gwas = fbgn_gwas.intersection(fbackground)
    #print("There are", len(fbgn_gwas), "significant genes (overlapping with background) for this GWAS results")

    # Create a mapping table
    mapping_table = fgeneset[['feature_id', 'feature_name']].drop_duplicates(subset=['feature_id'])
    mapping_table.set_index('feature_id', inplace=True)
    
    # List to store the results
    results = []
    for feature in tqdm(set(fgeneset.feature_id), desc="Processing features"):
        oddsratio, pvalue, overlap_len, phenotype_len, gwas_len  = enrichment([feature], fgeneset, fbackground, fbgn_gwas)
        results.append([feature, mapping_table.loc[feature, 'feature_name'], oddsratio, pvalue, overlap_len, phenotype_len, gwas_len])
    return pd.DataFrame(results, columns=['feature_id', 'feature_name', 'odds_ratio', 'p_value', 'nb_gene_overlap', 'nb_gene_phenotype', 'nb_gene_gwas']).sort_values(by='p_value')

def run(fstudy, fpheno, fsex, fgwas, fgenesetdb, fbackgrounddb, ffilter = 0.01):   
    # Loading data
    file_path = f"{fgwas}/S{fstudy}_{fpheno}_{fsex}.glm.linear.fdr.annot.tsv.gz"
    if not os.path.exists(file_path):
        file_path = f"{fgwas}/S{fstudy}_{fpheno}_{fsex}.glm.logistic.hybrid.fdr.annot.tsv.gz"
    data_gwas_results = pd.read_csv(file_path, sep = "\t")
    
    # Filtering
    data_gwas_results = data_gwas_results[data_gwas_results.P < ffilter]

    # Compute enrichment for each phenotype
    results = run_enrichment(data_gwas_results, fgenesetdb, fbackgrounddb)

    # Prepare result DataFrame
    results["nb_gene_background"] = len(fbackgrounddb)
    results["id_study"] = fstudy
    results["id_phenotype"] = fpheno
    results["sex"] = fsex
    results["path_gwas"] = fgwas
    results["nb_variant_significant"] = len(data_gwas_results)
    results["filter_threshold"] = ffilter
    return results

def main():
    # Parsing the command-line arguments
    parser = argparse.ArgumentParser(description="Geneset enrichment of GWAS results")
    parser.add_argument("directory_gwas", type=str, help="The path where all filtered results are stored")
    parser.add_argument("geneset_db", type=str, help="The Flybase database of phenotypes")
    parser.add_argument("ncsu_variant_annot", type=str, help="The dm3 variant annotation from NCSU")
    parser.add_argument("study_id", type=int, help="The study ID in the DB (for the desired phenotype)")
    parser.add_argument("phenotype_id", type=int, help="The phenotype ID in the DB")
    parser.add_argument("sex", type=str, help="The sex of the desired phenotype")
    parser.add_argument("filter_threshold", type=float, help="The filtering threshold")
    parser.add_argument("output_file", type=str, help="The path where all gene enrichment results are output")
    args = parser.parse_args()
    
    # Loading already processed dataset with all phenotypes/alleles/genes
    data_geneset = pd.read_csv(args.geneset_db, sep="\t")
    #print(f"Size of data_geneset dataset: {data_geneset.shape}") # Should be (587739 , 9)
    
    # Looking for background (all genes we could intersect with flybase)
    with open(args.ncsu_variant_annot, "r") as file:
        fbgn_background = [line.strip() for line in file]
    #print("Background is", len(fbgn_background), "genes") # Should be 12053
    
    # Run the main method
    results = run(args.study_id, args.phenotype_id, args.sex, args.directory_gwas, data_geneset, fbgn_background, args.filter_threshold)

    # Correct p-values by performing FDR correction using the Benjamini-Hochberg procedure
    rejected, pvals_corrected = fdrcorrection(results.p_value, alpha=0.05, method='indep')
    results['FDR'] = pvals_corrected
    
    # Write results to a TSV file
    results.to_csv(args.output_file, index=False, sep="\t")
    
if __name__ == "__main__":
    main()
