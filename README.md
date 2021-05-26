# Mega_variants_aggregator

* author: J.T. Baker

## Purpose:
___
cli tool to quickly gather all missense/nonsense variants associated with genes of interest

## Project Structure:
____
* **get_variant_snps.py**: main file for the program that contains the cli interface and the necessary functions to aggregate all of the mega probes

## Expected Input:
____
* **Annotation file directory**: Expects the user to provide the path to a directory that contains the annotation files for the probes on the mega array. These files should be excel files and should have columns: name, RsID, Chr, MapInfo, Alleles, Transcript, Gene(s), In-exon, Mutation(s) (capitalization and spelling need to be the same as what is mentioned)
* **Genes of Interest**: Tab separated text file that has a column named genes which list the genes that the user wishes to get probes for

## Notes on using the program:
___
* To find all necessary arguments that can be pass to the function run
```
python3 get_variant_snps.py --help
```

* An example of how to run the program is shown below. The user can replace each argument with the correct filepath

```
python3 get_variant_snps.py  annotation_file_path output_file_path genes_of_interest_file_path
```
