from typing import List, Dict, Optional
import pandas as pd
import os
import sys
import re
import numpy as np
from glob import glob
import typer
import datetime
from tqdm import tqdm

app = typer.Typer(add_completion=False)

class Filter:
    """class that will apply the regex filter to find either the missense/nonse variants or all of them"""
    def __init__(self, find_all_snps: bool) -> None:
        # This attribute will be tree or false
        self.snps_to_gather: bool = find_all_snps

    def filter_for_pathogenicity(self, variant_df: pd.DataFrame) -> pd.DataFrame:
        """function that will filter the provided dataframe for variants based on if an argument is passed. If the self attribute is True then all the variants will be returned otherwise just the pathogenic ones will be 
        Parameter
        _________
        variant_df : pd.DataFrame
            dataframe that has the information about the gene of interest such as what mega probes there are and what the mutation does
            
        Return
        ______
        pd.DataFrame
            dataframe that either only has missense/nonsense variants or has all the variants
        """
        # if the user chooses to keep all snps then
        if self.snps_to_gather:
            
            return variant_df[variant_df["Mutation(s)"].str.contains(r'Missense|Nonsense|Silent|Synonymous', na=False)]

        else:

            return variant_df[variant_df["Mutation(s)"].str.contains(r'Missense|Nonsense', na=False)]

def get_files(file_directory: str) -> List[str]:
    """Function that can get all the annotated excel files from the specified directory
    Parameters
    __________
    file_directory : str
        string that list the filepath to all the xlsx files that have the MEGA array probes
    
    Returns
    _______
    List[str]
        returns a list of strings that has all the file paths to the annotation files in it
    """
    cur_dir: str = os.getcwd()

    annotation_file_list: List[str] = []

    os.chdir(file_directory)

    # gathering all the files that have and .xlsx extension and iterating through them
    for file in glob("*.xlsx"):

        full_file_path: str = os.path.join(file_directory, file)

        # check to make sure that the file has the format ChrXX where XX is a number
        match_string: str = re.search(r'Chr\d\d', file)

        if match_string:

            annotation_file_list.append(full_file_path)
    
    # if there are no files detected in the provided directory then the program needs to end
    if len(annotation_file_list) == 0:

        print("There were no annotation files found in the specified directory. Please ensure that there are excel annotation files. Ending program now...")
        sys.exit(0)

    os.chdir(cur_dir)

    return annotation_file_list

def get_gene_list(gene_file: str) -> List[str]:
    """function to get a list of genes of interest
    Parameter
    _________
    gene_file : str
        string that contains the file path for the file contains genes of interest 
    
    Returns
    _______
    List[str]
        list of strings that has the genes from the file of interest
    """
    gene_df: pd.DataFrame = pd.read_csv(gene_file, sep="\t")
    
    if "gene" not in gene_df.columns:
        print("expected the input file with gene targets to have a column named gene. This column was not found")
        sys.exit(0)

    return gene_df.gene.values.tolist()

def filter(row: str, list: List[str]):

    if type(row) == datetime.datetime:
        return np.nan

    split_row: list = row.split(",")

    for element in split_row:

        if element in list:
            return ",".join(split_row)
            
    
    return np.nan

def remove_previous_file(output_name: str):
    """function to check if the output file exist from a previous run and deletes it if exist
    Parameter
    _________
    output_name : str
        string that list the output file path
    """
    try:
        os.remove(output_name)

    except FileNotFoundError:
        pass

def find_variant_snps(file_list: List[str], gene_list: List[str], output_path: str, aggregate_all_probes: bool):
    """Function to find the variant snps on the mega probe for a specific gene 
    Parameters
    __________
    file_list : List[str]
        list of files for each chromosome mega annotation file
    
    gene_list : List[str]
        list of all the genes of interest
    
    output_path : str
        string that list the output file path (including the filename) 

    aggregate_all_probes : bool
        boolean value for whether or not the user wants to keep all the probes or only the missense/nonsense ones
    """
    # looking for the gene in each of the files

    # creating a filter type
    filter_object: Filter = Filter(aggregate_all_probes)

    result = 0
    for i in tqdm(range(len(file_list))):
        
        file = file_list[i]

        if re.match(r".*Chr06.*", file): 

            file_df: pd.DataFrame = pd.read_excel(file, sheet_name="cleaned")
        
        else:
            file_df: pd.DataFrame = pd.read_excel(file)

        if "Gene(s)" not in file_df.columns:

            print("expected a column called Gene(s) to be in the file")
            sys.exit(1)

        file_df = file_df.dropna()
        
        file_df["Genes"] = file_df["Gene(s)"].apply(lambda row: filter(row, gene_list))
        
        filtered_df: pd.DataFrame = file_df.dropna()

        if not filtered_df.empty:

            filtered_df = filter_object.filter_for_pathogenicity(filtered_df)
            
            if re.match(r".*Chr06.*", file): 
                
                file_info_dict: Dict[str, Optional[List]] = {
                    "name":filtered_df.IlmnID.values.tolist(), 
                    "RsID":filtered_df["RS Name"].values.tolist(),
                    "Chr": filtered_df.Chr.values.tolist(),
                    "MapInfo": filtered_df.MapInfo.values.tolist(),
                    "Alleles": filtered_df.SNP.values.tolist(),
                    "Transcript":None,
                    "Gene(s)":filtered_df["Gene(s)"].values.tolist(),
                    "In-exon": filtered_df["In-exon"].values.tolist(),
                    "Mutation(s)": filtered_df["Mutation(s)"].values.tolist()
                    }
                
                file_info_df: pd.DataFrame = pd.DataFrame.from_dict(file_info_dict)

                file_info_df = file_info_df.drop(["Genes"], axis=1)

                filtered_df = file_info_df

            else:

                filtered_df = filtered_df.drop(["Genes"], axis=1)

                
        if result == 0:

            filtered_df.to_csv(output_path, sep="\t", mode="a+", index=False)

        else:

            filtered_df.to_csv(output_path, sep="\t", mode="a+", index=False, header=False)

        result += 1

@app.command()
def get_variants(
    annotation_file_dir: str = typer.Argument(
        ..., help="String that list the filepath to the directory that contains the annotationm files for each chromosome that list the variant probes that are on vanderbilts mega array"
    ), 
    output_filepath: str = typer.Argument(
        ..., help="String that list the output filepath including the file name to output the resulting file into "
    ), 
    gene_target_file: str = typer.Argument(
        ..., help="String that list the file path to a text file that has one column that has gene targets for the program. This file should have a column named 'gene'"
    ), 
    gather_all_snps: bool = typer.Option(False, help="Argument to indicate if the program should return all variant probes whose mutation consequences are known or if the program should only return the missense/nonsense probes")
    ):
    """
    main function to generate a file that contains all pathogenic snps for genes of interests from the annotated mega files
    """
    print("Aggregating variants of the mega array for the genes of interest...\n")
    print(f"Annotation File Directory: {annotation_file_dir}")
    print(f"Gene Targets File: {gene_target_file}")
    print(f"Output File: {output_filepath}\n")

    remove_previous_file(output_filepath)

    # getting a list that has all of the annotation files

    file_list: List[str] = get_files(annotation_file_dir)

    gene_list: List[str] = get_gene_list(gene_target_file)

    # getting all the snps for a specific variant
    find_variant_snps(file_list, gene_list, output_filepath, gather_all_snps)

if __name__ == '__main__':
    app()