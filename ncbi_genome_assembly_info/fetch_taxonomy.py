import argparse
import csv
import logging
import sys
import time
from typing import List, Optional, Tuple

from Bio import Entrez


def setup_logging(level: str) -> None:
    """
    Configure the logging settings.

    Args:
        level (str): Logging level as a string (e.g., 'DEBUG', 'INFO').
    """
    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        print(f"Invalid log level: {level}")
        sys.exit(1)
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extract NCBI TaxID, common names, and assembly information for a list of species."
    )
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Path to the input text file containing species names (one per line).'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='species_taxonomy.csv',
        help='Path to the output CSV file. Default is species_taxonomy.csv'
    )
    parser.add_argument(
        '-e', '--email',
        type=str,
        required=True,
        help='Your email address (required by NCBI Entrez).'
    )
    parser.add_argument(
        '-k', '--api_key',
        type=str,
        default='',
        help='Your NCBI API key (optional, but recommended for higher rate limits).'
    )
    parser.add_argument(
        '-l', '--log_level',
        type=str,
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help='Set the logging level. Default is INFO.'
    )
    parser.add_argument(
        '-d', '--delay',
        type=float,
        default=0.4,
        help='Delay between API requests in seconds. Default is 0.4.'
    )
    return parser.parse_args()


def read_species_names(file_path: str) -> List[str]:
    """
    Read species names from a text file, one per line.

    Args:
        file_path (str): Path to the input file.

    Returns:
        List[str]: List of species names.
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            species = [line.strip() for line in file if line.strip()]
        logging.debug(f"Read {len(species)} species names from {file_path}")
        return species
    except FileNotFoundError:
        logging.error(f"Input file not found: {file_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading input file: {e}")
        sys.exit(1)


def get_taxid(species_name: str) -> Optional[str]:
    """
    Retrieve the TaxID for a given species name using Entrez.esearch.

    Args:
        species_name (str): The scientific name of the species.

    Returns:
        Optional[str]: The TaxID if found, else None.
    """
    try:
        handle = Entrez.esearch(db="taxonomy", term=species_name, retmode="xml")
        results = Entrez.read(handle)
        handle.close()
        if results["IdList"]:
            taxid = results["IdList"][0]
            logging.debug(f"Found TaxID {taxid} for species '{species_name}'")
            return taxid
        else:
            logging.warning(f"No TaxID found for species '{species_name}'")
            return None
    except Exception as e:
        logging.error(f"Error fetching TaxID for '{species_name}': {e}")
        return None


def get_common_name(taxid: str) -> Optional[str]:
    """
    Retrieve the common name for a given TaxID using Entrez.efetch.

    Args:
        taxid (str): The NCBI Taxonomy ID.

    Returns:
        Optional[str]: The common name if available, else the scientific name.
    """
    if not taxid.isdigit():
        logging.error(f"Invalid TaxID '{taxid}'. It should contain only digits.")
        return None

    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if records and 'ScientificName' in records[0]:
            record = records[0]
            if 'OtherNames' in record and 'CommonName' in record['OtherNames']:
                common_name = record['OtherNames']['CommonName']
                # Ensure common_name is a string
                if isinstance(common_name, list):
                    common_name = '; '.join(common_name)
                logging.debug(f"Found common name '{common_name}' for TaxID {taxid}")
                return common_name
            else:
                scientific_name = record['ScientificName']
                logging.info(f"No common name for TaxID {taxid}; using scientific name '{scientific_name}'")
                return scientific_name
        else:
            logging.warning(f"No records found for TaxID {taxid}")
            return None
    except Exception as e:
        logging.error(f"Error fetching common name for TaxID '{taxid}': {e}")
        return None


def get_reference_genome(taxid: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Retrieve the GenBank accession and assembly name for the reference genome of a given TaxID.

    Args:
        taxid (str): The NCBI Taxonomy ID.

    Returns:
        Tuple[Optional[str], Optional[str]]: (GenBank Accession, Assembly Name)
    """
    if not taxid.isdigit():
        logging.error(f"Invalid TaxID '{taxid}'. It should contain only digits.")
        return (None, None)

    try:
        # Search for assemblies for the given TaxID, prioritizing reference genomes
        search_term = f'txid{taxid}[Organism:exp]'
        handle = Entrez.esearch(db="assembly", term=search_term, retmode="xml", retmax=100)
        search_results = Entrez.read(handle)
        handle.close()

        if not search_results["IdList"]:
            logging.info(f"No assembly found for TaxID {taxid}")
            return (None, None)

        assembly_ids = search_results["IdList"]
        logging.debug(f"Found Assembly IDs {assembly_ids} for TaxID {taxid}")

        # Fetch assembly summaries for all assemblies
        handle = Entrez.esummary(db="assembly", id=",".join(assembly_ids), retmode="xml")
        summaries = Entrez.read(handle)
        handle.close()

        if not summaries or 'DocumentSummarySet' not in summaries:
            logging.warning(f"No assembly records found for TaxID {taxid}")
            return (None, None)

        assembly_list = summaries['DocumentSummarySet']['DocumentSummary']

        # Look for reference genomes
        reference_assemblies = [
            asm for asm in assembly_list if asm.get('RefSeq_category') == 'reference genome'
        ]

        # If no reference genome is found, consider representative genome
        if not reference_assemblies:
            reference_assemblies = [
                asm for asm in assembly_list if asm.get('RefSeq_category') == 'representative genome'
            ]

        # Select the first reference genome found
        if reference_assemblies:
            assembly = reference_assemblies[0]
            logging.debug(f"Selected Reference Assembly ID {assembly['AssemblyAccession']} for TaxID {taxid}")
        else:
            # If no reference or representative genome, select the latest assembly
            assembly = assembly_list[0]
            logging.debug(f"No reference genome found. Selected Assembly ID {assembly['AssemblyAccession']} for TaxID {taxid}")

        # Extract GenBank Accession and Assembly Name
        genbank_accession = assembly.get('AssemblyAccession', 'Not Available')
        assembly_name = assembly.get('AssemblyName', 'Not Available')

        logging.debug(f"GenBank Accession: {genbank_accession}, Assembly Name: {assembly_name} for TaxID {taxid}")

        return (genbank_accession, assembly_name)

    except Exception as e:
        logging.error(f"Error fetching assembly for TaxID '{taxid}': {e}")
        return (None, None)


def write_to_csv(output_file: str, data: List[Tuple[str, str, str, str, str]]) -> None:
    """
    Write the taxonomy and assembly data to a CSV file.

    Args:
        output_file (str): Path to the output CSV file.
        data (List[Tuple[str, str, str, str, str]]): List of tuples containing
            species name, TaxID, common name, GenBank Accession, and Assembly Name.
    """
    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['Species Name', 'TaxID', 'Common Name', 'GenBank Accession', 'Assembly Name']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for species, taxid, common_name, genbank_acc, assembly_name in data:
                writer.writerow({
                    'Species Name': species,
                    'TaxID': taxid if taxid else 'Not Found',
                    'Common Name': common_name if common_name else 'Not Available',
                    'GenBank Accession': genbank_acc if genbank_acc else 'Not Available',
                    'Assembly Name': assembly_name if assembly_name else 'Not Available'
                })
        logging.info(f"Data successfully written to {output_file}")
    except Exception as e:
        logging.error(f"Error writing to CSV file '{output_file}': {e}")
        sys.exit(1)


def main() -> None:
    """
    Main function to execute the taxonomy and assembly extraction process.
    """
    args = parse_arguments()
    setup_logging(args.log_level)

    # Ensure the email is a single, properly formatted string
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
        logging.debug("NCBI API key set.")
    else:
        logging.debug("No NCBI API key provided.")

    species_list = read_species_names(args.input)
    logging.info(f"Total species to process: {len(species_list)}")

    results = []

    for idx, species in enumerate(species_list, start=1):
        logging.info(f"Processing ({idx}/{len(species_list)}): {species}")
        taxid = get_taxid(species)
        if taxid:
            common_name = get_common_name(taxid)
            genbank_acc, assembly_name = get_reference_genome(taxid)
        else:
            common_name = None
            genbank_acc = None
            assembly_name = None

        results.append((
            species,
            taxid if taxid else 'Not Found',
            common_name if common_name else 'Not Available',
            genbank_acc if genbank_acc else 'Not Available',
            assembly_name if assembly_name else 'Not Available'
        ))

        # Respect NCBI's rate limiting
        logging.debug(f"Sleeping for {args.delay} seconds to respect rate limits.")
        time.sleep(args.delay)

    write_to_csv(args.output, results)
    logging.info("Data extraction complete.")


if __name__ == "__main__":
    main()
