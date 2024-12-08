import argparse
import csv
import logging
import sys
import time
from typing import List, Optional, Tuple, Dict

from Bio import Entrez

VALID_TAXONOMY_RANKS = {
    'superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum',
    'subphylum', 'superclass', 'class', 'subclass', 'infraclass',
    'superorder', 'order', 'suborder', 'infraorder', 'parvorder',
    'superfamily', 'family', 'subfamily', 'tribe', 'subtribe',
    'genus', 'subgenus', 'species', 'subspecies', 'varietas', 'forma'
}


def setup_logging(level: str) -> None:
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
    parser = argparse.ArgumentParser(
        description="Extract NCBI TaxID, common names, assembly information, and additional taxonomy levels for a list of species."
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
    parser.add_argument(
        '-t', '--tax_levels',
        type=str,
        default='family,subfamily',
        help='Comma-separated taxonomy levels to include (e.g., family,subfamily). Default is family,subfamily.'
    )
    parser.add_argument(
        '-a', '--assembly_preference',
        type=str,
        default='auto',
        choices=['auto', 'refseq', 'genbank'],
        help='How to select assemblies: auto (default), refseq (GCF_), or genbank (GCA_).'
    )
    args = parser.parse_args()

    # Process taxonomy levels
    tax_levels = [level.strip().lower() for level in args.tax_levels.split(',') if level.strip()]
    invalid_levels = [level for level in tax_levels if level not in VALID_TAXONOMY_RANKS]
    if invalid_levels:
        print(f"Invalid taxonomy levels: {', '.join(invalid_levels)}")
        print(f"Valid taxonomy levels are: {', '.join(sorted(VALID_TAXONOMY_RANKS))}")
        sys.exit(1)
    args.tax_levels = tax_levels

    return args


def read_species_names(file_path: str) -> List[str]:
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


def get_taxonomy_levels(taxid: str, desired_levels: List[str]) -> Dict[str, str]:
    taxonomy_info = {level.capitalize(): 'Not Available' for level in desired_levels}
    if not taxid.isdigit():
        logging.error(f"Invalid TaxID '{taxid}'. It should contain only digits.")
        return taxonomy_info

    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if records:
            lineage = records[0].get('LineageEx', [])
            for rank_info in lineage:
                rank = rank_info.get('Rank', '').lower()
                name = rank_info.get('ScientificName', '')
                if rank in desired_levels:
                    taxonomy_info[rank.capitalize()] = name
            # Also check the current rank (species)
            current_rank = records[0].get('Rank', '').lower()
            current_name = records[0].get('ScientificName', '')
            if current_rank in desired_levels:
                taxonomy_info[current_rank.capitalize()] = current_name
            logging.debug(f"Retrieved taxonomy levels for TaxID {taxid}: {taxonomy_info}")
        else:
            logging.warning(f"No taxonomy records found for TaxID {taxid}")
    except Exception as e:
        logging.error(f"Error fetching taxonomy levels for TaxID '{taxid}': {e}")
    
    return taxonomy_info


def main_numeric_part(accession: str) -> str:
    # GCF_902806735.2 -> 902806735
    # GCA_902806735.1 -> 902806735
    # Extract the main numeric part before the dot
    if '_' in accession:
        numeric_part = accession.split('_', 1)[1]
        return numeric_part.split('.', 1)[0]
    return accession


def select_assembly(candidates: List[Dict], preference: str, has_refseq_reference: bool) -> Optional[Dict]:
    # Apply the hierarchy:
    # 1. reference genome
    reference = [asm for asm in candidates if asm.get('RefSeq_category') == 'reference genome']
    if reference:
        return reference[0]

    # 2. representative genome
    representative = [asm for asm in candidates if asm.get('RefSeq_category') == 'representative genome']
    if representative:
        return representative[0]

    # 3. Assemblies that have a RefSeq counterpart (FtpPath_RefSeq)
    with_refseq_path = [asm for asm in candidates if asm.get('FtpPath_RefSeq')]
    if with_refseq_path:
        return with_refseq_path[0]

    # 4. If genbank mode and we know there's a refseq reference-like assembly in the group,
    #    we ideally pick a GCA that matches that group. If we're here, it means we didn't find
    #    a direct counterpart, but we might still just pick from what we have.
    # (This step is now less critical because we're pre-grouping.)

    # 5. fallback: first available
    if candidates:
        return candidates[0]
    return None


def get_reference_genome(taxid: str, preference: str) -> Tuple[Optional[str], Optional[str]]:
    if not taxid.isdigit():
        logging.error(f"Invalid TaxID '{taxid}'. It should contain only digits.")
        return (None, None)

    try:
        search_term = f'txid{taxid}[Organism:exp]'
        handle = Entrez.esearch(db="assembly", term=search_term, retmode="xml", retmax=100)
        search_results = Entrez.read(handle)
        handle.close()

        if not search_results["IdList"]:
            logging.info(f"No assembly found for TaxID {taxid}")
            return (None, None)

        assembly_ids = search_results["IdList"]
        logging.debug(f"Found Assembly IDs {assembly_ids} for TaxID {taxid}")

        handle = Entrez.esummary(db="assembly", id=",".join(assembly_ids), retmode="xml")
        summaries = Entrez.read(handle)
        handle.close()

        if not summaries or 'DocumentSummarySet' not in summaries:
            logging.warning(f"No assembly records found for TaxID {taxid}")
            return (None, None)

        assembly_list = summaries['DocumentSummarySet']['DocumentSummary']

        # Exclude suppressed assemblies
        assembly_list = [asm for asm in assembly_list if asm.get('Status', '').lower() != 'suppressed']

        if not assembly_list:
            logging.info(f"No assemblies available (after excluding suppressed) for TaxID {taxid}.")
            return (None, None)

        # Group assemblies by numeric portion
        groups: Dict[str, Dict[str, List[Dict]]] = {}
        # Also track if we have a refseq reference/representative genome in a group
        refseq_reference_groups = set()

        for asm in assembly_list:
            acc = asm.get('AssemblyAccession', '')
            numeric_id = main_numeric_part(acc)
            if numeric_id not in groups:
                groups[numeric_id] = {'GCF': [], 'GCA': []}
            if acc.startswith('GCF_'):
                groups[numeric_id]['GCF'].append(asm)
                # Check if reference or representative
                if asm.get('RefSeq_category') in ['reference genome', 'representative genome']:
                    refseq_reference_groups.add(numeric_id)
            elif acc.startswith('GCA_'):
                groups[numeric_id]['GCA'].append(asm)

        # Now select based on preference
        # If preference is auto:
        #   - Consider all assemblies (both GCF and GCA) in each group and pick the best
        # If preference is refseq:
        #   - Only consider GCF assemblies in each group
        # If preference is genbank:
        #   - If we know there's a refseq reference-like assembly in a group, try to pick a GCA from that group.
        #     If that fails, pick any GCA from any group.

        selected = None

        if preference == 'auto':
            # Consider all assemblies together
            all_candidates = []
            for numeric_id, data in groups.items():
                all_candidates.extend(data['GCF'])
                all_candidates.extend(data['GCA'])
            selected = select_assembly(all_candidates, preference='auto', has_refseq_reference=False)

        elif preference == 'refseq':
            # Consider only GCF assemblies
            # Combine all GCF assemblies from all groups
            gcf_candidates = []
            for numeric_id, data in groups.items():
                gcf_candidates.extend(data['GCF'])
            if gcf_candidates:
                selected = select_assembly(gcf_candidates, preference='refseq', has_refseq_reference=False)

        elif preference == 'genbank':
            # If a group has a known refseq reference, try to pick GCA from that group first
            # to ensure we match a known refseq reference/rep genome
            # Prioritize groups that have refseq references
            prioritized_candidates = []
            fallback_candidates = []
            for numeric_id, data in groups.items():
                if numeric_id in refseq_reference_groups:
                    # Try from this group first
                    asm = select_assembly(data['GCA'], preference='genbank', has_refseq_reference=True)
                    if asm:
                        prioritized_candidates.append(asm)
                else:
                    # Will consider if we don't find any from a refseq reference group
                    asm = select_assembly(data['GCA'], preference='genbank', has_refseq_reference=False)
                    if asm:
                        fallback_candidates.append(asm)

            if prioritized_candidates:
                selected = prioritized_candidates[0]
            elif fallback_candidates:
                selected = fallback_candidates[0]

        if not selected:
            return (None, None)

        genbank_accession = selected.get('AssemblyAccession', 'Not Available')
        assembly_name = selected.get('AssemblyName', 'Not Available')
        logging.debug(f"Selected assembly for TaxID {taxid}: {genbank_accession}, {assembly_name}")
        return (genbank_accession, assembly_name)

    except Exception as e:
        logging.error(f"Error fetching assembly for TaxID '{taxid}': {e}")
        return (None, None)


def write_to_csv(output_file: str, data: List[Dict[str, str]], tax_levels: List[str]) -> None:
    try:
        fieldnames = ['Species Name', 'TaxID', 'Common Name', 'GenBank Accession', 'Assembly Name']
        taxonomy_columns = [level.capitalize() for level in tax_levels]
        fieldnames.extend(taxonomy_columns)

        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for entry in data:
                row = {
                    'Species Name': entry.get('Species Name', 'Not Available'),
                    'TaxID': entry.get('TaxID', 'Not Found'),
                    'Common Name': entry.get('Common Name', 'Not Available'),
                    'GenBank Accession': entry.get('GenBank Accession', 'Not Available'),
                    'Assembly Name': entry.get('Assembly Name', 'Not Available')
                }
                for level in taxonomy_columns:
                    row[level] = entry.get(level, 'Not Available')
                writer.writerow(row)
        logging.info(f"Data successfully written to {output_file}")
    except Exception as e:
        logging.error(f"Error writing to CSV file '{output_file}': {e}")
        sys.exit(1)


def main() -> None:
    args = parse_arguments()
    setup_logging(args.log_level)

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
            genbank_acc, assembly_name = get_reference_genome(taxid, args.assembly_preference)
            taxonomy_info = get_taxonomy_levels(taxid, args.tax_levels)
        else:
            common_name = None
            genbank_acc = None
            assembly_name = None
            taxonomy_info = {level.capitalize(): 'Not Available' for level in args.tax_levels}

        entry = {
            'Species Name': species,
            'TaxID': taxid if taxid else 'Not Found',
            'Common Name': common_name if common_name else 'Not Available',
            'GenBank Accession': genbank_acc if genbank_acc else 'Not Available',
            'Assembly Name': assembly_name if assembly_name else 'Not Available'
        }
        # Add taxonomy levels
        for level in args.tax_levels:
            key = level.capitalize()
            entry[key] = taxonomy_info.get(key, 'Not Available')

        results.append(entry)

        # Respect NCBI's rate limiting
        logging.debug(f"Sleeping for {args.delay} seconds to respect rate limits.")
        time.sleep(args.delay)

    write_to_csv(args.output, results, args.tax_levels)
    logging.info("Data extraction complete.")


if __name__ == "__main__":
    main()