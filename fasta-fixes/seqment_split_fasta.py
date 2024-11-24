#!/usr/bin/env python3

import argparse
import logging
import sys
import os
from typing import Tuple, List
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Split a FASTA file into two based on header string matching using Biopython."
    )
    parser.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="Path to the input FASTA file."
    )
    parser.add_argument(
        "-s", "--string",
        type=str,
        required=True,
        help="String to search for in the FASTA headers."
    )
    parser.add_argument(
        "-o1", "--output_match",
        type=str,
        default="matched.fasta",
        help="Output FASTA file for matching headers. Default: matched.fasta"
    )
    parser.add_argument(
        "-o2", "--output_non_match",
        type=str,
        default="non_matched.fasta",
        help="Output FASTA file for non-matching headers. Default: non_matched.fasta"
    )
    parser.add_argument(
        "-c", "--case_sensitive",
        action="store_true",
        help="Enable case-sensitive matching. Default is case-insensitive."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging."
    )
    return parser.parse_args()

def setup_logging(verbose: bool) -> None:
    """
    Configure logging settings.

    Args:
        verbose (bool): If True, set logging level to DEBUG.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        stream=sys.stdout
    )

def confirm_overwrite(files: List[str]) -> bool:
    """
    Check if any of the specified files exist and prompt the user for overwrite confirmation.

    Args:
        files (List[str]): List of file paths to check.

    Returns:
        bool: True if it's safe to proceed (either files don't exist or user confirmed overwrite), False otherwise.
    """
    existing_files = [f for f in files if os.path.exists(f)]
    if not existing_files:
        logging.debug("No existing output files detected.")
        return True  # No files exist, safe to proceed

    logging.warning("The following output file(s) already exist:")
    for f in existing_files:
        logging.warning(f"  - {f}")

    while True:
        response = input("Do you want to overwrite them? [y/N]: ").strip().lower()
        if response in ('y', 'yes'):
            logging.info("User confirmed to overwrite existing files.")
            return True
        elif response in ('n', 'no', ''):
            logging.info("User declined to overwrite existing files. Exiting.")
            return False
        else:
            print("Please respond with 'y' or 'n'.")

def open_output_files(output_match: str, output_non_match: str) -> Tuple[str, str]:
    """
    Validate output file paths and confirm overwriting if necessary.

    Args:
        output_match (str): Path for matched output.
        output_non_match (str): Path for non-matched output.

    Returns:
        Tuple[str, str]: Tuple containing output paths.
    """
    files_to_check = [output_match, output_non_match]
    if not confirm_overwrite(files_to_check):
        sys.exit(0)  # Exit gracefully if user does not want to overwrite

    logging.debug(f"Output match file: {output_match}")
    logging.debug(f"Output non-match file: {output_non_match}")
    return output_match, output_non_match

def validate_fasta(input_path: str) -> bool:
    """
    Validate that the input file is in FASTA format.

    Args:
        input_path (str): Path to the input FASTA file.

    Returns:
        bool: True if valid FASTA format, False otherwise.
    """
    try:
        with open(input_path, 'r') as infile:
            records = list(SeqIO.parse(infile, 'fasta'))
            if not records:
                logging.error("No FASTA records found in the input file.")
                return False
        logging.info("Input file is a valid FASTA format.")
        return True
    except Exception as e:
        logging.error(f"Error validating FASTA format: {e}")
        return False

def split_fasta(
    input_path: str,
    search_string: str,
    output_match: str,
    output_non_match: str,
    case_sensitive: bool
) -> None:
    """
    Split the FASTA file into matched and non-matched based on header.

    Args:
        input_path (str): Path to the input FASTA file.
        search_string (str): String to search for in headers.
        output_match (str): Path for matched sequences.
        output_non_match (str): Path for non-matched sequences.
        case_sensitive (bool): Flag for case-sensitive matching.
    """
    try:
        logging.info(f"Reading input FASTA file: {input_path}")
        matched_records: List[SeqRecord] = []
        non_matched_records: List[SeqRecord] = []

        with open(input_path, 'r') as infile:
            for record in SeqIO.parse(infile, 'fasta'):
                header = record.description
                if case_sensitive:
                    match = search_string in header
                else:
                    match = search_string.lower() in header.lower()

                if match:
                    matched_records.append(record)
                    logging.debug(f"Matched: {header}")
                else:
                    non_matched_records.append(record)
                    logging.debug(f"Non-matched: {header}")

        logging.info(f"Writing {len(matched_records)} matched records to {output_match}")
        with open(output_match, 'w') as matched_out:
            SeqIO.write(matched_records, matched_out, 'fasta')

        logging.info(f"Writing {len(non_matched_records)} non-matched records to {output_non_match}")
        with open(output_non_match, 'w') as non_matched_out:
            SeqIO.write(non_matched_records, non_matched_out, 'fasta')

        logging.info("FASTA splitting completed successfully.")

    except FileNotFoundError:
        logging.error(f"Input file not found: {input_path}")
        sys.exit(1)
    except IOError as e:
        logging.error(f"I/O error({e.errno}): {e.strerror}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        sys.exit(1)

def main() -> None:
    """
    Main function to orchestrate the FASTA splitting process.
    """
    args = parse_arguments()
    setup_logging(args.verbose)

    # Validate input FASTA format
    if not validate_fasta(args.input):
        logging.error("Input file validation failed. Please provide a valid FASTA file.")
        sys.exit(1)

    output_match, output_non_match = open_output_files(args.output_match, args.output_non_match)
    split_fasta(
        input_path=args.input,
        search_string=args.string,
        output_match=output_match,
        output_non_match=output_non_match,
        case_sensitive=args.case_sensitive
    )

if __name__ == "__main__":
    main()