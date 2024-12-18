# README

## Species Taxonomy and Assembly Data Retrieval Scripts

This repository contains two Python scripts designed to help you fetch taxonomy information and genome assembly details for a list of species and prepare the data for use with the [`sanger-tol/insdcdownload`](https://github.com/sanger-tol/insdcdownload) tool.

### Scripts Included:

1. **`fetch_taxonomy.py`**: Retrieves NCBI TaxIDs, common names, and genome assembly information (including GenBank accession numbers and assembly names) for a list of species.

2. **`process_species_assemblies.py`**: Processes the output from `fetch_taxonomy.py` to generate a file suitable for use with the `sanger-tol/insdcdownload` tool.

---

## Prerequisites

- **Python 3.6 or later**: Ensure you have Python installed on your system.
- **Biopython**: This library is required for interacting with NCBI's Entrez API.
- **NCBI API Key (optional but recommended)**: Increases your request rate limit when accessing NCBI's Entrez API.

## Installation

1. **Clone the Repository** (if applicable):

   ```bash
   git clone https://github.com/druvus/misc_scripts.git
   cd misc_scripts
   ```

2. **Install Biopython**:

   ```bash
   pip install biopython
   ```

3. **Download the Scripts**:

   Ensure you have both `fetch_taxonomy.py` and `process_species_assemblies.py` in your working directory.

---

## Usage

### 1. `fetch_taxonomy.py`

This script fetches NCBI TaxIDs, common names, and genome assembly information for a list of species.

#### **Command-Line Arguments**

- `-i` or `--input`: Path to the input text file containing species names (one per line).
- `-o` or `--output`: Path to the output CSV file. Default is `species_taxonomy.csv`.
- `-e` or `--email`: Your email address (required by NCBI Entrez).
- `-k` or `--api_key`: Your NCBI API key (optional but recommended).
- `-l` or `--log_level`: Set the logging level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). Default is `INFO`.
- `-d` or `--delay`: Delay between API requests in seconds. Default is `0.4`.

#### **Example Usage**

```bash
python fetch_taxonomy.py -i species_names.txt -o species_taxonomy.csv -e "your_email@example.com" -k YOUR_NCBI_API_KEY
```

#### **Steps**

1. **Prepare the Input File** (`species_names.txt`):

   Create a text file with one species name per line. Example:

   ```
   Homo sapiens
   Escherichia coli
   Mus musculus
   ```

2. **Run the Script**:

   ```bash
   python fetch_taxonomy.py -i species_names.txt -o species_taxonomy.csv -e "your_email@example.com" -k YOUR_NCBI_API_KEY
   ```

3. **Check the Output**:

   After execution, `species_taxonomy.csv` will contain the taxonomy and assembly details.

### 2. `process_species_assemblies.py`

This script processes the output from `fetch_taxonomy.py` to generate a file suitable for use with `sanger-tol/insdcdownload`.

#### **Command-Line Arguments**

- `-i` or `--input`: Path to the input CSV file generated by `fetch_taxonomy.py`.
- `-o` or `--output`: Path to the output file.
- `-p` or `--prefix`: Prefix to be combined with species names to form the output directory (`outdir`).
- `-l` or `--log_level`: Set the logging level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). Default is `INFO`.

#### **Example Usage**

```bash
python process_species_assemblies.py -i species_taxonomy.csv -o insdc_download_input.csv -p darwin/data/mammals
```

#### **Steps**

1. **Ensure You Have the Output from `fetch_taxonomy.py`**:

   The input file (`species_taxonomy.csv`) should be generated from the previous script.

2. **Run the Script**:

   ```bash
   python process_species_assemblies.py -i species_taxonomy.csv -o insdc_download_input.csv -p darwin/data/mammals
   ```

3. **Check the Output**:

   The script will generate `insdc_download_input.csv` with columns:

   - `outdir`: Combined prefix and species name (spaces replaced by underscores).
   - `assembly_name`: The assembly name from NCBI.
   - `assembly_accession`: The GenBank accession number.

---

## Output Files

### 1. `species_taxonomy.csv` (from `fetch_taxonomy.py`)

Contains the following columns:

- `Species Name`
- `TaxID`
- `Common Name`
- `GenBank Accession`
- `Assembly Name`

**Sample Content**:

```csv
Species Name,TaxID,Common Name,GenBank Accession,Assembly Name
Homo sapiens,9606,human,GCA_000001405.28,GRCh38.p13
Escherichia coli,562,E. coli,GCA_000005845.2,ASM584v2
Mus musculus,10090,house mouse,GCF_000001635.27,GRCm39
```

### 2. `insdc_download_input.csv` (from `process_species_assemblies.py`)

Contains the following columns:

- `outdir`
- `assembly_name`
- `assembly_accession`

**Sample Content**:

```csv
outdir,assembly_name,assembly_accession
darwin/data/mammals/Homo_sapiens,GRCh38.p13,GCA_000001405.28
darwin/data/bacteria/Escherichia_coli,ASM584v2,GCA_000005845.2
darwin/data/mammals/Mus_musculus,GRCm39,GCF_000001635.27
```

---

## Examples

### Example 1: Fetching Taxonomy and Assembly Data

```bash
python fetch_taxonomy.py -i species_list.txt -o species_taxonomy.csv -e "user@example.com" -k YOUR_NCBI_API_KEY
```

- **species_list.txt**:

  ```
  Homo sapiens
  Escherichia coli
  Mus musculus
  ```

- **Command Explanation**:

  - `-i species_list.txt`: Input file containing species names.
  - `-o species_taxonomy.csv`: Output CSV file.
  - `-e "user@example.com"`: Your email address.
  - `-k YOUR_NCBI_API_KEY`: Your NCBI API key.

### Example 2: Processing Assemblies for `insdcdownload`

```bash
python process_species_assemblies.py -i species_taxonomy.csv -o insdc_download_input.csv -p darwin/data
```

- **Command Explanation**:

  - `-i species_taxonomy.csv`: Input CSV file from `fetch_taxonomy.py`.
  - `-o insdc_download_input.csv`: Output file for `insdcdownload`.
  - `-p darwin/data`: Prefix for the `outdir`.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
