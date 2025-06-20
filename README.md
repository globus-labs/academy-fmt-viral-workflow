# academy-fmt-viral-workflow
FMT virus detection workflow in Academy -- Naomi Summer 2025 Project

# Downloading All Databases and Environments Needed

## Genomad
- Github: https://github.com/apcamargo/genomad/blob/main/README.md

### Conda Install
    conda create -n genomad_env -c conda-forge -c bioconda genomad

### Testing
    conda activate genomad_env
    genomad --help

### Database
    genomad download-database /path/to/databases/directory

---

## CheckV
- Bitbucket: https://bitbucket.org/berkeleylab/checkv/src/master/README.md#markdown-header-installation
- Database portal: https://portal.nersc.gov/CheckV/

### Database
    wget -P /path/to/databases/directory https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
    tar -xzvf checkv-db-v1.5.tar.gz

### Conda Install
    conda create -n checkv_env -c conda-forge -c bioconda checkv=1.0.1 -y

### Testing
    conda activate checkv_env
    checkv --help

---

## AVrC Database
- Github: https://github.com/aponsero/Aggregated_Viral_Catalogue/blob/main/README.md
- Database: https://zenodo.org/records/11426065

### Download
    wget -O /xdisk/bhurwitz/databases/AVrC_allrepresentatives.fasta.gz "https://zenodo.org/records/11426065/files/AVrC_allrepresentatives.fasta.gz?download=1"
    gunzip AVrC_allrepresentatives.fasta.gz
Make a database directory (titled `AVrC`) and put the fasta file inside that directory.

---

## Mmseqs
- Github: https://github.com/soedinglab/MMseqs2

### Conda Install
    conda create -n mmseqs2_env -c bioconda mmseqs2=13.45111

---

## Seqtk

### Conda Install
    conda create -n seqtk_env -c bioconda seqtk

---

## FaSplit
- https://anaconda.org/bioconda/ucsc-fasplit

### Conda Install
    conda create -n fasplit_env -c bioconda ucsc-fasplit

### Testing
    conda activate fasplit_env
    faSplit

---

## Blast
- https://bioconda.github.io/recipes/blast/README.html

### Conda Install
    conda create --name blast_env -c bioconda blast

### Testing
    conda activate blast_env
    blastn -version
