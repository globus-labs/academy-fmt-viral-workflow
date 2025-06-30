import os
import subprocess
import shutil
import parsl
from parsl.app.app import python_app, bash_app
from parsl.configs.local_threads import config
from parsl.data_provider.files import File

print(parsl.__version__)

parsl.load(config)

# === Turn config file into a dictionary of variables ===

@python_app
def make_config(config_file):
    config = {}
    with open(config_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('export '):
                line = line[len('export '):]  # Remove 'export '
            if '=' in line:
                key, val = line.split('=', 1)
                config[key.strip()] = val.strip().strip('"').strip("'")
    return config

@python_app
def read_sample_ids(sample_ids_file):
    # Read sample IDs from the file
    with open(sample_ids_file, "r") as f:
        sample_ids = [line.strip() for line in f if line.strip()]
    return sample_ids


@python_app
def unzip_fasta(spades_gz, unzipped_spades_path):

    """
    Check and unzip the FASTA files for a list of sample IDs using subprocess.
    """
    import subprocess
    import os
    from parsl.data_provider.files import File

    # Skip if unzipped file already exists
    if os.path.exists(unzipped_spades_path):
        print(f"[INFO] File already exists: {unzipped_spades}")
        return

    # Attempt to unzip the file
    try:
        subprocess.run(["gzip", "-d", spades_gz], check=True)
        print(f"[INFO] Successfully unzipped: {spades_gz}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to unzip {spades_gz}: {e}")

    unzipped_spades = File(unzipped_spades_path)
    return unzipped_spades

# === GeNomad ===
@python_app
def run_genomad(unzipped_spades,genomad_output_dir,db):
    """
    Runs Genomad end-to-end for a list of sample IDs using subprocess.Popen for parallel execution.
 """
    import subprocess
    import os

    # === Genomad command ====
    cmd = [
        "conda", "run", "-n", "genomad_env",
        "genomad", "end-to-end", "--cleanup","--restart",
        unzipped_spades, genomad_output_dir, db
    ]
    
    # Run the command
    subprocess.run(cmd, check=True)

    genomad_virus = os.path.join(genomad_output_dir, "contigs_summary", "contigs_virus.fna")
    return genomad_virus

# === CheckV on GeNomad ===

@python_app
def run_checkv_genomad(checkv_parser, parse_length, work_dir, unzipped_spades, genomad_virus, checkv_output_dir, parse_input, selection_csv,checkvdb):

    """
    Runs CheckV and post-processing steps for multiple samples using subprocess.Popen for parallel execution.
    """
    import subprocess
    import os

    # Ensure output directory exists
    os.makedirs(checkv_output_dir, exist_ok=True)

    # === CheckV Command ===
    cmd_checkv = [
        "conda", "run", "-n", "checkv_env",
        "checkv", "end_to_end", genomad_virus, checkv_output_dir, "-t", "4",
        "-d", checkvdb
    ]

    # === Parse output with R Script ===
    cmd_parser = [
        "conda", "run", "-n", "r_env", 
        "Rscript", checkv_parser,    
        "-i", parse_input,
        "-l", parse_length,
        "-o", selection_csv            
    ]

    # === Subset FASTA with seqtk ===
    cmd_seqtk = [
        "conda", "run", "-n", "seqtk_env",
        "seqtk", "subseq", unzipped_spades, selection_csv
    ]

    # Run the commands
    subprocess.run(cmd_checkv, check=True)
    subprocess.run(cmd_parser, check=True)
    subprocess.run(cmd_seqtk, check=True)

    # === Return to working directory ===
    os.chdir(work_dir)
    print(f"[Done] Returned to working directory: {work_dir}")
    
    subset_spades = os.path.join(checkv_output_dir, "subset_spades.fasta")
    return subset_spades

# === Dereplication ===

@python_app
def run_dereplicate(subset_spades, cluster_res, tmp_dir, input_fasta, cleaned_fasta, out_derep):

    """
    Runs mmseqs2 dereplication and awk processing for multiple samples in parallel using subprocess.Popen.
    """
    import subprocess
    import os

    # === Create output directory ===
    os.makedirs(cluster_dir, exist_ok=True)

    # === Run mmseqs easy-cluster for dereplication ===       
    print(f"[mmseqs] Clustering {subset_spades} for sample {sample_id}")
    cmd_mmseqs = [
        "conda", "run", "-n", "mmseqs2_env",
        "mmseqs", "easy-cluster", subset_spades, cluster_res, tmp_dir,
        "--min-seq-id", "0.99", "-c", "0.90", "--cov-mode", "1"
    ]
        
    subprocess.run(cmd_mmseqs, check=True)

    print(f"[awk] Removing duplicate headers for {input_fasta} in sample {sample_id}")
    cmd_awk = (
        r"""awk '/^>/{if($0!=prev){print; prev=$0}} !/^>/' """
        + input_fasta +
        f" > {cleaned_fasta}"
    )

    subprocess.run(cmd_awk, check=True)

    derep_fasta = os.path.join(out_derep, "dereplicated.fasta")
    return derep_fasta 
   
# === Clustering ===

@python_app
def run_cluster_all(out_derep, out_cluster, derep_fasta, cluster_res, tmp_dir, rep_seq_src, rep_seq_dst):

    """
    Clusters all dereplicated sequences using mmseqs2 in a conda environment.
    """
    import subprocess
    import os 
    # === Combine all dereplicated FASTA files ===

    print(f"[combine] Writing all cleaned_clusterRes_all_seqs.fasta files to {derep_fasta}")
    with open(derep_fasta, 'w') as outfile:
        for root, _, files in os.walk(out_derep):
            for file in files:
                if file.endswith("cleaned_clusterRes_all_seqs.fasta"):
                    file_path = os.path.join(root, file)
                    print(f" - Adding: {file_path}")
                    with open(file_path) as infile:
                        shutil.copyfileobj(infile, outfile)

    # === Run mmseqs easy-cluster for clustering ===

    print(f"[mmseqs] Running clustering on {derep_fasta}")
    os.makedirs(out_cluster, exist_ok=True)
    subprocess.run([
        "conda", "run", "-n", "mmseqs2_env",
        "mmseqs", "easy-cluster", derep_fasta, cluster_res, tmp_dir,
        "--min-seq-id", "0.95", "-c", "0.75", "--cov-mode", "1"
    ], check=True)

    # === Copy representative sequences to query folder ===

    print(f"[copy] Copying {rep_seq_src} to {rep_seq_dst}")
    os.makedirs(rep_seq_dst, exist_ok=True)
    shutil.copy(rep_seq_src, os.path.join(rep_seq_dst, "clusterRes_rep_seq.fasta"))

    print(f"[done] Clustering and copy complete.")

    query_dir = rep_seq_dst
    return query_dir

# === Make BLAST DB ===

@python_app
def make_blast_db(db_dir, max_db_size, db_list_path):

    """
    Scans DB_DIR for all .fasta files and builds BLAST databases using apptainer and makeblastdb.
    """
    import subprocess
    import os 

    print(f"Working in: {db_dir}")
    os.chdir(db_dir)

    # === Generate db-list file using os.walk ===

    with open(db_list_path, "w") as db_list:
        for root, dirs, files in os.walk("."):
            for file in files:
                if file.endswith(".fasta"):
                    rel_path = os.path.join(root, file).lstrip("./")
                    db_list.write(rel_path + "\n")

    # === Check if db-list exists and is non-empty ===

    if not os.path.exists(db_list_path) or os.path.getsize(db_list_path) == 0:
        print(f"Cannot find or empty database list: {db_list_path}")
        return

    # === Process each FASTA file ===

    with open(db_list_path) as f:
        for i, line in enumerate(f, start=1):
            db_file = line.strip()
            db_name = os.path.basename(db_file)

            print(f"{i:5d}: {db_name}")

            # Build and run makeblastdb command via apptainer

            cmd = [
                "conda", "run", "-n", "blast_env",
                "makeblastdb",
                "-title", db_name,
                "-out", db_name,
                "-in", db_file,
                "-dbtype", "nucl",
                "-max_file_sz", str(max_db_size)
            ]

            subprocess.run(cmd, check=True)

    # === Print finish time ===

    result = subprocess.run(["date"], capture_output=True, text=True)
    print(f"Finished {result.stdout.strip()}")

# === Launch Blast === 

@python_app
def run_launch_blast(split_size, results_dir, files_list_path, query_dir, db_dir, blast_results_dir, blast_type, eval_param, out_fmt, max_target_seqs, merge_results_dir):

    """
    Splits FASTA files, runs BLAST per split file, and merges results using Python functions.
    """
    import subprocess
    import os 
    # === Reinitialize results and query directory ===
    if os.path.exists(results_dir):
        shutil.rmtree(results_dir)
    os.makedirs(results_dir)
    
    # === List all .fasta files in FASTA_DIR ===
    with open(files_list_path, "w") as file_list:
        for root, dirs, files in os.walk(query_dir):
            for file in files:
                if file.endswith(".fasta"):
                    full_path = os.path.join(root, file)
                    file_list.write(os.path.relpath(full_path, query_dir) + "\n")

    with open(files_list_path) as f:
        fasta_files = [line.strip() for line in f.readlines()]

    if not fasta_files:
        print("No FASTA files found.")
        return

    print(f"Found {len(fasta_files)} FASTA files.")

    for i, file_rel in enumerate(fasta_files, 1):
        file_path = os.path.join(query_dir, file_rel)
        file_name = os.path.basename(file_path)
        print(f"\n[{i}] Processing {file_name}")

        out_dir = os.path.join(results_dir, file_name)
        split_dir = os.path.join(query_dir, "fa_split")

        # === Reset output directory ===

        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.makedirs(split_dir)

        # === Split the FASTA file using faSplit ===

        print(f"Splitting {file_name} into chunks of size {split_size}")
 
        cmd = [
            "conda", "run", "-n", "fasplit_env",
            "faSplit", "about", file_path, str(split_size), f"{split_dir}/"
        ]

        subprocess.run(cmd, check=True)

        # === List all split files ===

        split_files = sorted([
            f for f in os.listdir(split_dir) if f.endswith(".fa")
        ])
        num_split = len(split_files)
        print(f"{num_split} split files created")

        # === Save list of split files for BLAST step ===

        split_files_list = os.path.join(split_dir, "split_files_list.txt")
        with open(split_files_list, "w") as f:
            for sfile in split_files:
                f.write(sfile + "\n")

        # === Run BLAST on each split file ===
        for task_id in range(num_split):
            print(f"Running BLAST for split {task_id}")
            db_list_path = run_blast(task_id, db_dir, split_files_list, split_dir, blast_results_dir, blast_type, eval_param, out_fmt, max_target_seqs)

        # === Merge results to GFF ===
        hits_file = merge_blast(merge_results_dir, db_list_path, file_name)

    print("\n All BLAST jobs completed.")
    return hits_file


# === Blast ===

@python_app
def run_blast(task_id, db_dir, split_files_list, split_dir, results_dir, blast_type, eval_param, out_fmt, max_target_seqs):

    """
    Runs the BLAST command for each split file against each database.
    """
    import subprocess
    import os 
    # Read the SPLIT_FILES_LIST to get the split file names
    with open(split_files_list, 'r') as f:
        split_files = [line.strip() for line in f.readlines()]

    # Select the split file based on the task_id
    split_file = split_files[task_id]

    # Create results directory if it doesn't exist
    os.makedirs(results_dir, exist_ok=True)

    # === Iterate over databases ===
    db_list_path = os.path.join(db_dir, "db-list")
    with open(db_list_path, 'r') as f:
        databases = [line.strip() for line in f.readlines()]

    # Iterate over each database
    for i, db in enumerate(databases, 1):
        print(f"{i:5d}: {db}")

        # Prepare output directories and paths
        results_by_db = os.path.join(results_dir, db, os.path.basename(split_file))
        os.makedirs(results_by_db, exist_ok=True)

        blast_out = os.path.join(results_by_db, f"{split_file}.blastout")
        blast_db = os.path.join(db_dir, db)

        # Run BLAST command
        blast_cmd = [
            "conda", "run", "-n", "blast_env", blast_type,
            "-num_threads", "48",
            "-db", blast_db,
            "-query", os.path.join(split_dir, split_file),
            "-out", blast_out,
            "-evalue", str(eval_param),
            "-outfmt", str(out_fmt),
            "-max_target_seqs", str(max_target_seqs)
        ]
        # Execute the BLAST command
        subprocess.run(blast_cmd, check=True)

    print("Finished BLAST processing.")
    return db_list_path

# === Merge Blast ===

@python_app
def merge_blast(results_dir, db_list_path, file_name):

    """
    Merges BLAST results and converts them to GFF format.
    """
    import subprocess
    import os 
    # Create the results directory (remove prior runs and create new)
    os.makedirs(merge_results_dir, exist_ok=True)

    # === Read the db-list and process each database ===
    with open(db_list_path, 'r') as f:
        databases = [line.strip() for line in f.readlines()]

    # Iterate over each database
    for i, db in enumerate(databases, 1):
        print(f"{i:5d}: {db}")

        # Prepare paths for results
        results_by_db = os.path.join(results_dir, db)
        os.makedirs(results_by_db, exist_ok=True)

        blast_results = os.path.join(results_by_db, f"{file_name}.txt")
        blast_gff = os.path.join(results_by_db, f"{file_name}.gff")
        blast_out = os.path.join(work_dir, "results", "05C_blast", db, file_name)

        # Merge BLAST results into a single file
        print(f"Merging BLAST results for {db}")
        with open(blast_results, 'w') as outfile:
            for result_file in os.listdir(blast_out):
                result_path = os.path.join(blast_out, result_file)
                with open(result_path, 'r') as infile:
                    outfile.write(infile.read())

        # Convert to GFF format using Python
        print(f"Converting BLAST results to GFF format for {db}")
        with open(blast_results, 'r') as infile, open(blast_gff, 'w') as outfile:
            for line in infile:
                fields = line.strip().split('\t')
                if len(fields) > 7:  # Ensure valid BLAST output line
                    gff_line = f"{fields[0]}\tblast\tgene\t{fields[6]}\t{fields[7]}\t.\t.\t.\tID=Gene{fields[6]};Name={fields[1]}\n"
                    outfile.write(gff_line)
        
    print("Finished processing BLAST results.")
    hits_file = os.path.join(results_dir, "AVrC_allrepresentativesfasta", "clusterRes_rep_seq.fasta.txt")
    return hits_file

# === Annotation ===

@python_app
def annotate_blast(hits_file, annotations_dir, output_dir, script_path, pctid, length):

    """
    Loops over all annotation CSV files and runs the annotation script using subprocess.run.
    Replicates shell logic: looping over files, constructing arguments, and running script.
    """
    import subprocess
    import os
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Remove unwanted files
    for file in os.listdir(annotations_dir):
        if file == '.DS_Store' or file.startswith('._'):
            os.remove(os.path.join(annotations_dir, file))

    # Debug print of annotation files
    print("Files in annotation directory:")
    for f in os.listdir(annotations_dir):
        print(os.path.join(annotations_dir, f))

    # Loop over each annotation file
    for file in os.listdir(annotations_dir):
        ann_path = os.path.join(annotations_dir, file)
        if os.path.isfile(ann_path):
            out_file = os.path.join(output_dir, f"annotated_{file}")
            cmd = [
                script_path,
                "-b", hits_file,
                "-a", ann_path,
                "-o", out_file,
                "-p", pctid,
                "-l", length
            ]
            print("Running:", " ".join(cmd))
            subprocess.run(cmd, check=True)


def main():
    # === Load configuration ===
    config_path = os.path.join(os.getcwd(), "config_py.sh")
    config_file = File(config_path)
    config = make_config(config_file)

    sample_ids_file = File(os.path.join(config['XFILE_DIR'], config['XFILE']))
    sample_ids = read_sample_ids(sample_ids_file)
    
    for sample_id in sample_ids:
        # === Define variables for unzipping ===
        spades_gz = File(os.path.join(config['SPADES_DIR'], sample_id, "contigs.fasta.gz"))
        unzipped_spades_path = os.path.join(config['SPADES_DIR'], sample_i
d, "contigs.fasta")
        # === Unzip ===
        unzipped_spades = unzip_fasta(spades_gz, unzipped_spades_path)

        # === Define variables for Genomad ===
        genomad_output_dir = os.path.join(config['OUT_GENOMAD'], sample_id)
        db = config["GENOMAD_DB"]
        # === Run Genomad ===
        genomad_virus = run_genomad(unzipped_spades,genomad_output_dir,db)
        
        # === Define variables for CheckV ===
        checkv_parser = config["CHECKV_PARSER"]
        parse_length = str(config["PARSE_LENGTH"])
        work_dir = config["WORK_DIR"]
        checkv_output_dir = os.path.join(config["OUT_CHECKV_GENOMAD"], sample_id)
        parse_input = os.path.join(checkv_output_dir, "contamination.tsv")
        selection_csv = os.path.join(checkv_output_dir, "selection2_viral.csv")
        checkvdb = config["CHECKVDB"] 
        # === Run CheckV ===
        subset_spades = run_checkv_genomad(checkv_parser, parse_length, work_dir, unzipped_spades, genomad_virus, checkv_output_dir, parse_input, selection_csv,checkvdb)

        # === Define variables for Dereplication ===
        cluster_res = os.path.join(config["OUT_DEREP"], sample_id, "clusterRes")
        tmp_dir = os.path.join(config["OUT_DEREP"], sample_id, "tmp")
        input_fasta = f"{cluster_res}_all_seqs.fasta"
        cleaned_fasta = os.path.join(config["OUT_DEREP"], sample_id, "cleaned_clusterRes_all_seqs.fasta")
        out_derep = config["OUT_DEREP"]
        # === Dereplicate ===
        derep_fasta = run_dereplicate(subset_spades, cluster_res, tmp_dir, input_fasta, cleaned_fasta, out_derep)

    # === Define variables for clustering ===
    out_cluster = config["OUT_CLUSTER"]
    work_dir = config["WORK_DIR"]
    cluster_res = os.path.join(out_cluster, "clusterRes")
    tmp_dir = os.path.join(out_cluster, "tmp")
    rep_seq_src = os.path.join(out_cluster, "clusterRes_rep_seq.fasta")
    rep_seq_dst = os.path.join(work_dir, "query")
    # === Cluster ===
    query_dir = run_cluster_all(out_derep, out_cluster, derep_fasta, cluster_res, tmp_dir, rep_seq_src, rep_seq_dst)

    # === Define variables to Make BLASTDB ===
    db_dir = config["DB_DIR"]
    max_db_size = config["MAX_DB_SIZE"]
    db_list_path = os.path.join(db_dir, "db-list")
    # Make BLASTDB
    make_blast_db(db_dir, max_db_size, db_list_path)

    # === Define BLAST variables ===
    prog = "05B_launchblast"
    fasta_dir = config["FASTA_DIR"]
    split_size = config["FA_SPLIT_FILE_SIZE"]
    results_dir = os.path.join(work_dir, "results_testing", prog)
    files_list_path = os.path.join(fasta_dir, "fasta-files")
         # RUN BLAST
    blast_results_dir = os.path.join(work_dir, "results_testing", "05C_blast")
    blast_type = config["BLAST_TYPE"]
    eval_param = config["EVAL"]
    out_fmt = config["OUT_FMT"]
    max_target_seqs = config["MAX_TARGET_SEQS"]
         # MERGE BLAST 
    merge_results_dir = os.path.join(work_dir, "results_testing", "05D_mergeblast")
    # === Launch BLAST ===
    hits_file = run_launch_blast(split_size, results_dir, files_list_path, query_dir, db_dir, blast_results_dir, blast_type, eval_param, out_fmt, max_target_seqs, merge_results_dir)

    # === Define variables for annotation ===
    annotations_dir = config['ANNOTATIONS']
    out_annotate = config['OUTPUT']
    script_path = os.path.join(config['SCRIPT_DIR'], "solution1_manual.py")
    pctid = config['PCTID']
    length = config['LENGTH']
    # === Annotate ===
    annotate_blast(hits_file, annotations_dir, out_annotate, script_path, pctid, length)


if __name__ == "__main__":
    main()

