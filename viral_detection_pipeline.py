import os
import subprocess
import shutil

# === Turn config file into a dictionary of variables ===

def make_config(config_path="./config_py.sh"):
    config = {}
    with open(config_path) as f:
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

def unzip_fasta(sample_ids, config):

    """
    Check and unzip the FASTA files for a list of sample IDs using subprocess.
    """
    for sample_id in sample_ids:
        spades_gz = os.path.join(config['SPADES_DIR'], sample_id, "contigs.fasta.gz")
        unzipped_spades = os.path.join(config['SPADES_DIR'], sample_id, "contigs.fasta")

        # Skip if unzipped file already exists
        if os.path.exists(unzipped_spades):
            print(f"[INFO] File already exists: {unzipped_spades}")
            continue

        # If .gz file doesn't exist, issue a warning
        if not os.path.exists(spades_gz):
            print(f"[WARNING] Gzipped file not found, assuming already unzipped: {spades_gz}")
            continue

        # Attempt to unzip the file
        try:
            subprocess.run(["gzip", "-d", spades_gz], check=True)
            print(f"[INFO] Successfully unzipped: {spades_gz}")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Failed to unzip {spades_gz}: {e}")

# === GeNomad ===

def run_genomad(sample_ids, config):

    """
    Runs Genomad end-to-end for a list of sample IDs using subprocess.Popen for parallel execution.

    Args:
        sample_ids (list): List of sample ID strings to process.
        config (dict): Configuration dict with 'SPADES_DIR', 'OUT_GENOMAD', 'GENOMAD_DB'.
 """
    processes = []

    for sample_id in sample_ids:
        # === Define paths from config ===
        unzipped_spades = os.path.join(config['SPADES_DIR'], sample_id, "contigs.fasta")
        output_dir = os.path.join(config['OUT_GENOMAD'], sample_id)
        db = config["GENOMAD_DB"]

        # === Genomad command ====
        cmd = [
            "conda", "run", "-n", "genomad_env",
            "genomad", "end-to-end", "--cleanup","--restart",
            unzipped_spades, output_dir, db
        ]

        # Run the command
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        processes.append((sample_id, proc))

    # Wait for all processes to complete
    for sample_id, proc in processes:
        stdout, stderr = proc.communicate()

        if proc.returncode == 0:
            print(f"Genomad {sample_id} completed successfully.")
        else:
            print(f"Genomad {sample_id} failed with return code {proc.returncode}.")
            print(stderr.decode())

    print("All Genomad jobs finished.")


# === CheckV on GeNomad ===

def run_checkv_genomad(sample_ids, config):

    """
    Runs CheckV and post-processing steps for multiple samples using subprocess.Popen for parallel execution.
    """
    processes = []
    for sample_id in sample_ids:
        # === Define paths from config ===
        spades_dir = config["SPADES_DIR"]
        out_genomad = config["OUT_GENOMAD"]
        out_checkv = config["OUT_CHECKV_GENOMAD"]
        checkv_parser = config["CHECKV_PARSER"]
        parse_length = str(config["PARSE_LENGTH"])  
        work_dir = config["WORK_DIR"]
        unzipped_spades = f"{spades_dir}/{sample_id}/contigs.fasta"
        genomad_input = f"{out_genomad}/{sample_id}/contigs_summary/contigs_virus.fna"
        output_dir = f"{out_checkv}/{sample_id}"
        parse_input = f"{output_dir}/contamination.tsv"
        selection_csv = f"{output_dir}/selection2_viral.csv"

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # === CheckV Command ===
        cmd_checkv = [
            "conda", "run", "-n", "checkv_env",
            "checkv", "end_to_end", genomad_input, output_dir, "-t", "4",
            "-d", config["CHECKVDB"]
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
        proc_checkv = subprocess.Popen(cmd_checkv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        processes.append(("CheckV", sample_id, proc_checkv, None))

        # Wait for CheckV to finish before proceeding
        stdout, stderr = proc_checkv.communicate() 
        if proc_checkv.returncode == 0:
            print(f"[CheckV] {sample_id} completed successfully.")
        else:
            print(f"[CheckV] {sample_id} failed.")
            print(stderr.decode())
            continue  

        proc_parser = subprocess.Popen(cmd_parser, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output_dir)
        processes.append(("Parser", sample_id, proc_parser, None))

        # Wait for Parser to finish before proceeding
        stdout, stderr = proc_parser.communicate()  
        if proc_parser.returncode == 0:
            print(f"[Parser] {sample_id} completed successfully.")
        else:
            print(f"[Parser] {sample_id} failed.")
            print(stderr.decode())
            continue  
        
        output_file_path = os.path.join(output_dir, "subset_spades.fasta")
        output_file = open(output_file_path, "w")
        proc_seqtk = subprocess.Popen(cmd_seqtk, stdout=output_file, stderr=subprocess.PIPE, cwd=output_dir)
        processes.append(("seqtk", sample_id, proc_seqtk, output_file))

    # === Wait for seqtk processes to complete ===

    for tool, sample_id, proc, output_file in processes:
        stdout, stderr = proc.communicate()
        if output_file is not None:        
            output_file.close()
        if proc.returncode == 0:
            print(f"[{tool}] {sample_id} completed successfully.")
        else:
            print(f"[{tool}] {sample_id} failed.")
            print(stderr.decode())

    # === Return to working directory ===
    os.chdir(work_dir)
    print(f"[Done] Returned to working directory: {work_dir}")

# === Dereplication ===

def run_dereplicate(sample_ids, config):

    """
    Runs mmseqs2 dereplication and awk processing for multiple samples in parallel using subprocess.Popen.
    """
    processes = []

    for sample_id in sample_ids:
        # === Define paths from config ===
        out_checkv = config["OUT_CHECKV_GENOMAD"]
        out_derep = config["OUT_DEREP"]
        subset_spades = f"{out_checkv}/{sample_id}/subset_spades.fasta"
        cluster_dir = f"{out_derep}/{sample_id}"
        cluster_res = f"{cluster_dir}/clusterRes"
        tmp_dir = f"{cluster_dir}/tmp"
        input_fasta = f"{cluster_res}_all_seqs.fasta"
        cleaned_fasta = f"{cluster_dir}/cleaned_clusterRes_all_seqs.fasta"

        # === Create output directory ===
        os.makedirs(cluster_dir, exist_ok=True)

        # === Run mmseqs easy-cluster for dereplication ===       
        print(f"[mmseqs] Clustering {subset_spades} for sample {sample_id}")
        cmd_mmseqs = [
            "conda", "run", "-n", "mmseqs2_env",
            "mmseqs", "easy-cluster", subset_spades, cluster_res, tmp_dir,
            "--min-seq-id", "0.99", "-c", "0.90", "--cov-mode", "1"
        ]
        proc_mmseqs = subprocess.Popen(cmd_mmseqs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        processes.append(("mmseqs", sample_id, proc_mmseqs))

        # Wait for mmseqs to finish before proceeding
        stdout, stderr = proc_mmseqs.communicate()
        if proc_mmseqs.returncode == 0:
            print(f"[mmseqs] {sample_id} completed successfully.")
        else:
            print(f"[mmseqs] {sample_id} failed.")
            print(stderr.decode())

        print(f"[awk] Removing duplicate headers for {input_fasta} in sample {sample_id}")
        awk_cmd = (
            r"""awk '/^>/{if($0!=prev){print; prev=$0}} !/^>/' """
            + input_fasta +
            f" > {cleaned_fasta}"
        )
        proc_awk = subprocess.Popen(awk_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        processes.append(("awk", sample_id, proc_awk))

        # Wait for awk to finish before proceeding
        stdout, stderr = proc_awk.communicate()
        if proc_awk.returncode == 0:
            print(f"[awk] {sample_id} completed successfully.")
        else:
            print(f"[awk] {sample_id} failed.")
            print(stderr.decode())

        print("[Done] Dereplication completed for all samples.")

# === Clustering ===

def run_cluster_all(config):

    """
    Clusters all dereplicated sequences using mmseqs2 in a conda environment.
    """
    # === Define paths from config ===
    out_derep = config["OUT_DEREP"]
    out_cluster = config["OUT_CLUSTER"]
    work_dir = config["WORK_DIR"]
    derep_fasta = f"{out_derep}/dereplicated.fasta"
    cluster_res = f"{out_cluster}/clusterRes"
    tmp_dir = f"{out_cluster}/tmp"
    rep_seq_src = f"{out_cluster}/clusterRes_rep_seq.fasta"
    rep_seq_dst = f"{work_dir}/query"

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

# === Make BLAST DB ===

def make_blast_db(config):

    """
    Scans DB_DIR for all .fasta files and builds BLAST databases using apptainer and makeblastdb.
    """
    # === Define paths from config ===
    db_dir = config["DB_DIR"]
    max_db_size = config["MAX_DB_SIZE"]
    db_list_path = f"{db_dir}/db-list"

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

def run_launch_blast(config):

    """
    Splits FASTA files, runs BLAST per split file, and merges results using Python functions.
    """
    # === Define paths from config ===
    prog = "05B_launchblast"
    work_dir = config["WORK_DIR"]
    fasta_dir = config["FASTA_DIR"]
    split_size = config["FA_SPLIT_FILE_SIZE"]
    results_dir = f"{work_dir}/results_testing/{prog}"
    files_list_path = f"{fasta_dir}/fasta-files"
    query_dir = f"{work_dir}/query"

    # === Reinitialize results and query directory ===
    if os.path.exists(results_dir):
        shutil.rmtree(results_dir)
    os.makedirs(results_dir)
    
    # === List all .fasta files in FASTA_DIR ===
    with open(files_list_path, "w") as file_list:
        for root, dirs, files in os.walk(fasta_dir):
            for file in files:
                if file.endswith(".fasta"):
                    full_path = os.path.join(root, file)
                    file_list.write(os.path.relpath(full_path, fasta_dir) + "\n")

    with open(files_list_path) as f:
        fasta_files = [line.strip() for line in f.readlines()]

    if not fasta_files:
        print("No FASTA files found.")
        return

    print(f"Found {len(fasta_files)} FASTA files.")

    for i, file_rel in enumerate(fasta_files, 1):
        file_path = os.path.join(fasta_dir, file_rel)
        file_name = os.path.basename(file_path)
        print(f"\n[{i}] Processing {file_name}")

        out_dir = os.path.join(results_dir, file_name)
        split_dir = os.path.join(query_dir, "fa_split")
        config["SPLIT_DIR"] = split_dir

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

        split_files_list_path = os.path.join(split_dir, "split_files_list.txt")
        config["SPLIT_FILES_LIST"] = split_files_list_path
        with open(split_files_list_path, "w") as f:
            for sfile in split_files:
                f.write(sfile + "\n")

        # === Run BLAST on each split file ===
        for task_id in range(num_split):
            print(f"Running BLAST for split {task_id}")
            run_blast(config, task_id)

        # === Merge results to GFF ===
        config["FILE_NAME"] = file_name
        merge_blast(config)

    print("\n All BLAST jobs completed.")

# === Blast ===

def run_blast(config, slurm_array_task_id):

    """
    Runs the BLAST command for each split file against each database.
    """

    # === Define paths from config ===
    work_dir = config["WORK_DIR"]
    db_dir = config["DB_DIR"]
    split_files_list = config["SPLIT_FILES_LIST"]
    split_dir = config["SPLIT_DIR"]
    results_dir = f"{work_dir}/results_testing/05C_blast"
    blast_type = config["BLAST_TYPE"]
    eval_param = config["EVAL"]
    out_fmt = config["OUT_FMT"]
    max_target_seqs = config["MAX_TARGET_SEQS"]

    # Read the SPLIT_FILES_LIST to get the split file names
    with open(split_files_list, 'r') as f:
        split_files = [line.strip() for line in f.readlines()]

    # Select the split file based on the SLURM_ARRAY_TASK_ID
    split_file = split_files[slurm_array_task_id]

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

# === Merge Blast ===

def merge_blast(config):

    """
    Merges BLAST results and converts them to GFF format.
    """
    # === Define paths from config ===
    work_dir = config["WORK_DIR"]
    db_dir = config["DB_DIR"]
    results_dir = f"{work_dir}/results_testing/05D_mergeblast"
    db_list_path = f"{db_dir}/db-list"
    file_name = config["FILE_NAME"] 

    # Create the results directory (remove prior runs and create new)
    os.makedirs(results_dir, exist_ok=True)

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

# === Annotation ===

def annotate_blast(config):

    """
    Loops over all annotation CSV files and runs the annotation script using subprocess.run.
    Replicates shell logic: looping over files, constructing arguments, and running script.
    """
    # === Define paths from config === 
    hits_file = config['BLAST_HITS']
    annotations_dir = config['ANNOTATIONS']
    output_dir = config['OUTPUT']
    script_path = f"{config['SCRIPT_DIR']}/solution1_manual.py"
    pctid = config['PCTID']
    length = config['LENGTH']

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
    config = make_config("./config_py.sh")

    sample_ids_file = f"{config['XFILE_DIR']}/{config['XFILE']}"
    
    # Read sample IDs from the file
    with open(sample_ids_file, "r") as f:
        sample_ids = [line.strip() for line in f if line.strip()]
    
    unzip_fasta(sample_ids, config)

    # === Run each pipeline step in the correct order ===

    run_genomad(sample_ids, config)
    run_checkv_genomad(sample_ids, config)
    run_dereplicate(sample_ids, config)
    run_cluster_all(config)
    make_blast_db(config)
    run_launch_blast(config)
    annotate_blast(config)

if __name__ == "__main__":
    main()

