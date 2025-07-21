import os
import subprocess
import shutil
from academy.manager import Manager
import asyncio
from concurrent.futures import ThreadPoolExecutor
from academy.manager import Manager
from academy.exchange.local import LocalExchangeFactory




# === Turn config file into a dictionary of variables ===

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

def read_sample_ids(sample_ids_file):
    # Read sample IDs from the file
    with open(sample_ids_file, "r") as f:
        sample_ids = [line.strip() for line in f if line.strip()]
    return sample_ids


def unzip_fasta(spades_gz, unzipped_spades_path):

    """
    Check and unzip the FASTA files for a list of sample IDs using subprocess.
    """
    import subprocess
    import os

    # Skip if unzipped file already exists
    if os.path.exists(unzipped_spades_path):
        print(f"[INFO] File already exists: {unzipped_spades_path}")
        unzipped_spades = unzipped_spades_path
        return unzipped_spades

    # Attempt to unzip the file
    try:
        subprocess.run(["gzip", "-d", spades_gz], check=True)
        print(f"[INFO] Successfully unzipped: {spades_gz}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to unzip {spades_gz}: {e}")

    unzipped_spades = unzipped_spades_path
    return unzipped_spades

# === Subprocess Wrapper ===

import asyncio
import subprocess
'''
async def run_subprocess(cmd, **kwargs):
    return await asyncio.to_thread(subprocess.run, cmd, check=True, **kwargs)
'''                       
async def run_subprocess(cmd, **kwargs):
    try:
        return await asyncio.to_thread(subprocess.run, cmd, check=True, **kwargs)
    except subprocess.CalledProcessError as e:
        print("Command failed:", e.cmd)
        print("Return code:", e.returncode)
        print("Output:", e.output)
        print("Stderr:", e.stderr)
       

# === GeNomad Agent ===

from academy.agent import Agent, action
import os
import subprocess

class GeNomadAgent(Agent):
    @action
    async def run_genomad(self, unzipped_spades: str, genomad_output_dir: str, db: str) -> str:
        if not unzipped_spades or not os.path.exists(unzipped_spades):
            raise ValueError(f"Invalid input file: {unzipped_spades}")
        if not genomad_output_dir:
            raise ValueError(f"Invalid database path: {db}")

        cmd = [
            "conda", "run", "-n", "genomad_env",
            "genomad", "end-to-end", "--cleanup", "--restart", 
            unzipped_spades, genomad_output_dir, db
        ]
        await run_subprocess(cmd)
        genomad_virus =  os.path.join(genomad_output_dir, "contigs_summary", "contigs_virus.fna")
        return genomad_virus

# === CheckV Agent ===

import subprocess
import os

class CheckVAgent(Agent):
    @action
    async def run_checkv(
        self, checkv_parser: str, parse_length: str, work_dir: str,
        unzipped_spades: str, genomad_virus: str, checkv_output_dir: str,
        parse_input: str, selection_csv: str, checkvdb: str
    ) -> str:
        os.makedirs(checkv_output_dir, exist_ok=True)

        cmd_checkv = [
            "conda", "run", "-n", "checkv_env", "checkv", "end_to_end",
            genomad_virus, checkv_output_dir, "-t", "4", "-d", checkvdb
        ]
        cmd_parser = [
            "conda", "run", "-n", "r_env", "Rscript", checkv_parser,
            "-i", parse_input, "-l", parse_length, "-o", selection_csv
        ]
        cmd_seqtk = [
            "conda", "run", "-n", "seqtk_env", "seqtk", "subseq",
            unzipped_spades, selection_csv
        ]

        await run_subprocess(cmd_checkv)
        await run_subprocess(cmd_parser)

        subset_spades = os.path.join(checkv_output_dir, "subset_spades.fasta")
        with open(subset_spades, "w") as out_f:
            await run_subprocess(cmd_seqtk, stdout=out_f)
        os.chdir(work_dir)
        return subset_spades

# === Dereplication/Clustering Agent ===

import os
import shutil
import time
from academy.agent import Agent, action

class DereplicationClusteringAgent(Agent):

    @action
    async def run_dereplicate(
        self, sample_id: str, subset_spades: str, cluster_dir: str,
        cluster_res_derep: str, tmp_dir_derep: str, input_fasta: str,
        cleaned_fasta: str, out_derep: str) -> str:
        os.makedirs(cluster_dir, exist_ok=True)

        # Run dereplication
        cmd_mmseqs_derep = [
            "conda", "run", "-n", "mmseqs2_env",
            "mmseqs", "easy-cluster", subset_spades, cluster_res_derep, tmp_dir_derep,
            "--min-seq-id", "0.99", "-c", "0.90", "--cov-mode", "1"
        ]
        await run_subprocess(cmd_mmseqs_derep)

        # Clean headers with awk
        cmd_awk = (
            r"""awk '/^>/{if($0!=prev){print; prev=$0}} !/^>/' """
            + input_fasta + f" > {cleaned_fasta}"
        )
        await run_subprocess(cmd_awk, shell=True)

        # Touch done flag
        done_flag = os.path.join(out_derep, f"done_{sample_id}.flag")
        with open(done_flag, "w") as f:
            f.write("done\n")
        derep_fasta = os.path.join(out_derep, "dereplicated.fasta")
        return derep_fasta

    @action
    async def run_cluster(
        self, sample_ids: list[str], out_derep: str, derep_fasta: str,
        out_cluster: str, cluster_res_cluster: str, tmp_dir_cluster: str,
        rep_seq_src: str, rep_seq_dst: str) -> str:
        # Wait for all dereplication flags
        done_flags = [os.path.join(out_derep, f"done_{sid}.flag") for sid in sample_ids]
        print("[wait] Waiting for all dereplication steps to complete...")
        while True:
            done_count = sum(os.path.exists(flag) for flag in done_flags)
            if done_count == len(done_flags):
                break
            time.sleep(5)
        print("[done] All dereplication samples complete.")

        # Combine all cleaned_clusterRes_all_seqs.fasta files
        with open(derep_fasta, 'w') as outfile:
            for root, _, files in os.walk(out_derep):
                for file in files:
                    if file.endswith("cleaned_clusterRes_all_seqs.fasta"):
                        with open(os.path.join(root, file)) as infile:
                            shutil.copyfileobj(infile, outfile)

        # Run clustering
        os.makedirs(out_cluster, exist_ok=True)
        cmd_mmseqs_cluster = [
            "conda", "run", "-n", "mmseqs2_env",
            "mmseqs", "easy-cluster", derep_fasta, cluster_res_cluster, tmp_dir_cluster,
            "--min-seq-id", "0.95", "-c", "0.75", "--cov-mode", "1"
        ]
        await run_subprocess(cmd_mmseqs_cluster)

        # Copy rep seqs to query dir
        os.makedirs(rep_seq_dst, exist_ok=True)
        shutil.copy(rep_seq_src, os.path.join(rep_seq_dst, "clusterRes_rep_seq.fasta"))
        query_dir = rep_seq_dst
        return query_dir


# === Make BLAST DB ===

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

    return db_name 

# === BLAST Agent ===

import os
import shutil
from academy.agent import Agent, action

class BLASTAgent(Agent):

    async def _make_blast_db(self, db_dir: str, max_db_size: str, db_list_path: str):
        """
        Check if BLAST databases already exist, and run makeblastdb if they don't.
        """
        print(f"[BLAST DB] Checking in: {db_dir}")
        os.makedirs(db_dir, exist_ok=True)
        os.chdir(db_dir)

        # Generate db-list file using os.walk
        with open(db_list_path, "w") as db_list:
            for root, _, files in os.walk("."):
                for file in files:
                    if file.endswith(".fasta"):
                        rel_path = os.path.join(root, file).lstrip("./")
                        db_list.write(rel_path + "\n")

        if not os.path.exists(db_list_path) or os.path.getsize(db_list_path) == 0:
            print(f"[BLAST DB] Cannot find or empty database list: {db_list_path}")
            return

        # Process each FASTA file
        with open(db_list_path) as f:
            for i, line in enumerate(f, start=1):
                db_file = line.strip()
                db_name = os.path.splitext(os.path.basename(db_file))[0]
                db_prefix = os.path.join(db_dir, db_name)

                # Skip if already made
                if all(os.path.exists(f"{db_prefix}.{ext}") for ext in ["nhr", "nin", "nsq"]):
                    print(f"[BLAST DB] Database already exists: {db_name}")
                    continue

                print(f"[BLAST DB] Creating: {db_name}")
                cmd = [
                    "conda", "run", "-n", "blast_env",
                    "makeblastdb",
                    "-title", db_name,
                    "-out", db_prefix,
                    "-in", db_file,
                    "-dbtype", "nucl",
                    "-max_file_sz", str(max_db_size)
                ]
                subprocess.run(cmd, check=True)
        print("[BLAST DB] Setup complete.")
        return db_name

    @action
    async def run_full_blast(
        self, work_dir: str, split_size: int,
        results_dir: str, query_dir: str, db_dir: str,
        blast_results_dir: str, blast_type: str, eval_param: float,
        out_fmt: int, max_target_seqs: int, merge_results_dir: str,
        max_db_size: int, db_list_path: str) -> str:
        
        # === Make BLAST DB if needed ===
        db_name = await self._make_blast_db(db_dir, max_db_size, db_list_path)

        # === Reinitialize BLAST results directory ===
        if os.path.exists(results_dir):
            shutil.rmtree(results_dir)
        os.makedirs(results_dir, exist_ok=True)

        # === Find input FASTA file ===
        fasta_file = None
        for root, _, files in os.walk(query_dir):
            for file in files:
                if file.endswith(".fasta"):
                    fasta_file = os.path.join(root, file)
                    break
            if fasta_file:
                break
        if not fasta_file:
            raise FileNotFoundError("No .fasta file found in query directory.")

        file_name = os.path.basename(fasta_file)

        # === Split FASTA using faSplit ===

        split_dir = os.path.join(query_dir, "fa_split")

        if os.path.exists(split_dir):
            shutil.rmtree(split_dir)
        os.makedirs(split_dir)

        print(f"[faSplit] Splitting {file_name} into chunks of size {split_size}")
        cmd_split = [
            "conda", "run", "-n", "fasplit_env",
            "faSplit", "about", fasta_file, str(split_size), f"{split_dir}/"
        ]
        await run_subprocess(cmd_split)

        # === Save list of split files ===
        split_files = sorted([f for f in os.listdir(split_dir) if f.endswith(".fa")])
        num_split = len(split_files)
        print(f"[faSplit] {num_split} split files created.")

        split_files_list = os.path.join(split_dir, "split_files_list.txt")
        with open(split_files_list, "w") as f:
            for sfile in split_files:
                f.write(sfile + "\n")

        # === Run BLAST on each split file ===
        db_list_path = os.path.join(db_dir, "db-list")
        with open(db_list_path, 'r') as f:
            databases = [line.strip() for line in f.readlines()]

        for task_id, split_file in enumerate(split_files):
            for i, db in enumerate(databases, 1):
                print(f"[BLAST] Split {task_id}, DB {i}: {db}")
                results_by_db = os.path.join(blast_results_dir, db, split_file)
                os.makedirs(results_by_db, exist_ok=True)
                blast_out = os.path.join(results_by_db, f"{split_file}.blastout")
                blast_db = os.path.join(db_dir, db)

                cmd_blast = [
                    "conda", "run", "-n", "blast_env", blast_type,
                    "-num_threads", "48",
                    "-db", blast_db,
                    "-query", os.path.join(split_dir, split_file),
                    "-out", blast_out,
                    "-evalue", str(eval_param),
                    "-outfmt", str(out_fmt),
                    "-max_target_seqs", str(max_target_seqs)
                ]
                await run_subprocess(cmd_blast)

        # === Merge BLAST results and convert to GFF ===
        os.makedirs(merge_results_dir, exist_ok=True)

        for i, db in enumerate(databases, 1):
            print(f"[Merge] Processing DB {i}: {db}")
            results_by_db = os.path.join(merge_results_dir, db)
            os.makedirs(results_by_db, exist_ok=True)

            blast_out_dir = os.path.join(work_dir, "results", "05C_blast", db, file_name)
            blast_results = os.path.join(results_by_db, f"{file_name}.txt")
            blast_gff = os.path.join(results_by_db, f"{file_name}.gff")

            with open(blast_results, 'w') as outfile:
                for result_file in os.listdir(blast_out_dir):
                    with open(os.path.join(blast_out_dir, result_file), 'r') as infile:
                        outfile.write(infile.read())

            print(f"[Convert] Writing GFF for {db}")
            with open(blast_results, 'r') as infile, open(blast_gff, 'w') as outfile:
                for line in infile:
                    fields = line.strip().split('\t')
                    if len(fields) > 7:
                        gff_line = f"{fields[0]}\tblast\tgene\t{fields[6]}\t{fields[7]}\t.\t.\t.\tID=Gene{fields[6]};Name={fields[1]}\n"
                        outfile.write(gff_line)

        print("[BLAST] All steps complete.")
        hits_file = os.path.join(merge_results_dir, "AVrC_allrepresentatives.fasta", "clusterRes_rep_seq.fasta.txt")
        return hits_file

# === Annotation ===

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

    final = "Pipeline Complete."
    return final

import asyncio
from concurrent.futures import ThreadPoolExecutor
from academy.manager import Manager
from academy.exchange.local import LocalExchangeFactory

# Helper for per-sample pipeline
async def process_sample(sample_id, config, genomad_handle, checkv_handle, cluster_handle, first_sample_id):
    # === Unzip ===
    spades_gz = os.path.join(config['SPADES_DIR'], sample_id, "contigs.fasta.gz")
    unzipped_spades_path = os.path.join(config['SPADES_DIR'], sample_id, "contigs.fasta")
    unzipped_spades = unzip_fasta(spades_gz, unzipped_spades_path)

    # === GeNomad ===
    genomad_output_dir = os.path.join(config['OUT_GENOMAD'], sample_id)
    db = config["GENOMAD_DB"]
    genomad_virus = await(await genomad_handle.run_genomad(unzipped_spades, genomad_output_dir, db))

    # === CheckV ===
    checkv_parser = config["CHECKV_PARSER"]
    parse_length = str(config["PARSE_LENGTH"])
    work_dir = config["WORK_DIR"]
    checkv_output_dir = os.path.join(config["OUT_CHECKV_GENOMAD"], sample_id)
    parse_input = os.path.join(checkv_output_dir, "contamination.tsv")
    selection_csv = os.path.join(checkv_output_dir, "selection2_viral.csv")
    checkvdb = config["CHECKVDB"]
    subset_spades = await(await checkv_handle.run_checkv(
        checkv_parser, parse_length, work_dir, unzipped_spades, genomad_virus,
        checkv_output_dir, parse_input, selection_csv, checkvdb))

    # === Dereplication ===
    cluster_dir = os.path.join(config["OUT_DEREP"], sample_id)
    cluster_res_derep = os.path.join(cluster_dir, "clusterRes")
    tmp_dir_derep = os.path.join(cluster_dir, "tmp")
    input_fasta = f"{cluster_res_derep}_all_seqs.fasta"
    cleaned_fasta = os.path.join(cluster_dir, "cleaned_clusterRes_all_seqs.fasta")
    out_derep = config["OUT_DEREP"]
    derep_fasta = await(await cluster_handle.run_dereplicate(
        sample_id, subset_spades, cluster_dir, cluster_res_derep,
        tmp_dir_derep, input_fasta, cleaned_fasta, out_derep))

    # === Only return derep_fasta from the first sample ===
    if sample_id == first_sample_id:
        return derep_fasta
    else:
        return None

# === Main function ===
async def main():
    async with await Manager.from_exchange_factory(
        factory=LocalExchangeFactory(),
        executors=ThreadPoolExecutor()
    ) as manager:
        genomad_handle = await manager.launch(GeNomadAgent())
        checkv_handle = await manager.launch(CheckVAgent())
        cluster_handle = await manager.launch(DereplicationClusteringAgent())
        blast_handle = await manager.launch(BLASTAgent())

        # === Load configuration ===
        config_path = os.path.join(os.getcwd(), "config_py.sh")
        config = make_config(config_path)
        sample_ids_file = os.path.join(config['XFILE_DIR'], config['XFILE'])
        sample_ids = read_sample_ids(sample_ids_file)

        # === Run all samples in parallel ===
        first_sample_id = sample_ids[0]  # Choose one sample to return the derep_fasta

        per_sample_tasks = [
            asyncio.create_task(process_sample(sid, config, genomad_handle, checkv_handle, cluster_handle, first_sample_id))
            for sid in sample_ids
        ]

        results = await asyncio.gather(*per_sample_tasks)

        # Only one of them should return a non-None derep_fasta
        derep_fasta = next((r for r in results if r is not None), None)

        if derep_fasta is None:
            raise ValueError("Dereplicated FASTA not found from any sample.")

        # === Cluster ===
        out_cluster = config["OUT_CLUSTER"]
        work_dir = config["WORK_DIR"]
        cluster_res_cluster = os.path.join(out_cluster, "clusterRes")
        tmp_dir_cluster = os.path.join(out_cluster, "tmp")
        rep_seq_src = os.path.join(out_cluster, "clusterRes_rep_seq.fasta")
        rep_seq_dst = os.path.join(work_dir, "query")
        out_derep = config["OUT_DEREP"]
        query_dir = await(await cluster_handle.run_cluster(
            sample_ids, out_derep, derep_fasta, out_cluster,
            cluster_res_cluster, tmp_dir_cluster, rep_seq_src, rep_seq_dst))

        # === BLAST ===
        db_dir = config["DB_DIR"]
        max_db_size = config["MAX_DB_SIZE"]
        db_list_path = os.path.join(db_dir, "db-list")
        prog = "05B_launchblast"
        fasta_dir = config["FASTA_DIR"]
        split_size = config["FA_SPLIT_FILE_SIZE"]
        results_dir = os.path.join(work_dir, "results_testing", prog)
        files_list_path = os.path.join(fasta_dir, "fasta-files")
        blast_results_dir = os.path.join(work_dir, "results_testing", "05C_blast")
        blast_type = config["BLAST_TYPE"]
        eval_param = config["EVAL"]
        out_fmt = config["OUT_FMT"]
        max_target_seqs = config["MAX_TARGET_SEQS"]
        merge_results_dir = os.path.join(work_dir, "results_testing", "05D_mergeblast")
        hits_file = await(await blast_handle.run_full_blast(
            work_dir, split_size, results_dir, query_dir, db_dir,
            blast_results_dir, blast_type, eval_param, out_fmt,
            max_target_seqs, merge_results_dir, max_db_size, db_list_path))

        # === Annotation ===
        annotations_dir = config['ANNOTATIONS']
        out_annotate = config['OUTPUT']
        script_path = os.path.join(config['SCRIPT_DIR'], "solution1_manual.py")
        pctid = config['PCTID']
        length = config['LENGTH']
        final = annotate_blast(hits_file, annotations_dir, out_annotate, script_path, pctid, length)
        print(final)

if __name__ == "__main__":
    asyncio.run(main())


