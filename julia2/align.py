"""
Run alignments on the small indexes.
"""

from glob import glob
import logging
import os

logger = logging.getLogger("julia2.align")

import julia2.job_tracking as job_tracking
import julia2.utils as utils


def _load_sequence_ids(sequence_name_list):
    with open(sequence_name_list) as fh:
        return [os.path.basename(line.split(".fasta")[0]).strip() for line in fh if line.strip()]


def _canonical_index_sample_id(index_id):
    index_sample_id = index_id.split("_")[0]
    if "c74742" in index_sample_id:
        return "s020"
    if "c49446" in index_sample_id:
        return "s018"
    return index_sample_id


def _sample_lane(sample_id, project_config):
    numeric_id = int(sample_id.split("s")[1])
    if numeric_id <= project_config.num_samples_per_lane - 1:
        return 1
    return 2


def _reads_sample_ids(project_config):
    return [f"s{str(i).zfill(3)}" for i in range(1, project_config.num_samples + 1)]


def _should_align(index_id, reads_sample_id, mode, project_config):
    index_sample_id = _canonical_index_sample_id(index_id)

    if mode == "all":
        return True

    if mode == "true_auto":
        return f"{reads_sample_id}_" in index_id

    if mode == "allo":
        return project_config.sample_to_taxon_short[index_sample_id] != project_config.sample_to_taxon_short[reads_sample_id]

    if mode == "taxon_auto":
        return project_config.sample_to_taxon_short[index_sample_id] == project_config.sample_to_taxon_short[reads_sample_id]

    reads_lane = _sample_lane(reads_sample_id, project_config)
    index_lane = _sample_lane(index_sample_id, project_config)

    if mode == "same_lane":
        return index_lane == reads_lane

    if mode == "other_lane":
        return index_lane != reads_lane

    raise ValueError(f"Unsupported alignment mode: {mode}")


def _build_alignment_plan(sequence_name_list, mode, project_config):
    plan = []
    for index_id in _load_sequence_ids(sequence_name_list):
        for reads_sample_id in _reads_sample_ids(project_config):
            if _should_align(index_id, reads_sample_id, mode, project_config):
                plan.append((index_id, reads_sample_id))
    return plan


def _alignment_key(index_id, reads_sample_id):
    return ("alignment", index_id, reads_sample_id)


def _latest_alignment_rows(project_config):
    rows = job_tracking.load_job_records(project_config)
    latest = {}
    for row in rows:
        if row["job_type"] != "alignment":
            continue
        key = (row["job_type"], row["index_id"], row["reads_sample_id"])
        latest[key] = row
    return latest


def _filter_alignment_plan(plan, project_config, resume=False):
    latest = _latest_alignment_rows(project_config)
    to_submit = []
    skipped_completed = 0
    skipped_active = 0
    resubmitting_failed = 0

    for index_id, reads_sample_id in plan:
        row = latest.get(_alignment_key(index_id, reads_sample_id))
        if not resume or row is None:
            to_submit.append((index_id, reads_sample_id))
            continue

        state = row["state"] or "UNKNOWN"
        if state == "COMPLETED":
            skipped_completed += 1
            continue
        if state in job_tracking.ACTIVE_STATES or state == "SUBMITTED":
            skipped_active += 1
            continue
        if state in job_tracking.FAILED_STATES:
            resubmitting_failed += 1
            to_submit.append((index_id, reads_sample_id))
            continue

        to_submit.append((index_id, reads_sample_id))

    return {
        "to_submit": to_submit,
        "skipped_completed": skipped_completed,
        "skipped_active": skipped_active,
        "resubmitting_failed": resubmitting_failed,
    }


def _reset_alignment_records(plan, project_config):
    plan_keys = {_alignment_key(index_id, reads_sample_id) for index_id, reads_sample_id in plan}
    removed = job_tracking.remove_job_records(
        project_config,
        lambda row: (row["job_type"], row["index_id"], row["reads_sample_id"]) in plan_keys)
    return removed


def run_alignment(reads_sample_id, index_id, system_config, project_config):
    lane = _sample_lane(reads_sample_id, project_config)
    reads_sample_num = reads_sample_id.split("s")[1]
    logger.debug("Running alignment for reads %s and index %s", reads_sample_num, index_id)
    sbatch_template, cpus = utils.create_sbatch_template(system_config.slurm_settings,
                                                         project_config,
                                                         cpus=True,
                                                         align_index="ALIGN")
    dir_1_filename = glob(
        f"{project_config.project_dir}/raw_reads/lane{lane}-{reads_sample_id}*R1*")[0]

    sbatch_cmds = f"""
echo "index_{index_id} read_{reads_sample_id}"

bowtie2 -f --threads {cpus} -x {project_config.project_dir}/indexes/{index_id}_index/{index_id}_index -U {dir_1_filename} > {project_config.project_dir}/output/alignment_database_data/raw/index_{index_id}_read_{reads_sample_id}.sam
"""

    sbatch_text = f"""{sbatch_template}

{sbatch_cmds}
"""
    utils.run_slurm_job(sbatch_text,
                        f"align_index_{index_id}_reads_{reads_sample_num}",
                        project_config,
                        system_config,
                        job_type="alignment",
                        index_id=index_id,
                        reads_sample_id=reads_sample_id)


def _submit_alignment_plan(system_config,
                           project_config,
                           sequence_name_list,
                           mode,
                           resume=False,
                           reset=False):
    plan = _build_alignment_plan(sequence_name_list, mode, project_config)

    if reset:
        removed = _reset_alignment_records(plan, project_config)
        logger.info("Reset %s tracked alignment records for %s target jobs", removed, len(plan))

    summary = _filter_alignment_plan(plan, project_config, resume=resume)
    to_submit = summary["to_submit"]

    logger.info("Alignment plan for mode=%s target_jobs=%s submit=%s skipped_completed=%s skipped_active=%s resubmitting_failed=%s",
                mode,
                len(plan),
                len(to_submit),
                summary["skipped_completed"],
                summary["skipped_active"],
                summary["resubmitting_failed"])

    for index_id, reads_sample_id in to_submit:
        run_alignment(reads_sample_id, index_id, system_config, project_config)


def run_all_samples(system_config, project_config, sequence_name_list, resume=False, reset=False):
    _submit_alignment_plan(system_config, project_config, sequence_name_list, "all", resume=resume, reset=reset)


def run_all_true_auto_samples(system_config, project_config, sequence_name_list, resume=False, reset=False):
    _submit_alignment_plan(system_config, project_config, sequence_name_list, "true_auto", resume=resume, reset=reset)


def run_all_allo_samples(system_config, project_config, sequence_name_list, resume=False, reset=False):
    _submit_alignment_plan(system_config, project_config, sequence_name_list, "allo", resume=resume, reset=reset)


def run_all_taxon_auto_samples(system_config, project_config, sequence_name_list, resume=False, reset=False):
    _submit_alignment_plan(system_config, project_config, sequence_name_list, "taxon_auto", resume=resume, reset=reset)


def run_all_intra_lane_samples(system_config, project_config, sequence_name_list, resume=False, reset=False):
    _submit_alignment_plan(system_config, project_config, sequence_name_list, "same_lane", resume=resume, reset=reset)


def run_all_cross_lane_samples(system_config, project_config, sequence_name_list, resume=False, reset=False):
    _submit_alignment_plan(system_config, project_config, sequence_name_list, "other_lane", resume=resume, reset=reset)
