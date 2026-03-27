import pandas
import os

# -----------------------------------------------------------------------------
# Extract age information for all samples with WGS data, and exclude those 
# with age < 65.

# This query represents dataset "WGS_IPF" for domain "person" and was generated for All of Us Controlled Tier Dataset v8
dataset_25367808_person_sql = """
    SELECT
        person.person_id,
        person.gender_concept_id,
        p_gender_concept.concept_name as gender,
        person.birth_datetime as date_of_birth,
        person.race_concept_id,
        p_race_concept.concept_name as race,
        person.ethnicity_concept_id,
        p_ethnicity_concept.concept_name as ethnicity,
        person.sex_at_birth_concept_id,
        p_sex_at_birth_concept.concept_name as sex_at_birth,
        person.self_reported_category_concept_id,
        p_self_reported_category_concept.concept_name as self_reported_category 
    FROM
        `""" + os.environ["WORKSPACE_CDR"] + """.person` person 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_gender_concept 
            ON person.gender_concept_id = p_gender_concept.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_race_concept 
            ON person.race_concept_id = p_race_concept.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_ethnicity_concept 
            ON person.ethnicity_concept_id = p_ethnicity_concept.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_sex_at_birth_concept 
            ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_self_reported_category_concept 
            ON person.self_reported_category_concept_id = p_self_reported_category_concept.concept_id  
    WHERE
        person.PERSON_ID IN (SELECT
            distinct person_id  
        FROM
            `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` cb_search_person  
        WHERE
            cb_search_person.person_id IN (SELECT
                person_id 
            FROM
                `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` p 
            WHERE
                has_whole_genome_variant = 1 ) )"""

person_df = pandas.read_gbq(
    dataset_25367808_person_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook")

import datetime as dt

# data was extracted in 2025, so we can calculate age based on that.
cur_age = [ (2025 - dob.to_pydatetime().year) for dob in person_df['date_of_birth'] ]

person_df['age'] = cur_age

ov65_ids = set([str(x) for x in person_df[person_df['age'] >= 65]['person_id']])
print(len(ov65_ids))

# -----------------------------------------------------------------------------
# Extract ancestry information to exclude non-European samples.

# Predicted ancestry file
!gsutil -u $GOOGLE_PROJECT cp "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv" ./

# genetic ancestry
anc = pd.read_table("ancestry_preds.tsv")

# subject ids with predicted European ancestry
eur_ids = set(anc[anc['ancestry_pred'] == 'eur'].research_id)
eur_ids = set([str(x) for x in eur_ids])

# -----------------------------------------------------------------------------
# Fetch relatedness information to exclude related samples
!gsutil -u $GOOGLE_PROJECT cp "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/relatedness/relatedness.tsv" ./

rel = pd.read_table("relatedness.tsv")

# -----------------------------------------------------------------------------
# Fetch QC information to exclude flagged samples

!gsutil -u $GOOGLE_PROJECT cp "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv" ./

import pandas as pd

flag_df = pd.read_table("flagged_samples.tsv")
flagged_ids = set([ str(x) for x in set(flag_df.s) ])

# -----------------------------------------------------------------------------
# Identify who received antifibrotic medications (nintedanib or pirfenidone) to enrich for IPF cases.

# Medication concept ids for nintedanib and pirfenidone (45775206, 45775396)
# This query represents dataset "WGS_IPF" for domain "drug" and was generated for All of Us Controlled Tier Dataset v8
dataset_79846942_drug_sql = """
    SELECT
        d_exposure.person_id,
        d_exposure.drug_concept_id,
        d_standard_concept.concept_name as standard_concept_name,
        d_standard_concept.concept_code as standard_concept_code,
        d_standard_concept.vocabulary_id as standard_vocabulary,
        d_exposure.drug_exposure_start_datetime,
        d_exposure.drug_exposure_end_datetime,
        d_exposure.verbatim_end_date,
        d_exposure.drug_type_concept_id,
        d_type.concept_name as drug_type_concept_name,
        d_exposure.stop_reason,
        d_exposure.refills,
        d_exposure.quantity,
        d_exposure.days_supply,
        d_exposure.sig,
        d_exposure.route_concept_id,
        d_route.concept_name as route_concept_name,
        d_exposure.lot_number,
        d_exposure.visit_occurrence_id,
        d_visit.concept_name as visit_occurrence_concept_name,
        d_exposure.drug_source_value,
        d_exposure.drug_source_concept_id,
        d_source_concept.concept_name as source_concept_name,
        d_source_concept.concept_code as source_concept_code,
        d_source_concept.vocabulary_id as source_vocabulary,
        d_exposure.route_source_value,
        d_exposure.dose_unit_source_value 
    FROM
        ( SELECT
            * 
        FROM
            `""" + os.environ["WORKSPACE_CDR"] + """.drug_exposure` d_exposure 
        WHERE
            (
                drug_concept_id IN (SELECT
                    DISTINCT ca.descendant_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria_ancestor` ca 
                JOIN
                    (SELECT
                        DISTINCT c.concept_id       
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` c       
                    JOIN
                        (SELECT
                            CAST(cr.id as string) AS id             
                        FROM
                            `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` cr             
                        WHERE
                            concept_id IN (45775206, 45775396)             
                            AND full_text LIKE '%_rank1]%'       ) a 
                            ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                            OR c.path LIKE CONCAT('%.', a.id) 
                            OR c.path LIKE CONCAT(a.id, '.%') 
                            OR c.path = a.id) 
                    WHERE
                        is_standard = 1 
                        AND is_selectable = 1) b 
                        ON (ca.ancestor_id = b.concept_id)))  
                    AND (d_exposure.PERSON_ID IN (SELECT
                        distinct person_id  
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` cb_search_person  
                WHERE
                    cb_search_person.person_id IN (SELECT
                        person_id 
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` p 
                    WHERE
                        has_whole_genome_variant = 1 ) )
            )) d_exposure 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` d_standard_concept 
            ON d_exposure.drug_concept_id = d_standard_concept.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` d_type 
            ON d_exposure.drug_type_concept_id = d_type.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` d_route 
            ON d_exposure.route_concept_id = d_route.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.visit_occurrence` v 
            ON d_exposure.visit_occurrence_id = v.visit_occurrence_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` d_visit 
            ON v.visit_concept_id = d_visit.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` d_source_concept 
            ON d_exposure.drug_source_concept_id = d_source_concept.concept_id"""

drug_df = pandas.read_gbq(
    dataset_79846942_drug_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook")

print(set(drug_df.standard_concept_name))

print(drug_df.head())

antifib_ids = [str(x) for x in drug_df.person_id]
antifib_ids = set(antifib_ids)
len(antifib_ids)

# -----------------------------------------------------------------------------
# Identify who has pulmonary fibrosis based on diagnosis codes (J84.1)

# This query represents dataset "WGS_IPF" for domain "condition" and was generated for All of Us Controlled Tier Dataset v8
dataset_23273306_condition_sql = """
    SELECT
        c_occurrence.person_id,
        c_occurrence.condition_concept_id,
        c_standard_concept.concept_name as standard_concept_name,
        c_standard_concept.concept_code as standard_concept_code,
        c_standard_concept.vocabulary_id as standard_vocabulary,
        c_occurrence.condition_start_datetime,
        c_occurrence.condition_end_datetime,
        c_occurrence.condition_type_concept_id,
        c_type.concept_name as condition_type_concept_name,
        c_occurrence.stop_reason,
        c_occurrence.visit_occurrence_id,
        visit.concept_name as visit_occurrence_concept_name,
        c_occurrence.condition_source_value,
        c_occurrence.condition_source_concept_id,
        c_source_concept.concept_name as source_concept_name,
        c_source_concept.concept_code as source_concept_code,
        c_source_concept.vocabulary_id as source_vocabulary,
        c_occurrence.condition_status_source_value,
        c_occurrence.condition_status_concept_id,
        c_status.concept_name as condition_status_concept_name 
    FROM
        ( SELECT
            * 
        FROM
            `""" + os.environ["WORKSPACE_CDR"] + """.condition_occurrence` c_occurrence 
        WHERE
            (
                condition_concept_id IN (SELECT
                    DISTINCT c.concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` c 
                JOIN
                    (SELECT
                        CAST(cr.id as string) AS id       
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` cr       
                    WHERE
                        concept_id IN (45763750)       
                        AND full_text LIKE '%_rank1]%'      ) a 
                        ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                        OR c.path LIKE CONCAT('%.', a.id) 
                        OR c.path LIKE CONCAT(a.id, '.%') 
                        OR c.path = a.id) 
                WHERE
                    is_standard = 1 
                    AND is_selectable = 1) 
                OR  condition_source_concept_id IN (SELECT
                    DISTINCT c.concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` c 
                JOIN
                    (SELECT
                        CAST(cr.id as string) AS id       
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` cr       
                    WHERE
                        concept_id IN (35208072, 44824293, 44824294, 44825479, 44826684, 44827828, 44827829, 44832429, 44834773, 44835988, 45543275, 45548125, 45552899, 45562465, 45567274, 45581864, 45586679, 45591563, 45601135, 725357, 725358)       
                        AND full_text LIKE '%_rank1]%'      ) a 
                        ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                        OR c.path LIKE CONCAT('%.', a.id) 
                        OR c.path LIKE CONCAT(a.id, '.%') 
                        OR c.path = a.id) 
                WHERE
                    is_standard = 0 
                    AND is_selectable = 1)
            )  
            AND (
                c_occurrence.PERSON_ID IN (SELECT
                    distinct person_id  
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` cb_search_person  
                WHERE
                    cb_search_person.person_id IN (SELECT
                        person_id 
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` p 
                    WHERE
                        has_whole_genome_variant = 1 ) )
            )) c_occurrence 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` c_standard_concept 
            ON c_occurrence.condition_concept_id = c_standard_concept.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` c_type 
            ON c_occurrence.condition_type_concept_id = c_type.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.visit_occurrence` v 
            ON c_occurrence.visit_occurrence_id = v.visit_occurrence_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` visit 
            ON v.visit_concept_id = visit.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` c_source_concept 
            ON c_occurrence.condition_source_concept_id = c_source_concept.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` c_status 
            ON c_occurrence.condition_status_concept_id = c_status.concept_id"""

condition_df = pandas.read_gbq(
    dataset_23273306_condition_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook")

condition_df = condition_df[condition_df['source_concept_code'].notna()]

icd10_top = [ x[0:5] for x in condition_df.source_concept_code ]

condition_df['icd10_top'] = icd10_top

select_ctype_df = condition_df[condition_df['condition_type_concept_name'].isin({"Primary Condition", "EHR encounter diagnosis", "Discharge diagnosis", "EHR billing diagnosis", "EHR billing record", "EHR discharge record"})]
select_ctype_df = select_ctype_df[select_ctype_df['source_concept_code'].notna()]

any_pf_ids = set(select_ctype_df[select_ctype_df['icd10_top'] == 'J84.1'].person_id)
any_pf_ids = set([str(x) for x in any_pf_ids])

# -----------------------------------------------------------------------------
# Define cases and controls

# Enrich for cases by intersecting those with pulmonary fibrosis diagnosis codes 
# and those who received antifibrotic medications.
case_ids = any_pf_ids.intersection(antifib_ids)

# Removed QC-flagged samples
case_ids = case_ids.difference(flagged_ids)

# Controls
ctrl_ids = ov65_ids.difference(set([str(x) for x in condition_df.person_id]))
ctrl_ids = ctrl_ids.difference(flagged_ids)

# Remove one sample for each related pair
fixrel = rel[rel['i.s'].isin([ int(x) for x in case_ids ]) & rel['j.s'].isin([ int(x) for x in case_ids ])]
case_ids = case_ids.difference(set([str(x) for x in fixrel['i.s']]))

fixrel = rel[rel['i.s'].isin([ int(x) for x in case_ids ]) & rel['j.s'].isin([ int(x) for x in ctrl_ids ])]
ctrl_ids = ctrl_ids.difference(set([str(x) for x in fixrel['j.s']]))

fixrel = rel[rel['i.s'].isin([ int(x) for x in ctrl_ids ]) & rel['j.s'].isin([ int(x) for x in case_ids ])]
ctrl_ids = ctrl_ids.difference(set([str(x) for x in fixrel['i.s']]))

fixrel = rel[rel['i.s'].isin([ int(x) for x in ctrl_ids ]) & rel['j.s'].isin([ int(x) for x in ctrl_ids ])]
ctrl_ids = ctrl_ids.difference(set([str(x) for x in fixrel['j.s']]))

print(f"Number of cases: {len(case_ids)}")
print(f"Number of controls: {len(ctrl_ids)}")
