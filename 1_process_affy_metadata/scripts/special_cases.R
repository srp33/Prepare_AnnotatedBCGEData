
if (gseID == "GSE1561") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("title")) %>%
    dplyr::select(-starts_with("description"))
}

if (gseID == "GSE2034") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("title")) %>%
    dplyr::select(-starts_with("description"))
}

if (gseID == "GSE2603") {
  metadata <- metadata %>%
    filter(str_detect(title, "^B")) %>% # the other samples are from cell lines
    dplyr::select(-c("title", "name")) %>%
    dplyr::select(-starts_with("description"))
}

if (gseID == "GSE2990") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "ggi")) %>%
    mutate(across(where(is.character), ~replace(., . %in% c("KJ67", "KJ68", "KJ69", "KJX46", "KJX38", "KJ117"), NA)))
}

if (gseID == "GSE3744") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description")) %>%
    rename(tissue_source = characteristics)
}

if (gseID == "GSE4611") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "title"))
}

if (gseID == "GSE5327") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "lms_status")) %>%
    dplyr::rename(er_status = title) %>% 
    mutate(across(er_status, ~str_replace(., " human primary breast tumor ", ","))) %>%
    separate("er_status", c("er_status", "Patient_ID"), sep = ",")
}

if (gseID == "GSE5460") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("title")) %>%
    dplyr::select(-starts_with("description"))
}

if (gseID == "GSE5764") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description")) %>%
    separate("characteristics", c("tissue", "mastectomy", "postmenupausal", "subtype"), sep = cumsum(c(13, 25, 28, 40))) %>%
    mutate(mastectomy = str_trim(mastectomy)) %>%
    mutate(postmenupausal = str_trim(postmenupausal))
}

if (gseID == "GSE5847") {
  metadata <- metadata %>%
    mutate(description = sub(";", ",", description)) %>%
    separate("description", c("chemotherapy", "ER_status", "Her2Neu", "TNM_stage", "description"), sep = ", ") %>%
    mutate(Clinical_IBC = str_extract(description, "^Clinical.*")) %>%
    mutate(Cyclin_E = str_extract(description, "^Cyclin.*")) %>%
    mutate(across(chemotherapy, ~str_replace(., "Chemo: ", ""))) %>%
    mutate(across(ER_status, ~str_replace(., "ER: ", ""))) %>%
    mutate(across(Her2Neu, ~str_replace(., "Her2Neu: ", ""))) %>%
    mutate(across(TNM_stage, ~str_replace(., "Stage: ", ""))) %>%
    mutate(across(Clinical_IBC, ~str_replace(., "Clinical IBC:", ""))) %>%
    mutate(Clinical_IBC = str_trim(Clinical_IBC)) %>%
    mutate(across(Cyclin_E, ~str_replace(., "Cyclin E: ", ""))) %>%
    dplyr::select(-c("title", "er_status", "patient_id", "tnm_stage", "description")) 
}

if (gseID == "GSE6434") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("title")) %>%
    dplyr::select(-starts_with("description"))
}

if (gseID == "GSE7378") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "characteristics_4", "description")) %>%
    rename(er_status = characteristics)
}

if (gseID == "GSE7390") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description", "filename", "veridex_risk", "risknpi", "risksg"))
}

if (gseID == "GSE7904") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "description_1", "title")) %>%
    dplyr::rename(breast_cancer_subtype = characteristics) %>%
    mutate(across(breast_cancer_subtype, ~str_replace(., "NO", "Normal organelle"))) %>%
    mutate(across(breast_cancer_subtype, ~str_replace(., "NB", "Normal breast")))
    # mutate(across(breast_cancer_subtype, ~case_when(. == "NO" ~ "Normal organelle", . == "NB" ~ "Normal breast", TRUE ~ as.character(.))))
}

if (gseID == "GSE8977") {
  metadata <- metadata %>%
    dplyr::select(-"description") %>%
    rename(tissue_source = characteristics) %>%
    rename(tissue_type = title) %>%
    mutate(tissue_type = "Stroma")
}

if (gseID == "GSE9195") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("title")) %>%
    dplyr::select(-starts_with("description"))
}

if (gseID == "GSE9574") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("title")) %>%
    dplyr::select(-starts_with("description"))
}

if (gseID == "GSE10780") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description")) %>%
    dplyr::rename(tissue_type = characteristics)
}

if (gseID == "GSE10797") {
  metadata <- metadata %>%
    dplyr::select(-"description") %>%
    rename(tissue_source = characteristics) %>%
    rename(replicate = title) %>%
    mutate(across(replicate, ~str_replace(., "cancer_epithelial_rep", ""))) %>%
    mutate(across(replicate, ~str_replace(., "cancer_stroma_rep", "")))
}

if (gseID == "GSE10810") {
  metadata <- metadata %>%
    dplyr::select(-title) %>%
    separate("description", c("paired_status", "Patient_ID"), sep = "  ") %>%
    mutate(across(paired_status, ~str_replace(., "Control ", ""))) %>%
    mutate(across(paired_status, ~str_replace(., "Tumor ", ""))) %>%
    mutate(across(paired_status, ~str_replace(., "not paired", "No"))) %>%
    mutate(across(paired_status, ~str_replace(., "paired", "Yes")))
}

if (gseID == "GSE11121") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description", "storage"))
}

if (gseID == "GSE12093") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "characteristics", "description")) 
}

if (gseID == "GSE12276") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("title")) %>%
    dplyr::select(-starts_with("description"))
}

if (gseID == "GSE16391") {
  metadata <- metadata %>%
    rename(`tumor_size_<=_2cm_=_1_>_2cm_=2` = size) %>%
    rename(`treatment_Letrozol_=_0_Tamoxifen_=_1` = treatment) %>%
    rename(`node_node_negative_=_0_node_positive_=_1` = node) %>%
    rename(`local_therapy_BCS/RT_=_1_BCS/no_RT_=_2_Mx/RT_=_3` = local_therapy) %>%
    rename(tumor_grade = grade) %>%
    rename(`ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` = er_pgr) %>%
    rename(`post_menopausal_status_before_chemotherapy_=_1_after_chemotherapy_=_2` = post_menopausal_status) %>%
    dplyr::select(-starts_with(c("description", "post_menopausal_status"))) %>%
    dplyr::select(-c("title", "tissue", "cluster_id", "ggi"))
}

if (gseID == "GSE16446") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description", "final_analysis"))
}

if (gseID == "GSE16873") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "label_protocol_1", "description", "tissue"))
}

if (gseID == "GSE17705") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description", "tissue"))
}

if (gseID == "GSE17907") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("title")) %>%
    dplyr::select(-starts_with("description"))
}

if (gseID == "GSE18864") {
  metadata <- metadata %>%
    dplyr::select(-("description")) %>%
    separate("er_pr_her2_status", c("er_status", "pr_status", "her2_status"), sep = "/")
}

if (gseID == "GSE19615") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("title")) %>%
    dplyr::select(-starts_with("description")) %>%
    dplyr::select(-starts_with("tissue"))
}

if (gseID == "GSE19697") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description", "percent")) 
}

if (gseID == "GSE20086") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("title")) %>%
    dplyr::select(-starts_with("description"))
}

if (gseID == "GSE20181") {
  metadata <- metadata %>%
    dplyr::select(-starts_with(c("description", "tissue", "characteristics", "subject"))) %>% 
    dplyr::select(-(c("administration", "agent", "dose", "gender", "time"))) %>% 
    mutate(across(title , ~str_replace(., "female;breast tumor;", "")))  %>%
    mutate(across(title , ~str_replace(., ": ", ";"))) %>%
    mutate(gender = "female") %>%
    separate("title", c("Patient_ID", "temp", "Treatment_Response"), sep = ";") %>%
    separate("temp", c("agent", "dose", "route_of_administration", "time"), sep = ",")
}

if (gseID == "GSE20194") {
  metadata <- metadata %>%
    dplyr::select(-starts_with(c("title", "description", "tissue")))
}

if (gseID == "GSE20271") {
  metadata <- metadata %>%    
    dplyr::select(-starts_with(c("title", "description", "array", "surgery date"))) %>%
    dplyr::select(-c("dlda30_pred_1_pcr_0_rd", "dlda30_score"))
}

if (gseID == "GSE20437") {
  metadata <- metadata %>%
    dplyr::select(-starts_with(c("tissue", "description"))) %>%
    rename(histology = title) %>%
    mutate(across(histology, ~str_replace(., "reduction mammoplasty", "normal"))) %>%
    mutate(across(histology, ~str_replace(., "\\d+$", ""))) 
}

if (gseID == "GSE20685") {
  metadata <- metadata %>%
    dplyr::select(-starts_with(c("title", "tissue", "characteristics", "description", "subtype")))
}

if (gseID == "GSE20711") {
  metadata <- metadata %>%
    dplyr::select(-("title")) %>%
    dplyr::select(-("description")) %>%
    dplyr::select(-starts_with("methylation")) %>%
    dplyr::select(-starts_with("quality"))
}

if (gseID == "GSE21653") {
  metadata <- metadata %>%
    dplyr::select(-("title")) %>%
    dplyr::select(-("description")) %>%
    dplyr::select(-starts_with("tissue"))
}

if (gseID == "GSE21947") {
  metadata <- metadata %>%
    rename(histology = title) %>%
    mutate(er_status = ifelse(is.na(er_status), specimen, er_status)) %>%
    mutate(histology = str_sub(histology, 1, 39)) %>%
    mutate(er_status = str_sub(er_status, 1, 3)) %>%
    dplyr::select(-c("description", "breast_cancer_patient_id", "disease_state", "specimen", "tissue"))
}

if (gseID == "GSE22513") {
  metadata <- metadata %>%
    dplyr::select(-("title")) %>%
    dplyr::select(-("description")) %>%
    dplyr::select(-starts_with("tissue"))
}

if (gseID == "GSE24185") {
  metadata <- metadata %>%
    dplyr::select(-("title")) %>%
    dplyr::select(-starts_with("description")) %>%
    dplyr::select(-starts_with("disease"))
}

if (gseID == "GSE25055") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "sample_id", "tissue", "characteristics_23", "dlda30_prediction", "ggi_class", "set_class"))
}

if (gseID == "GSE25065") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "sample_id", "tissue", "dlda30_prediction", "ggi_class", "set_class"))
}

if (gseID == "GSE26910") {
  metadata <- metadata %>%
    dplyr::filter(grepl("breast", title)) %>%
    dplyr::select(-c("description", "tissue", "title"))
}

if (gseID == "GSE28796") {
  metadata <- metadata %>%
    dplyr::select(-starts_with(c("description", "title", "tissue")))
}

if (gseID == "GSE28821") {  #examine associated journal article. suggests TNBC
  metadata <- metadata %>%
    dplyr::select(-("title")) %>%
    dplyr::select(-starts_with("description")) %>%
    dplyr::select(-starts_with("tissue"))
}

if (gseID == "GSE31138") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "tissue")) %>%
    dplyr::select(-starts_with("description"))
}

if (gseID == "GSE31192") {
  metadata <- metadata %>%
    separate("title", c("patient_id", "A", "er_status", "cell_type_1", "cell_type_2"), sep = ", ") %>%
    rename(cell_type_3 = cell_type) %>%
    dplyr::select(-c("A", "description", "description_1", "tissue"))
}

if (gseID == "GSE31519") {
  metadata <- metadata %>%
    dplyr::select(-starts_with(c("title", "description"))) %>%
    rename(breast_cancer_subtype = tissue) %>%
    rename_with(~str_replace_all(., "_0", "")) %>%
    rename_with(~str_replace_all(., "_1", "")) %>%
    rename_with(~str_replace_all(., "2", ""))

  #biopsy type (1: surgical, 2: core needle)
  metadata <- metadata %>%
    mutate(across(biopsy_type, ~str_replace(., "surgical, 2: core needle\\)\\:", ""))) %>%
    mutate(across(biopsy_type, ~str_replace(., "1", "surgical"))) %>%
    mutate(across(biopsy_type, ~str_replace(., "2", "core needle")))

  #event (1: yes, 0: no)
  metadata <- metadata %>%
    mutate(across(event, ~str_replace(., "yes, 0: no\\)\\:", ""))) %>%
    mutate(across(event, ~str_replace(., "0", "no"))) %>%
    mutate(across(event, ~str_replace(., "1", "yes")))

  #grade (12: G1 or G2, 3: G3)
  metadata <- metadata %>%
    mutate(across(grade, ~str_replace(., "G1 or G2, 3\\: G3\\)\\:", ""))) %>%
    mutate(across(grade, ~str_replace(., "12", "G1 or G2"))) %>%
    mutate(across(grade, ~str_replace(., "3", "G3")))

  #lymph node status (0: negative, 1: positive)
  metadata <- metadata %>%
    mutate(across(lymph_node_status, ~str_replace(., "negative, 1\\: positive\\)\\: ", ""))) %>%
    mutate(across(lymph_node_status, ~str_replace(., "0", "negative"))) %>%
    mutate(across(lymph_node_status, ~str_replace(., "1", "positive")))

  #tumor size (1: up to 1 cm, 2: >1cm)
  metadata <- metadata %>%
    mutate(across(tumor_size, ~str_replace(., "up to 1 cm, 2\\: \\>1cm\\)\\:", ""))) %>%
    mutate(across(tumor_size, ~str_replace(., "1", "up to 1 cm"))) %>%
    mutate(across(tumor_size, ~str_replace(., "2", "greater than 1 cm")))
}

if (gseID == "GSE32518") {
  metadata <- metadata %>%
    rename(treatment_status = description) %>%
    mutate(across(treatment_status, ~str_replace(., "Gene expression data from breast cancer FNA biopsy, ", ""))) %>%
    mutate(across(treatment_status, ~str_replace(., "Gene expression data from breast cancer CBX biopsy, ", ""))) %>%
    dplyr::select(-c("title", "tissue", "er_by_gene"))
}

if (gseID == "GSE32646") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "tissue"))
}

if (gseID == "GSE33692") {
  metadata <- metadata %>%
    dplyr::select(-("title")) %>%
    dplyr::select(-starts_with("data_processing"))
}

if (gseID == "GSE45255") {
  metadata <- metadata %>%
    dplyr::select(-starts_with(c("title", "description"))) %>%
    rename(`endocrine_0=no_1=yes` = characteristics_9) %>%
    mutate(across(`endocrine_0=no_1=yes`, ~str_replace(., "characteristics\\: endocrine\\? \\(0\\=no, 1\\=yes\\):", "")))
}

if (gseID == "GSE46184") {
  metadata <- metadata %>%
    dplyr::select(-("title")) %>%
    dplyr::select(-("description"))
}

if (gseID == "GSE48390") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description", "tissue", "predicted_risk"))
}

if (gseID == "GSE50948") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "patid", "bgus_ct"))
}

if (gseID == "GSE57968") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "tissue_type"))
}

if (gseID == "GSE58644") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description", "paquet_et_al", "suderman_et_al"))
}

if (gseID == "GSE58984") {
  metadata <- metadata %>%
    dplyr::select(-("title")) %>%
    dplyr::select(-("description"))
}

#GSE59772 is triple negative according to GEO website
if (gseID == "GSE59772") {
  metadata <- metadata %>%
  rename(replicate = title) %>%
  dplyr::select(-("description"))
}

if (gseID == "GSE76275") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "characteristics_18", "description", "tissue", "set", "tnbc_subtype"))
}

if (gseID == "GSE81838") {
  metadata <- metadata %>%
    dplyr::select(-("title"))
}

if (gseID == "GSE86374") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description"))
}

if (gseID == "GSE90521") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description", "tumor"))
}

if (gseID == "GSE111662") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description", "organ", "subject_id", "tissue_abbrevation"))
}

if (gseID == "GSE118432") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "characteristics_4", "description"))
}

if (gseID == "GSE120129") {
  metadata <- metadata %>%
    dplyr::select(- ("title"))
}

if (gseID == "GSE167213") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description", "tissue"))
}

if (gseID == "GSE8193") {
  metadata <- metadata %>%
    dplyr::select(-"description") %>%
    rename(er_status = characteristics) %>%
    rename(pr_status = characteristics_1) %>%
    rename(Patient_ID = title)
}

if (gseID == "GSE10281") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "tissue"))%>%
    mutate(title = str_replace(title, ";.+", ""))  %>%
    mutate(title = str_replace(title, "patient", "patient "))  %>%
    separate(treatment, c("agent", "dose", "route_of_administration", "time"), sep = ",") %>%
    rename(Patient_ID = title)
}

if (gseID == "GSE11001") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description"))
}

if (gseID == "GSE12763") {
    metadata <- metadata %>%      
      dplyr::select(-c("title", "description")) %>%
      rename(subtype = characteristics)
}

if (gseID == "GSE13787") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "title"))
}

if (gseID == "GSE14017") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "characteristics")) %>%
    rename(distant_metastasis = description)
}

if (gseID == "GSE14018") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "characteristics")) %>%
    rename(distant_metastasis = description)
}

if (gseID == "GSE18728") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "tissue_type", "upin", "include_in_bl", "include_in_paired")) %>%
    rename(Patient_ID = title)
}

if (gseID == "GSE21422") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "tissue", "title"))
}

if (gseID == "GSE22093") {
  metadata <- metadata %>%
    dplyr::select(-c("array_qc", "tissue", "postchemo_tumor_size")) %>%
    dplyr::select(-starts_with(c("description","er_positive_vs_negative_by_esr1_mrna"))) %>% 
    rename(Patient_ID = title) %>%
    mutate(across(Patient_ID , ~str_replace(., "breast cancer ", ""))) %>%
    rename(er_positive_vs_negative_by_esr1_mrna_gene_expression_probe = characteristics_4) %>%
    mutate(across(er_positive_vs_negative_by_esr1_mrna_gene_expression_probe , ~str_replace(., "probe: ", "probe "))) %>%
    mutate(across(er_positive_vs_negative_by_esr1_mrna_gene_expression_probe , 
                  ~str_replace(., "er positive vs negative by esr1 mrna gene expression", ""))) %>%
    mutate(across(er_positive_vs_negative_by_esr1_mrna_gene_expression_probe , 
                   ~str_replace(., "\\(probe 205225_at\\): ", "")))
}

if (gseID == "GSE22544") {
  metadata <- metadata %>%
    dplyr::select(-c("title", "description", "tissue")) %>%
    rename(Patient_ID = patient)
}

if (gseID == "GSE23720") {
  metadata <- metadata %>%
    dplyr::select(-"title")  
}

if (gseID == "GSE25407") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "characteristics_3", "title")) %>%
    rename(Patient_ID = patient)
}

if (gseID == "GSE27562") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("description")) %>%
    dplyr::select(-c("title", "data_set"))
}

if (gseID == "GSE29431") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("characteristics_")) %>%
    dplyr::select(- "title") %>%
    mutate(disease_state = str_replace(disease_state, "none", "normal")) 
}

if (gseID == "GSE31448") {
  metadata <- metadata %>%
    dplyr::select(-starts_with("characteristics_")) %>%
    dplyr::select(-c("tissue", "description", "p53_ihc_status", "er_ihc_status", "erbb2_ht_ihc_status", "erbb2_tab_ihc_status",
    "erbb2p_ihc_status", "foxa1_ihc_status", "age_of_diagnosis_year", "egfr_ihc_status",
    "grade_sbr", "igf1r_ihc_status", "ki67_ihc_status", "mfsdel_month", "mfs", "top2a_ihc_status", "type",
    "pr_ihc_status")) %>%
    rename(er_ihc_1 = er_ihc) %>%
    rename(pr_ihc_1 = pr_ihc) %>%
    rename(Patient_ID = title) %>%
    rename(erbb2_ihc_status_1 = erbb2_ihc_status) %>%
    rename(age_at_diagnosis_1 = age_at_diagnosis) %>%
    mutate(tumor_status = ifelse(str_detect(Patient_ID, ": NB"), "Normal", "Tumor")) %>%
    mutate(age_at_diagnosis = coalesce(age_at_diagnosis_1, age_at_diagnosis_2)) %>%
    mutate(er_ihc = coalesce(er_ihc_1, er_ihc_2)) %>%
    mutate(pr_ihc = coalesce(pr_ihc_1, pr_ihc_2, prihc)) %>%
    mutate(erbb2_ihc_status = coalesce(erbb2, erbb2_ihc_status_1)) %>%
    dplyr::select(-c("er_ihc_1", "er_ihc_2", "pr_ihc_1", "pr_ihc_2", "age_at_diagnosis_1", "age_at_diagnosis_2", "erbb2",
                     "prihc","erbb2_ihc_status_1"))
}

if (gseID == "GSE41194") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "tissue")) %>%
    rename(Patient_ID = title)
}

if (gseID == "GSE41196") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "tissue")) %>%
    rename(Patient_ID = title)
}

if (gseID == "GSE41197") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "tissue")) %>%
    rename(Patient_ID = title)
}

if (gseID == "GSE42568") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "tissue")) %>%
    separate(title, c("tissue_type", "Patient_ID"), sep = ",")
}

if (gseID == "GSE43365") {
  metadata <- metadata %>%
    dplyr::select(-c("description", "ggi", "gg")) %>%
    rename(Patient_ID = title)
}

if (gseID == "GSE47109") {
  metadata <- metadata %>%
    rename(Patient_ID = title)
}

if (gseID == "GSE161533") {
  metadata <- metadata %>%
    rename(Patient_ID = title)
}

if (gseID == "GSE61304") {
  metadata <- metadata %>%
    rename(tumor_source = title) %>%
    rename(histology_group = diagnosis) %>%
    mutate(across(tumor_source, ~str_replace(., ",", ""))) %>%
    mutate(across(tumor_source, ~str_replace(., "\\d+$", ""))) 
}

if (gseID == "GSE23988") { 
  metadata <- metadata %>%
    dplyr::select(-"array_qc") %>%
    rename(Patient_ID = title) 
}

if (gseID == "GSE93332") { 
  metadata <- metadata %>%
    dplyr::select(-c("tissue", "title"))
}
