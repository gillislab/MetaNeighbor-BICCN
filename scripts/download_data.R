
# C = Chromium 10X, SS = SmartSeq
data_url = list(
    scCv2="http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v2_AIBS/",
    scCv3="http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/10X_cells_v3_AIBS/",
    snCv2="http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v2_AIBS/",
    snCv3Z="http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_AIBS/",
    snCv3M="http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10X/processed/analysis/10X_nuclei_v3_Broad/",
    scSS="http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/SMARTer/processed/analysis/SMARTer_cells_MOp/",
    snSS="http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/SMARTer/processed/analysis/SMARTer_nuclei_MOp/"
)

common_files = list(
    ca_file="cluster.annotation.csv",
    cm_file="cluster.membership.csv",
    qc_file="QC.csv"
)

C_files = list(
    data_file="umi_counts.h5",
    metadata_file="sample_metadata.csv"
)

SS_files = list(
    data_file="exon.counts.csv.gz",
    metadata_file="sample_metadata.csv.gz"
)

joint_url="http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10X/processed/analysis/RNASeq_integrated/"

# Download all data
for (dataset_name in names(data_url)) {
    url = data_url[dataset_name]
    dir.create(dataset_name, showWarnings=FALSE)
    data_files = if (grepl("SS", dataset_name)) {SS_files} else {C_files}
    for (filename in c(data_files, common_files)) {
        download.file(file.path(data_url[dataset_name], filename), file.path(dataset_name, filename))
    }
}

# Download joint annotation
dir_name = "joint_annotation"
dir.create(dir_name, showWarnings=FALSE)
for (filename in common_files[c("ca_file", "cm_file")]) {
    download.file(file.path(joint_url, filename), file.path(dir_name, filename))
}