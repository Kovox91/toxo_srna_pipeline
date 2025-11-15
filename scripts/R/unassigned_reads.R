# Simple script to summarise the unasigned reads from the mito .bam files
setwd("/media/molgenuser/DATA/helicase_toxo_srna_pipeline/scripts/R")
files <- list.files("../../out/counts",
                    pattern = "*mito_featureCounts.txt.summary",
                    full.names = TRUE)

countSummary <- read.delim(files[1], header = T) 
for (file in files[2:length(files)]){
  counts <- read.delim(file, header = T)
  countSummary <- left_join(countSummary, counts, by = "Status")
}

coldata <- read.csv("../../samples.csv")

countSum <- pivot_longer(countSummary,
                         2:ncol(countSummary),
                         names_to = "sample") %>% 
  filter(value != 0) %>% 
  mutate(sample = base::gsub("out.filtered.", "", sample),
         sample = base::gsub("_mito_filtered.bam", "", sample)) %>% 
  pivot_wider(names_from = "Status",
              values_from = "value") %>% 
  mutate(sample = as.integer(sample),
         rel = Unassigned_NoFeatures / Assigned) %>% 
  left_join(coldata, by = "sample") %>% 
  select(sample_name, con, Assigned, Unassigned_NoFeatures, rel)
