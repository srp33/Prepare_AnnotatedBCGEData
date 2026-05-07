
sort_samples = function(x) {
  return(sapply(x, sort_sample))
}

sort_sample = function(x) {
  parts = str_split(x, "______")[[1]]
  parts = sort(parts)

  return(paste(parts[1], parts[2], sep = "______"))
}

data = read_tsv("/Data/merged_doppelgang_results.tsv.gz") %>%
  mutate(sample1 = str_replace(sample1, "GSE\\d+:", "")) %>%
  mutate(sample2 = str_replace(sample2, "GSE\\d+:", "")) %>%
  filter(sample1 != sample2) %>%
  distinct() %>%
  unite(sample_merged, sample1, sample2, sep = "______") %>%
  mutate(sample_merged = sort_samples(sample_merged)) %>%
  separate(sample_merged, into = c("sample1", "sample2"), sep = "______") %>%
  distinct()

ggplot(data, aes(x = expr.similarity)) +
  geom_histogram() +
  theme_bw() +
  xlab("Similarity score") +
  ylab("Count") +
  geom_vline(xintercept = 0.94, col = "red", linetype = "dashed")

ggsave("doppelgangr_similarity_scores.pdf")

data <- filter(data, expr.similarity > 0.94) %>%
  select(sample1, sample2) %>%
  distinct() %>%
  arrange(sample1)

# group_by(data, sample1) %>%
#   summarize(count = n()) %>%
#   View()