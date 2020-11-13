library(tidyverse)

files <- list.files(pattern = "*step100.tsv", full.names = TRUE, recursive = TRUE)

id <- files %>%
	str_remove_all("(^\\./|overlap_|\\.tsv)") %>%
	str_split("/") %>%
	map_chr(paste, collapse = "_")


valid_chr <- c(as.character(1:22), "X", "Y", "MT")
df <- files %>%
	setNames(nm = id) %>%
	map(read_tsv, col_names = FALSE, col_types = "cnncncnncnnncnncncn") %>%
	map(~ .x %>%
		count(X1, X16) %>%
		filter(n > 10)) %>%
	bind_rows(.id = "source") %>%
	mutate(X1 = factor(X1, levels = valid_chr))

df %>%
	ggplot(aes(X1, X16, size = n)) +
	geom_point() +
	scale_size_area() +
	theme_minimal() +
	facet_wrap(~ source) +
	labs(x = "Chromosome", y = "Repeat element")
ggsave("overlap_count.pdf", height = 10, width = 16)


df %>%
	ggplot(aes(X1, X16, fill = n)) +
	geom_tile() +
	theme_minimal() +
	facet_wrap(~ source) +
	labs(x = "Chromosome", y = "Repeat element") +
	scale_fill_gradient(low = "white", high = "orange")
