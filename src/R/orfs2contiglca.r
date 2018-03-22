#!/usr/bin/env Rscript

# orfs2contiglca.r
#
# Uses the taxonomy of orfs to assign taxonomy to contigs.
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))

SCRIPT_VERSION = "0.1"
DEFAULT_RANKS = "domain,phylum,class,order,family,genus,species"

# Get arguments
option_list = list(
  make_option(
    c("--ranks"), default = DEFAULT_RANKS,
    help="Comma separated list of taxonomic ranks, default %default"
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  ),
  make_option(
    c("-V", "--version"), action="store_true", default=FALSE, 
    help="Print program version and exit"
  )
)
opt = parse_args(OptionParser(option_list=option_list))

# opt <- list(options = list(verbose = TRUE, ranks = 'domain,phylum,class,order,family,genus,species', version = FALSE), args = c('orfs2contiglca.00.tsv'))
if ( opt$options$version ) {
  write(SCRIPT_VERSION, stdout())
  quit('no')
}

logmsg = function(msg, llevel='INFO') {
  if ( opt$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}
logmsg(sprintf("Reading %s", opt$args[1]))

orfs <- read_tsv(opt$args[1], col_types = cols(.default = col_character()))

# Contigs with only one orf, can be assigned immediately
contigs <- orfs %>%
  semi_join(orfs %>% group_by(contig) %>% count() %>% filter(n == 1), by = 'contig') %>%
  select(-orf)

# Delete rows we dealt with above
orfs <- orfs %>% anti_join(contigs, by = 'contig')

ranks <- str_split(opt$options$ranks, ',')[[1]] #%>% map(quo)
rank_quos <- ranks %>% map(quo)
rank_syms <- ranks %>% map(rlang::sym)
lca <- function(ds, rs) {
  rank_counts <- list()
  for ( r in rs ) {
    #print(r)
    rank_counts[[r]] <- ds %>% filter(!is.na(!!r)) %>% distinct(!!r) %>% nrow() #%>% print()
  }
  return(rank_counts)
}
rc <- lca(data, rank_syms) %>% print()

data %>% distinct(contig) %>% pull(contig) %>% walk(print)

logmsg("Done")
