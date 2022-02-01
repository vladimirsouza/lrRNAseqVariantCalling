### download the reference GRCh38.p13 for all chromossomes
### and create a single fasta file with all chromossomes


setwd("/home/vbarbo/project_2021/paper_analysis/reference/genome")

system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/")
ind <- readLines("index.html")
ind <- sub(".+href=\"", "", ind)
ind <- sub("\">chr.+$", "", ind)
ind <- ind [ grep("\\.fna.gz$", ind) ]

cmd <- gettextf("wget %s", ind)
lapply(cmd, system)

files <- sub(".+/", "", ind)

cmd <- gettextf("gunzip %s", files)
lapply(cmd, system)

files <- sub("\\.gz$", "", files)
cmd <- gettextf( "cat %s > GRCh38.p13_all_chr.fasta",
                  paste(files, collapse=" ") )
system(cmd)

