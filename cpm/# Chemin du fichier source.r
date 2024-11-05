# Chemin du fichier source
fichier_source <- "/media/audrey/EXTERNAL_USB/differential-expression_workflow-master/00_FASQC/DKVD9.R2.fastq"

# Chemin du fichier de destination
fichier_destination <- "/home/audrey/Documents/data_fastq/DKVD9_R2.fastq"

# Lire les 800 premières lignes du fichier source
donnees <- read.csv(fichier_source, nrows = 798)

# Écrire les données dans le fichier de destination
write.csv(donnees, fichier_destination, row.names = FALSE)
