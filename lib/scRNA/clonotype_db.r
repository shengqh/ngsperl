### TCR database search ###

#load dbs
mcpas <- read.csv(parFile2)
tbadb <- read.csv(parFile3)
vdjdb <- read.csv(parFile4)

#load clonotypes file
samplename <- outFile
clon <- read.csv(parFile1)

##############################################################
#split TCR sequences into alpha and beta columns
cdr_aa <- strsplit(as.character(clon$cdr3s_aa), ";")

alpha1 <- data.frame(clon[, 1])
alpha2 <- data.frame(clon[, 1])
beta1 <- data.frame(clon[, 1])
beta2 <- data.frame(clon[, 1])

for (i in 1:nrow(clon)) {
  cdr <- cdr_aa[[i]]
  if (length(grep("TRA", cdr)) == 1) {
    alpha1[i, 2] <- cdr[1]
  } else if (length(grep("TRA", cdr)) == 2) {
    alpha1[i, 2] <- cdr[1]
    alpha2[i, 2] <- cdr[2]
  }
  
  if (length(grep("TRB", cdr)) == 1) {
    beta1[i, 2] <- cdr[grep("TRB", cdr)[1]]
  } else if (length(grep("TRB", cdr)) == 2) {
    beta1[i, 2] <- cdr[grep("TRB", cdr)[1]]
    beta2[i, 2] <- cdr[grep("TRB", cdr)[2]]
  }
}

clon2 <-
  data.frame(
    clon[, 1:4],
    cdr3s_aa_TRA1 = gsub("TRA:", "", as.character(alpha1[, 2])),
    cdr3s_aa_TRA2 = if (ncol(alpha2) == 2) {
      gsub("TRA:", "", as.character(alpha2[, 2]))
    } else {
      rep("NA", nrow(clon))
    },
    cdr3s_aa_TRB1 = gsub("TRB:", "", as.character(beta1[, 2])),
    cdr3s_aa_TRB2 = if (ncol(beta2) == 2) {
      gsub("TRB:", "", as.character(beta2[, 2]))
    } else {
      rep("NA", nrow(clon))
    },
    stringsAsFactors = F
  )
clon2[is.na(clon2)] <- "0"

#merge with databases
#mcpas
mcpas_hsa <- mcpas[which(mcpas$Species == "Human"), ]
mcpas_hsa_slim <-
  data.frame(mcpas_hsa[, c(1, 2, 5)], stringsAsFactors = F)
mcpas_hsa[] <- lapply(mcpas_hsa, as.character)

mcpas_list <- list()
for (k in 1:4) {
  cdr3s_aa <-
    c("cdr3s_aa_TRA1",
      "cdr3s_aa_TRA2",
      "cdr3s_aa_TRB1",
      "cdr3s_aa_TRB2")
  if (k == 1 | k == 2) {
    clon2_mcpas1 <-
      merge(
        clon2,
        mcpas_hsa[, c(1, 5)],
        by.x = cdr3s_aa[k],
        by.y = "CDR3.alpha.aa",
        sort = F
      )
    colnames(clon2_mcpas1)[9] <- "Pathology_McPAS-TCR"
    clon2_mcpas1_unique <-
      clon2_mcpas1[!duplicated(clon2_mcpas1[c(1, 9)]), ]
    if (nrow(clon2_mcpas1_unique) > 0) {
      clon2_mcpas1_unique["cdr3_match_McPAS-TCR"] <-
        paste0("TRA:", clon2_mcpas1_unique[, 1])
    }
  } else if (k == 3 | k == 4) {
    clon2_mcpas1 <-
      merge(
        clon2,
        mcpas_hsa[, c(2, 5)],
        by.x = cdr3s_aa[k],
        by.y = "CDR3.beta.aa",
        sort = F
      )
    colnames(clon2_mcpas1)[9] <- "Pathology_McPAS-TCR"
    clon2_mcpas1_unique <-
      clon2_mcpas1[!duplicated(clon2_mcpas1[c(1, 9)]), ]
    if (nrow(clon2_mcpas1_unique) > 0) {
      clon2_mcpas1_unique["cdr3_match_McPAS-TCR"] <-
        paste0("TRB:", clon2_mcpas1_unique[, 1])
    }
  }
  mcpas_list[[k]] <- clon2_mcpas1_unique
}

mcpas_final <- do.call("rbind", mcpas_list)
if (nrow(mcpas_final) > 0) {
  mcpas_final <- mcpas_final[, c(2:5, 9, 10)]
  mcpas_final[] <-
    sapply(mcpas_final, function(x)
      gsub("^0$", "-", x))
}
#write.csv(mcpas_final, paste0(date, "_", samplename, "_McPAS-TCR.csv"), row.names=F)

#tbadb
tbadb[] <- lapply(tbadb, as.character)

tbadb_list <- list()
for (k in 1:4) {
  cdr3s_aa <-
    c("cdr3s_aa_TRA1",
      "cdr3s_aa_TRA2",
      "cdr3s_aa_TRB1",
      "cdr3s_aa_TRB2")
  if (k == 1 | k == 2) {
    clon2_tbadb <-
      merge(clon2,
            tbadb[, c(2, 8)],
            by.x = cdr3s_aa[k],
            by.y = "CDR3.alpha.aa",
            sort = F)
    colnames(clon2_tbadb)[9] <- "Disease.name_TBAdb"
    clon2_tbadb_unique <-
      clon2_tbadb[!duplicated(clon2_tbadb[c(1, 9)]), ]
    if (nrow(clon2_tbadb_unique) > 0) {
      clon2_tbadb_unique["cdr3_match_TBAdb"] <-
        paste0("TRA:", clon2_tbadb_unique[, 1])
    }
  } else if (k == 3 | k == 4) {
    clon2_tbadb <-
      merge(clon2,
            tbadb[, c(2, 9)],
            by.x = cdr3s_aa[k],
            by.y = "CDR3.beta.aa",
            sort = F)
    colnames(clon2_tbadb)[9] <- "Disease.name_TBAdb"
    clon2_tbadb_unique <-
      clon2_tbadb[!duplicated(clon2_tbadb[c(1, 9)]), ]
    if (nrow(clon2_tbadb_unique) > 0) {
      clon2_tbadb_unique["cdr3_match_TBAdb"] <-
        paste0("TRB:", clon2_tbadb_unique[, 1])
    }
  }
  tbadb_list[[k]] <- clon2_tbadb_unique
}

tbadb_final <- do.call("rbind", tbadb_list)
if (nrow(tbadb_final) > 0) {
  tbadb_final <- tbadb_final[, c(2:5, 9, 10)]
  tbadb_final[] <-
    sapply(tbadb_final, function(x)
      gsub("^0$", "-", x))
}
#write.csv(tbadb_final, paste0(date,  "_", samplename, "_tbadb-TCR.csv"), row.names=F)

#vdjdb
vdjdb[] <- lapply(vdjdb, as.character)
vdjdb_tra <-
  vdjdb[which(vdjdb$gene == "TRA") &&
          which(vdjdb$species == "HomoSapiens"),]
vdjdb_trb <-
  vdjdb[which(vdjdb$gene == "TRB") &&
          which(vdjdb$species == "HomoSapiens"),]

vdjdb_list <- list()
for (k in 1:4) {
  cdr3s_aa <-
    c("cdr3s_aa_TRA1",
      "cdr3s_aa_TRA2",
      "cdr3s_aa_TRB1",
      "cdr3s_aa_TRB2")
  if (k == 1 | k == 2) {
    clon2_vdjdb <-
      merge(
        clon2,
        vdjdb_tra[, c(3, 12)],
        by.x = cdr3s_aa[k],
        by.y = "cdr3",
        sort = F
      )
    colnames(clon2_vdjdb)[9] <- "antigen.species_vdjdb"
    clon2_vdjdb_unique <-
      clon2_vdjdb[!duplicated(clon2_vdjdb[c(1, 9)]), ]
    if (nrow(clon2_vdjdb_unique) > 0) {
      clon2_vdjdb_unique["cdr3_match_vdjdb"] <-
        paste0("TRA:", clon2_vdjdb_unique[, 1])
    }
  } else if (k == 3 | k == 4) {
    clon2_vdjdb <-
      merge(
        clon2,
        vdjdb_trb[, c(3, 12)],
        by.x = cdr3s_aa[k],
        by.y = "cdr3",
        sort = F
      )
    colnames(clon2_vdjdb)[9] <- "antigen.species_vdjdb"
    clon2_vdjdb_unique <-
      clon2_vdjdb[!duplicated(clon2_vdjdb[c(1, 9)]), ]
    if (nrow(clon2_vdjdb_unique) > 0) {
      clon2_vdjdb_unique["cdr3_match_vdjdb"] <-
        paste0("TRB:", clon2_vdjdb_unique[, 1])
    }
  }
  vdjdb_list[[k]] <- clon2_vdjdb_unique
}

vdjdb_final <- do.call("rbind", vdjdb_list)
if (nrow(vdjdb_final) > 0) {
  vdjdb_final <- vdjdb_final[, c(2:5, 9, 10)]
  vdjdb_final[] <-
    sapply(vdjdb_final, function(x)
      gsub("^0$", "-", x))
}
#write.csv(vdjdb_final, paste0(date,  "_", samplename, "_vdjdb-TCR.csv"), row.names=F)

# put results together in one single file
final <-
  rbind(mcpas_final[, 1:4], tbadb_final[, 1:4], vdjdb_final[, 1:4])
final <- final[!duplicated(final$clonotype_id), ]
if (nrow(mcpas_final) > 0) {
  final <-
    merge(final,
          mcpas_final[, c(1, 5, 6)],
          by = "clonotype_id",
          all.x = T,
          sort = F)
}
if (nrow(tbadb_final) > 0) {
  final <-
    merge(final,
          tbadb_final[, c(1, 5, 6)],
          by = "clonotype_id",
          all.x = T,
          sort = F)
}
if (nrow(vdjdb_final) > 0) {
  final <-
    merge(final,
          vdjdb_final[, c(1, 5, 6)],
          by = "clonotype_id",
          all.x = T,
          sort = F)
}

#merge with other sample and cell columns
final <- merge(final, clon[,c(1,6:10)], by="clonotype_id")

final <- final[order(final$frequency, decreasing = T), ]
final[is.na(final)] <- "-"
write.csv(
  final,
  paste0(samplename, ".clonotype_db.csv"),
  row.names = F,
  quote = F
)
