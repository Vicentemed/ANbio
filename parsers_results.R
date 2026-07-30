# ============================================================
# Parsers de resultados BV-BRC por herramienta
# ============================================================

# ---- FastqUtils / MultiQC ----
# files_content: named list  filename -> text string
# Returns: named list  sample -> list(sample, q30, gc, reads_m, dup, avg_len)
parse_fastqutils <- function(files_content) {
  result <- list()

  for (fname in names(files_content)) {
    content <- files_content[[fname]]
    if (is.null(content) || !nchar(trimws(content))) next

    # *_fastqc.html — reporte HTML de FastQC (salida real de BV-BRC FastqUtils)
    if (grepl("_fastqc\\.html$", fname, ignore.case = TRUE) ||
        (grepl("fastqc", fname, ignore.case = TRUE) &&
         grepl("FastQC Report|<html", content, ignore.case = TRUE))) {
      one <- gsub("[\r\n]", " ", content)
      grab <- function(lbl) {
        m <- regmatches(one, regexpr(paste0(lbl, "</td>[[:space:]]*<td>[^<]*"), one))
        if (length(m)) trimws(sub(".*<td>", "", m)) else NA_character_
      }
      fq_name <- grab("Filename"); if (is.na(fq_name)) fq_name <- fname
      sname <- sub("_L[0-9]+_R[12].*$", "", basename(fq_name), ignore.case = TRUE)
      sname <- sub("\\.(fastq|fq)(\\.gz)?$", "", sname, ignore.case = TRUE)
      total <- suppressWarnings(as.numeric(gsub("[^0-9]", "", grab("Total Sequences"))))
      gcv   <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", grab("%GC"))))
      slen  <- grab("Sequence length")
      avlen <- suppressWarnings(as.numeric(sub(".*[-–]", "", slen %||% "")))
      # % duplicados: FastQC reporta "Total Deduplicated Percentage" en los datos
      ddup  <- regmatches(one, regexpr("Total Deduplicated Percentage[^0-9]*([0-9.]+)", one))
      dupv  <- if (length(ddup)) {
        v <- suppressWarnings(as.numeric(regmatches(ddup, regexpr("[0-9.]+$", ddup))))
        if (!is.na(v)) round(100 - v, 1) else NA
      } else NA
      reads_m <- if (!is.na(total %||% NA)) round(total / 1e6, 2) else NA
      prev <- result[[sname]]
      if (is.null(prev)) {
        result[[sname]] <- list(sample = sname, q30 = NA, gc = gcv,
                                reads_m = reads_m, dup = dupv, avg_len = avlen)
      } else {
        # Combinar R1 + R2 de la misma muestra
        pr <- prev$reads_m %||% NA
        result[[sname]]$reads_m <- if (!is.na(pr) && !is.na(reads_m)) round(pr + reads_m, 2)
                                   else (if (!is.na(pr)) pr else reads_m)
        result[[sname]]$gc      <- if (!is.na(prev$gc %||% NA) && !is.na(gcv %||% NA))
                                     round((prev$gc + gcv) / 2, 1) else (prev$gc %||% gcv)
        result[[sname]]$avg_len <- suppressWarnings(max(prev$avg_len %||% NA, avlen,
                                                        na.rm = TRUE))
        if (is.na(result[[sname]]$dup %||% NA)) result[[sname]]$dup <- dupv
      }
      next
    }

    # multiqc_general_stats.txt — fuente principal
    if (grepl("multiqc_general_stats", fname, ignore.case = TRUE)) {
      lines <- strsplit(content, "\n")[[1]]
      lines <- lines[nchar(trimws(lines)) > 0]
      if (length(lines) < 2) next
      hdr <- strsplit(trimws(lines[1]), "\t")[[1]]

      for (row_line in lines[-1]) {
        vals <- strsplit(trimws(row_line), "\t")[[1]]
        if (!length(vals)) next
        n   <- min(length(hdr), length(vals))
        row <- setNames(as.list(vals[seq_len(n)]), hdr[seq_len(n)])

        sname    <- row$Sample %||% vals[1]
        pct_dup  <- suppressWarnings(as.numeric(
          row[["FastQC_mqc-generalstats-fastqc-percent_duplicates"]] %||%
          row[["percent_duplicates"]] %||% NA))
        pct_gc   <- suppressWarnings(as.numeric(
          row[["FastQC_mqc-generalstats-fastqc-percent_gc"]] %||%
          row[["percent_gc"]] %||% NA))
        total    <- suppressWarnings(as.numeric(
          row[["FastQC_mqc-generalstats-fastqc-total_sequences"]] %||%
          row[["total_sequences"]] %||% NA))
        avg_len  <- suppressWarnings(as.numeric(
          row[["FastQC_mqc-generalstats-fastqc-avg_sequence_length"]] %||%
          row[["avg_sequence_length"]] %||% NA))

        result[[sname]] <- list(
          sample  = sname,
          q30     = NA,
          gc      = pct_gc,
          reads_m = if (!is.na(total %||% NA)) round(total / 1e6, 2) else NA,
          dup     = pct_dup,
          avg_len = avg_len
        )
      }
    }

    # multiqc_fastqc.txt — puede contener campo pct_q30
    if (grepl("multiqc_fastqc\\.txt$", fname, ignore.case = TRUE)) {
      lines <- strsplit(content, "\n")[[1]]
      lines <- lines[nchar(trimws(lines)) > 0]
      if (length(lines) < 2) next
      hdr <- strsplit(trimws(lines[1]), "\t")[[1]]
      for (row_line in lines[-1]) {
        vals <- strsplit(trimws(row_line), "\t")[[1]]
        if (!length(vals)) next
        n   <- min(length(hdr), length(vals))
        row <- setNames(as.list(vals[seq_len(n)]), hdr[seq_len(n)])
        sname   <- row$Sample %||% vals[1]
        q30_val <- suppressWarnings(as.numeric(
          row[["pct_q30"]] %||% row[["percent_q30"]] %||% row[["q30"]] %||% NA))
        if (!is.null(result[[sname]]) && !is.na(q30_val %||% NA))
          result[[sname]]$q30 <- q30_val
      }
    }

    # fastqc_data.txt — statistics basicas de un solo sample
    if (grepl("fastqc_data\\.txt$", fname, ignore.case = TRUE)) {
      lines    <- strsplit(content, "\n")[[1]]
      in_basic <- FALSE; bstats <- list()
      for (l in lines) {
        if (grepl(">>Basic Statistics", l)) { in_basic <- TRUE; next }
        if (in_basic && grepl(">>END", l))  break
        if (in_basic && !grepl("^#", l) && nchar(trimws(l)) > 0) {
          parts <- strsplit(trimws(l), "\t")[[1]]
          if (length(parts) >= 2)
            bstats[[trimws(parts[1])]] <- trimws(parts[length(parts)])
        }
      }
      sname <- bstats$Filename %||% fname
      if (is.null(result[[sname]])) {
        reads <- suppressWarnings(
          as.numeric(gsub(",", "", bstats[["Total Sequences"]] %||% "NA")))
        result[[sname]] <- list(
          sample  = sname,
          reads_m = if (!is.na(reads %||% NA)) round(reads / 1e6, 2) else NA,
          gc      = suppressWarnings(as.numeric(bstats[["%GC"]] %||% NA)),
          q30 = NA, dup = NA, avg_len = NA
        )
      }
    }
  }
  result
}

# ---- TaxonomicClassification / Kraken2 report ----
# content: texto del archivo .report  (tab-delimitado: pct cov assign rank taxid name)
parse_taxonomy <- function(content) {
  if (is.null(content) || !nchar(trimws(content))) return(NULL)
  lines <- strsplit(content, "\n")[[1]]
  lines <- lines[nchar(trimws(lines)) > 0]

  rows <- Filter(Negate(is.null), lapply(lines, function(l) {
    parts <- trimws(strsplit(l, "\t")[[1]])
    if (length(parts) < 5) return(NULL)
    list(
      pct   = suppressWarnings(as.numeric(gsub("[%[:space:]]", "", parts[1]))),
      rank  = trimws(parts[4]),
      taxid = trimws(parts[5]),
      name  = if (length(parts) >= 6)
                trimws(paste(parts[6:length(parts)], collapse = " "))
              else trimws(parts[5])
    )
  }))

  uncl <- Filter(function(r) grepl("unclassified", r$name, ignore.case = TRUE), rows)
  spp  <- Filter(function(r) !is.na(r$pct %||% NA) && r$rank %in% c("S", "S1"), rows)
  if (!length(spp))
    spp <- Filter(function(r) !is.na(r$pct %||% NA) && r$rank == "G", rows)
  spp <- spp[order(sapply(spp, `[[`, "pct"), decreasing = TRUE)]

  # Género dominante: Kraken2 deja la mayoría de las lecturas a nivel de género
  # cuando la especie no es discriminable (p.ej. complejo A. calcoaceticus/
  # baumannii). Reportar sólo el % de especie subestima la identificación.
  gen  <- Filter(function(r) !is.na(r$pct %||% NA) && identical(r$rank, "G"), rows)
  gen  <- gen[order(sapply(gen, `[[`, "pct"), decreasing = TRUE)]
  pct_gen <- if (length(gen)) round(gen[[1]]$pct, 2) else NA

  # Dominancia relativa de la especie principal frente a la segunda:
  # mide cuán inequívoca es la asignación dentro del género.
  dom <- if (length(spp) >= 2 && (spp[[2]]$pct %||% 0) > 0)
           round(spp[[1]]$pct / spp[[2]]$pct, 1)
         else if (length(spp)) Inf else NA

  list(
    top_especie   = if (length(spp)) trimws(spp[[1]]$name) else "Desconocido",
    confianza     = if (length(spp)) round(spp[[1]]$pct, 2) else NA,
    top_genero    = if (length(gen)) trimws(gen[[1]]$name) else NA,
    pct_genero    = pct_gen,
    dominancia    = dom,
    contaminacion = if (length(uncl)) round(uncl[[1]]$pct, 2) else 0,
    top5          = head(spp, 5)
  )
}

# ---- Ensamblado a partir del contigs.fasta (salida real de BV-BRC) ----
# BV-BRC GenomeAssembly2 (unicycler) NO genera report.tsv de QUAST; entrega
# <output>_contigs.fasta con encabezados del tipo:
#   >assembly_X_contig_1 length 320186 coverage 22.6 normalized_cov 1.01
# De ahí se calculan todas las métricas de la guía LESP.
parse_assembly_fasta <- function(content) {
  if (is.null(content) || !nchar(trimws(content))) return(NULL)
  lines <- strsplit(content, "\n")[[1]]
  hdr_i <- grep("^>", lines)
  if (!length(hdr_i)) return(NULL)
  hdrs  <- lines[hdr_i]

  # Longitud declarada en el encabezado; si falta, contar bases de la secuencia
  lens <- suppressWarnings(as.numeric(
    sub(".*\\blength[ _]?([0-9]+).*", "\\1", hdrs)))
  if (any(is.na(lens))) {
    bounds <- c(hdr_i, length(lines) + 1L)
    lens <- vapply(seq_along(hdr_i), function(k) {
      seg <- lines[(bounds[k] + 1L):(bounds[k + 1L] - 1L)]
      sum(nchar(gsub("[^A-Za-z]", "", seg)))
    }, numeric(1))
  }
  covs <- suppressWarnings(as.numeric(
    sub(".*\\bcoverage[ _]?([0-9.]+).*", "\\1", hdrs)))
  covs <- covs[!is.na(covs)]

  seqs  <- paste(lines[-hdr_i], collapse = "")
  bases <- gsub("[^ACGTNacgtn]", "", seqs)
  n_gc  <- nchar(gsub("[^GCgc]", "", bases))
  n_n   <- nchar(gsub("[^Nn]",   "", bases))
  total <- sum(lens, na.rm = TRUE)

  # N50: contig en el que se acumula el 50% del ensamblado
  srt <- sort(lens, decreasing = TRUE)
  n50 <- srt[which(cumsum(srt) >= total / 2)[1]]

  list(
    contigs     = length(lens),
    n50_kb      = if (!is.na(n50)) round(n50 / 1000, 1) else NA,
    total_mb    = round(total / 1e6, 2),
    gc          = if (nchar(bases) > 0) round(100 * n_gc / nchar(bases), 1) else NA,
    largest     = as.integer(max(lens, na.rm = TRUE)),
    profundidad = if (length(covs)) round(stats::median(covs), 1) else NA,
    n_ambiguas  = n_n
  )
}

# ---- Assembly2 / QUAST report ----
# content: texto del report.tsv o report.txt de QUAST
parse_assembly <- function(content) {
  if (is.null(content) || !nchar(trimws(content))) return(NULL)
  # Si nos pasaron el FASTA de contigs, calcular métricas directamente
  if (grepl("^\\s*>", substr(content, 1, 200))) return(parse_assembly_fasta(content))
  lines <- strsplit(content, "\n")[[1]]

  stats <- list()
  for (l in lines) {
    parts <- strsplit(trimws(l), "\t")[[1]]
    if (length(parts) >= 2) {
      key        <- trimws(parts[1])
      stats[[key]] <- trimws(parts[length(parts)])
    }
  }

  gn <- function(...) {
    for (k in c(...)) {
      v <- suppressWarnings(as.numeric(gsub(",", "", stats[[k]] %||% "")))
      if (!is.na(v)) return(v)
    }
    NA_real_
  }

  n50   <- gn("N50", "N50 (>= 0 bp)")
  total <- gn("Total length", "Total length (>= 0 bp)")

  list(
    contigs  = as.integer(gn("# contigs", "# contigs (>= 0 bp)")),
    n50_kb   = if (!is.na(n50))   round(n50   / 1000, 1) else NA,
    total_mb = if (!is.na(total)) round(total / 1e6,  2) else NA,
    gc       = gn("GC (%)"),
    largest  = as.integer(gn("Largest contig"))
  )
}

# ---- Annotation RAST ----
# files_content: named list  filename -> content string
parse_annotation <- function(files_content) {
  result <- list(cds = NA, rrna = NA, trna = NA, funcion_pct = NA, hypo = NA,
                 completitud = NA, contaminacion = NA)

  for (fname in names(files_content)) {
    content <- files_content[[fname]]
    if (is.null(content) || !nchar(trimws(content))) next

    # ---- *.features.txt (salida real de BV-BRC GenomeAnnotation) ----
    # Columnas: feature_id | location | type | function | aliases | protein_md5
    if (grepl("\\.features\\.txt$", fname, ignore.case = TRUE)) {
      ln <- strsplit(content, "\n")[[1]]
      ln <- ln[nchar(trimws(ln)) > 0]
      if (length(ln) < 2) next
      parts <- strsplit(ln[-1], "\t")
      tipo  <- vapply(parts, function(p) if (length(p) >= 3) trimws(p[3]) else "", "")
      func  <- vapply(parts, function(p) if (length(p) >= 4) trimws(p[4]) else "", "")
      result$cds  <- sum(toupper(tipo) == "CDS")
      # BV-BRC marca todos los RNA con type="rna"; el tipo concreto (tRNA/rRNA)
      # va en la columna 'function' (p.ej. "tRNA-Trp-CCA", "16S rRNA").
      es_rna <- tolower(tipo) %in% c("rna", "rrna", "trna", "ncrna")
      result$rrna <- sum(es_rna & grepl("rRNA|ribosomal RNA", func, ignore.case = TRUE))
      result$trna <- sum(es_rna & grepl("tRNA", func, ignore.case = TRUE))
      cds_f <- func[tipo == "CDS"]
      result$hypo <- sum(grepl("hypothetical", cds_f, ignore.case = TRUE))
      if (result$cds > 0)
        result$funcion_pct <- round(100 * (1 - result$hypo / result$cds), 1)
      next
    }

    # ---- genome_quality_details.txt ----
    if (grepl("genome_quality_details", fname, ignore.case = TRUE)) {
      gl <- function(lbl) {
        m <- regmatches(content, regexpr(paste0(lbl, ":[ \t]*[0-9.]+"), content))
        if (length(m)) suppressWarnings(as.numeric(sub(".*:[ \t]*", "", m))) else NA_real_
      }
      result$completitud   <- gl("Completeness")
      result$contaminacion <- gl("Contamination")
      next
    }

    # genome JSON
    if (grepl("\\.json$", fname, ignore.case = TRUE)) {
      parsed <- tryCatch(
        jsonlite::fromJSON(content, simplifyVector = TRUE),
        error = \(e) NULL)
      if (!is.null(parsed)) {
        result$cds  <- parsed$cds_count          %||% result$cds
        result$rrna <- parsed$rrna_count         %||% result$rrna
        result$trna <- parsed$trna_count         %||% result$trna
        result$hypo <- parsed$hypothetical_count %||% result$hypo
      }
      next
    }

    # texto plano (genome_stats, feature_counts, summary…)
    lines <- strsplit(content, "\n")[[1]]
    for (l in lines) {
      parts <- trimws(strsplit(l, "[\t:]")[[1]])
      parts <- parts[nchar(parts) > 0]
      if (length(parts) < 2) next
      key <- tolower(parts[1])
      val <- suppressWarnings(as.numeric(gsub(",", "", parts[2])))
      if (is.na(val)) next
      if (grepl("^cds|^coding",       key)) result$cds  <- val
      if (grepl("rrna|ribosomal rna", key)) result$rrna <- val
      if (grepl("trna|transfer rna",  key)) result$trna <- val
      if (grepl("hypothetical",       key)) result$hypo <- val
    }
  }

  if (!is.na(result$cds %||% NA) && !is.null(result$hypo) &&
      !is.na(result$hypo %||% NA) && result$cds > 0)
    result$funcion_pct <- round((1 - result$hypo / result$cds) * 100, 1)

  result
}

# ---- amr-sir.txt — fenotipo AMR predicho por BV-BRC (GenomeAnnotation) ----
# Columnas: Antibiotic | Prediction | SR | F1 Avg | F1 CI Low | F1 CI High
# Devuelve data.frame  Antibiotico | Prediccion (S/I/R) | Confianza_F1
parse_amr_sir <- function(content) {
  if (is.null(content) || !nchar(trimws(content))) return(NULL)
  ln <- strsplit(content, "\n")[[1]]
  ln <- ln[nchar(trimws(ln)) > 0]
  if (length(ln) < 2) return(NULL)
  hdr <- trimws(strsplit(ln[1], "\t")[[1]])
  i_ab <- which(tolower(hdr) == "antibiotic")
  i_sr <- which(toupper(hdr) == "SR")
  i_f1 <- which(grepl("^F1 Avg$", hdr, ignore.case = TRUE))
  if (!length(i_ab) || !length(i_sr)) return(NULL)
  rows <- Filter(Negate(is.null), lapply(ln[-1], function(l) {
    p <- strsplit(l, "\t")[[1]]
    if (length(p) < max(i_ab, i_sr)) return(NULL)
    ab <- trimws(p[i_ab])
    if (!nchar(ab)) return(NULL)
    data.frame(
      Antibiotico  = gsub("_", "/", ab),
      Prediccion   = toupper(trimws(p[i_sr])),
      Confianza_F1 = if (length(i_f1) && length(p) >= i_f1)
                       round(suppressWarnings(as.numeric(p[i_f1])), 3) else NA_real_,
      stringsAsFactors = FALSE)
  }))
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# ---- MetagenomicReadMapping / CARD AMR ----
# content: texto del archivo TSV de resultados AMR
# Clase de antibiótico a partir del nombre del gen AMR (nomenclatura CARD)
amr_gene_class <- function(gen) {
  g <- tolower(gen %||% "")
  # CARD incluye entradas descriptivas ("... conferring resistance to X" o
  # mutaciones en gyrA/parC/rpoB); clasificar por la clase mencionada.
  if (grepl("fluoroquinolone|quinolone", g))                       return("Fluoroquinolonas")
  if (grepl("macrolide", g))                                       return("Macrólidos")
  if (grepl("aminoglycoside", g))                                  return("Aminoglucósidos")
  if (grepl("beta-lactam|carbapenem|cephalosporin|penicillin", g))  return("Betalactámicos")
  if (grepl("tetracycline", g))                                    return("Tetraciclinas")
  if (grepl("colistin|polymyxin", g))                              return("Polimixinas")
  if (grepl("rifampic|rifamycin", g))                              return("Rifamicinas")
  if (grepl("vancomycin|glycopeptide", g))                         return("Glicopéptidos")
  if (grepl("phenicol", g))                                        return("Fenicoles")
  if (grepl("sulfonamide|sulphonamide", g))                        return("Sulfonamidas")
  if (grepl("trimethoprim", g))                                    return("Trimetoprim")
  if (grepl("\\bgyr[ab]\\b|\\bpar[ce]\\b", g))                     return("Fluoroquinolonas")
  if (grepl("^(aph|aac|ant|aad|arm[a-z]?|rmt|str[ab])", g))        return("Aminoglucósidos")
  if (grepl("^(bla|oxa|ndm|kpc|imp|vim|ctx|tem|shv|ampc|cmy|pdc|adc)", g))
                                                                    return("Betalactámicos")
  if (grepl("^tet", g))                                            return("Tetraciclinas")
  if (grepl("^(cat|cml|flo|cmr)", g))                              return("Fenicoles")
  if (grepl("^sul", g))                                            return("Sulfonamidas")
  if (grepl("^dfr", g))                                            return("Trimetoprim")
  if (grepl("^(erm|mph|msr|mef|ere)", g))                          return("Macrólidos")
  if (grepl("^(qnr|gyr|par|qep|aac\\(6.\\)-ib-cr)", g))            return("Fluoroquinolonas")
  if (grepl("^van", g))                                            return("Glicopéptidos")
  if (grepl("^(mcr|lpx|pmr|arn|eptа|ept)", g))                     return("Polimixinas")
  if (grepl("^fos", g))                                            return("Fosfomicina")
  if (grepl("^(rpo|rif|arr)", g))                                  return("Rifamicinas")
  if (grepl("^(lnu|lsa|cfr|opt|poxt)", g))                         return("Oxazolidinonas/Lincosamidas")
  if (grepl("^(ade|mex|acr|mtr|oqx|smev|tolc|mdt|emr|nor|abe)", g)) return("Bomba de expulsión (multidrogas)")
  "Otros"
}

# ---- KMA (.res) — salida real de BV-BRC MetagenomicReadMapping con CARD ----
# Formato: #Template<TAB>Score<TAB>Expected<TAB>Template_length<TAB>
#          Template_Identity<TAB>Template_Coverage<TAB>...
# Template viene como  "CARD|<accession> <gen> [<organismo>]"
parse_kma_res <- function(content) {
  if (is.null(content) || !nchar(trimws(content))) return(NULL)
  lines <- strsplit(content, "\n")[[1]]
  lines <- lines[nchar(trimws(lines)) > 0]
  if (length(lines) < 2) return(NULL)
  hdr   <- strsplit(sub("^#", "", lines[1]), "\t")[[1]]
  hdr   <- trimws(hdr)
  i_id  <- which(hdr == "Template_Identity")
  i_cov <- which(hdr == "Template_Coverage")
  i_dep <- which(hdr == "Depth")

  rows <- Filter(Negate(is.null), lapply(lines[-1], function(l) {
    p <- strsplit(l, "\t")[[1]]
    if (length(p) < 3) return(NULL)
    tpl <- trimws(p[1])
    tpl <- sub("^CARD\\|", "", tpl)
    tpl <- sub("^\\S+\\s+", "", tpl)          # quitar accession
    org <- regmatches(tpl, regexpr("\\[[^]]+\\]", tpl))
    gen <- trimws(sub("\\s*\\[[^]]*\\]\\s*$", "", tpl))
    if (!nchar(gen)) return(NULL)
    num <- function(i) if (length(i) && length(p) >= i)
      suppressWarnings(as.numeric(trimws(p[i]))) else NA_real_
    data.frame(
      Gen         = gen,
      Identidad   = num(i_id),
      Cobertura   = num(i_cov),
      Profundidad = num(i_dep),
      Antibiotico = amr_gene_class(gen),
      Organismo   = if (length(org)) gsub("^\\[|\\]$", "", org) else "",
      stringsAsFactors = FALSE)
  }))
  if (!length(rows)) return(NULL)
  df <- do.call(rbind, rows)
  df[order(-(df$Identidad %||% 0)), , drop = FALSE]
}

parse_resistome <- function(content) {
  if (is.null(content) || !nchar(trimws(content))) return(NULL)
  # Autodetectar formato KMA (salida real de BV-BRC MetagenomicReadMapping)
  if (grepl("^#Template\\b", content) || grepl("Template_Identity", content))
    return(parse_kma_res(content))
  lines <- strsplit(content, "\n")[[1]]
  lines <- lines[nchar(trimws(lines)) > 0]
  if (length(lines) < 2) return(NULL)

  hdr  <- strsplit(trimws(lines[1]), "\t")[[1]]
  rows <- Filter(Negate(is.null), lapply(lines[-1], function(l) {
    vals <- strsplit(trimws(l), "\t")[[1]]
    if (!length(vals)) return(NULL)
    n <- min(length(hdr), length(vals))
    setNames(as.list(vals[seq_len(n)]), hdr[seq_len(n)])
  }))
  if (!length(rows)) return(NULL)

  gf <- function(r, ...) {
    for (k in c(...)) if (!is.null(r[[k]]) && nchar(r[[k]] %||% "") > 0) return(r[[k]])
    ""
  }

  df <- data.frame(
    Gen         = sapply(rows, gf, "ARO_term","gene","Gene","gene_name","AMR_gene","Name"),
    Identidad   = suppressWarnings(as.numeric(sapply(rows, gf,
                    "identity","percent_identity","Identity","% Identity"))),
    Cobertura   = suppressWarnings(as.numeric(sapply(rows, gf,
                    "coverage","percent_coverage","Coverage","% Coverage"))),
    Antibiotico = sapply(rows, gf,
                    "drug_class","Drug_class","AMR_Family","antibiotic","Drug Class"),
    stringsAsFactors = FALSE
  )
  df[nchar(df$Gen) > 0, ]
}

# ---- GenomeComparison — pangenoma (paso 9) ----
# genome_comparison.txt: fila 1 = nombres de genoma, fila 2 = cabecera real.
# Cada fila es un gen del genoma de referencia; las columnas
# comp_genome_N_hit indican si hay hit bidireccional en el genoma N.
parse_genome_comparison <- function(content) {
  if (is.null(content) || !nchar(trimws(content))) return(NULL)
  ln <- strsplit(content, "\n")[[1]]
  ln <- ln[nchar(trimws(ln)) > 0]
  if (length(ln) < 3) return(NULL)
  hdr  <- trimws(strsplit(ln[2], "\t")[[1]])
  hits <- grep("^comp_genome_[0-9]+_hit$", hdr)
  if (!length(hits)) return(NULL)
  datos <- ln[-c(1, 2)]
  p <- strsplit(datos, "\t")
  pres <- vapply(p, function(x)
    sum(vapply(hits, function(i)
      length(x) >= i && nchar(trimws(x[i])) > 0, logical(1))), integer(1))
  n_comp <- length(hits)
  list(
    total      = length(datos),                 # genes del genoma de referencia
    core       = sum(pres == n_comp),           # presentes en todos
    accesorios = sum(pres > 0 & pres < n_comp), # en algunos
    unicos     = sum(pres == 0),                # solo en la referencia
    n_genomas  = n_comp + 1L
  )
}

# ---- MLST (Galaxy mlst tool o BV-BRC CoreGenomeMLST) ----
# content: texto TSV  Filename\tScheme\tST\tallele1\t...
parse_mlst_tsv <- function(content) {
  if (is.null(content) || !nchar(trimws(content))) return(NULL)
  lines <- strsplit(content, "\n")[[1]]
  lines <- lines[nchar(trimws(lines)) > 0]
  # Descartar cabecera si empieza con "File" o "#"
  data_lines <- lines[!grepl("^#|^File|^Sample|^Filename", lines, ignore.case = TRUE)]
  if (!length(data_lines)) data_lines <- lines[-1]

  rows <- Filter(Negate(is.null), lapply(data_lines, function(l) {
    p <- strsplit(trimws(l), "\t")[[1]]
    if (length(p) < 3) return(NULL)
    data.frame(
      Muestra = sub("_contigs$", "",
                  sub("\\.(fa|fasta|fna|fastq)(\\.gz)?$", "",
                      basename(trimws(p[1])), ignore.case = TRUE)),
      ST      = trimws(p[3]),
      Alelos  = if (length(p) > 3) paste(trimws(p[4:length(p)]), collapse = "-") else "",
      Esquema = trimws(p[2]),
      stringsAsFactors = FALSE
    )
  }))
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# ---- Árbol filogenético Newick (IQ-TREE, BV-BRC PhylogeneticTree) ----
# Devuelve la cadena Newick directamente para mostrar / guardar en input
parse_newick <- function(content) {
  if (is.null(content) || !nchar(trimws(content))) return(NULL)
  lines <- strsplit(content, "\n")[[1]]
  # Ignorar líneas de comentario
  nwk <- paste(lines[!grepl("^#|^\\[", lines)], collapse = "")
  nwk <- trimws(nwk)
  if (nchar(nwk) == 0) return(NULL)
  nwk
}
