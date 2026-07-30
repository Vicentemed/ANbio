# ============================================================
# EBI JDispatcher REST API — MUSCLE, Clustal Omega, MAFFT
# API publica, sin cuenta requerida
# ============================================================

EBI_BASE <- "https://www.ebi.ac.uk/Tools/services/rest"

# ---- Enviar job de alineamiento multiple ----
ebi_submit <- function(tool, sequences_fasta, email = "lesp.ags@salud.gob.mx",
                        extra_params = list()) {
  url <- paste0(EBI_BASE, "/", tool, "/run")

  params <- c(
    list(email = email, sequence = sequences_fasta),
    extra_params
  )

  tryCatch({
    resp <- httr2::request(url) |>
      httr2::req_method("POST") |>
      httr2::req_body_form(!!!params) |>
      httr2::req_timeout(30) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    status <- httr2::resp_status(resp)
    body   <- trimws(httr2::resp_body_string(resp))

    if (status == 200) {
      list(success = TRUE, job_id = body, tool = tool)
    } else {
      list(success = FALSE, error = paste0("HTTP ", status, ": ", substr(body, 1, 300)))
    }
  }, error = function(e) {
    list(success = FALSE, error = paste("Error de envio EBI:", conditionMessage(e)))
  })
}

# ---- Estado del job ----
ebi_status <- function(tool, job_id) {
  url <- paste0(EBI_BASE, "/", tool, "/status/", job_id)

  tryCatch({
    resp <- httr2::request(url) |>
      httr2::req_timeout(10) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    if (httr2::resp_status(resp) == 200) {
      status_txt <- trimws(httr2::resp_body_string(resp))
      # Posibles valores: RUNNING, FINISHED, ERROR, FAILURE, NOT_FOUND, QUEUED
      list(success = TRUE, status = status_txt)
    } else {
      list(success = FALSE, status = "ERROR")
    }
  }, error = function(e) {
    list(success = FALSE, status = "ERROR", error = conditionMessage(e))
  })
}

# ---- Descargar resultado ----
ebi_result <- function(tool, job_id, result_type = "aln-fasta") {
  url <- paste0(EBI_BASE, "/", tool, "/result/", job_id, "/", result_type)

  tryCatch({
    resp <- httr2::request(url) |>
      httr2::req_timeout(30) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    if (httr2::resp_status(resp) == 200) {
      list(success = TRUE, content = httr2::resp_body_string(resp))
    } else {
      list(success = FALSE, error = paste("No se pudo descargar resultado:", httr2::resp_status(resp)))
    }
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

# ---- Listar tipos de resultado disponibles ----
ebi_result_types <- function(tool, job_id) {
  url <- paste0(EBI_BASE, "/", tool, "/resulttypes/", job_id)

  tryCatch({
    resp <- httr2::request(url) |>
      httr2::req_timeout(10) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    if (httr2::resp_status(resp) == 200) {
      parsed <- httr2::resp_body_json(resp)
      types  <- sapply(parsed$types, function(t) t$identifier)
      list(success = TRUE, types = types)
    } else {
      list(success = TRUE, types = c("aln-fasta", "aln-clustalw", "phylotree"))
    }
  }, error = function(e) {
    list(success = TRUE, types = c("aln-fasta", "aln-clustalw", "phylotree"))
  })
}

# ============================================================
# Funciones por herramienta especifica
# ============================================================

# MUSCLE — alineamiento multiple de secuencias
ebi_muscle <- function(sequences_fasta, email = "lesp.ags@salud.gob.mx") {
  ebi_submit("muscle", sequences_fasta, email,
             extra_params = list(format = "fasta", order = "aligned"))
}

# Clustal Omega — alineamiento multiple rapido
ebi_clustalo <- function(sequences_fasta, email = "lesp.ags@salud.gob.mx",
                          stype = "dna") {
  ebi_submit("clustalo", sequences_fasta, email,
             extra_params = list(stype = stype, outfmt = "fa"))
}

# MAFFT — alineamiento multiple preciso
ebi_mafft <- function(sequences_fasta, email = "lesp.ags@salud.gob.mx") {
  ebi_submit("mafft", sequences_fasta, email,
             extra_params = list(format = "fasta", order = "aligned"))
}

# Guardar resultado de alineamiento a archivo local
ebi_save_result <- function(tool, job_id, dest_dir, result_type = "aln-fasta") {
  res <- ebi_result(tool, job_id, result_type)
  if (!res$success) return(res)

  ext <- switch(result_type,
    "aln-fasta"    = ".fasta",
    "aln-clustalw" = ".aln",
    "phylotree"    = ".nwk",
    ".txt"
  )
  out_file <- file.path(dest_dir, paste0(tool, "_", job_id, ext))
  writeLines(res$content, out_file)
  list(success = TRUE, file = out_file, content = res$content)
}
