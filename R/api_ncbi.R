# ============================================================
# NCBI BLAST REST API — busqueda de similitud de secuencias
# API publica, sin cuenta (recomendado poner email)
# ============================================================

NCBI_BLAST_URL <- "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

# ---- Enviar secuencia a BLAST ----
ncbi_blast_submit <- function(sequence, program = "blastn", database = "nt",
                               email = "lesp.ags@salud.gob.mx") {
  tryCatch({
    resp <- httr2::request(NCBI_BLAST_URL) |>
      httr2::req_method("POST") |>
      httr2::req_body_form(
        CMD         = "Put",
        PROGRAM     = program,
        DATABASE    = database,
        QUERY       = sequence,
        FORMAT_TYPE = "XML2",
        HITLIST_SIZE = "50",
        EMAIL       = email,
        tool        = "ANbio_LESP_App"
      ) |>
      httr2::req_timeout(30) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    status <- httr2::resp_status(resp)
    body   <- httr2::resp_body_string(resp)

    if (status == 200) {
      # Extraer RID del HTML de respuesta
      rid_match <- regmatches(body, regexpr("RID = ([A-Z0-9]+)", body))
      rtoe_match <- regmatches(body, regexpr("RTOE = ([0-9]+)", body))

      rid  <- if (length(rid_match) > 0) sub("RID = ", "", rid_match) else NULL
      rtoe <- if (length(rtoe_match) > 0) as.integer(sub("RTOE = ", "", rtoe_match)) else 30

      if (!is.null(rid)) {
        list(success = TRUE, rid = rid, estimated_seconds = rtoe, program = program, database = database)
      } else {
        list(success = FALSE, error = "No se pudo obtener el RID de NCBI BLAST")
      }
    } else {
      list(success = FALSE, error = paste0("HTTP ", status))
    }
  }, error = function(e) {
    list(success = FALSE, error = paste("Error BLAST:", conditionMessage(e)))
  })
}

# ---- Consultar estado del job ----
ncbi_blast_status <- function(rid) {
  tryCatch({
    resp <- httr2::request(NCBI_BLAST_URL) |>
      httr2::req_url_query(
        CMD           = "Get",
        FORMAT_OBJECT = "SearchInfo",
        RID           = rid
      ) |>
      httr2::req_timeout(15) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    if (httr2::resp_status(resp) == 200) {
      body <- httr2::resp_body_string(resp)

      # Detectar estado en la respuesta HTML
      if (grepl("Status=READY", body))   return(list(success = TRUE, status = "FINISHED"))
      if (grepl("Status=WAITING", body)) return(list(success = TRUE, status = "RUNNING"))
      if (grepl("Status=UNKNOWN", body)) return(list(success = TRUE, status = "NOT_FOUND"))
      if (grepl("ThereAreHits=yes", body)) return(list(success = TRUE, status = "FINISHED"))

      list(success = TRUE, status = "RUNNING")
    } else {
      list(success = FALSE, status = "ERROR")
    }
  }, error = function(e) {
    list(success = FALSE, status = "ERROR", error = conditionMessage(e))
  })
}

# ---- Descargar resultados BLAST ----
ncbi_blast_results <- function(rid, format_type = "Text") {
  tryCatch({
    resp <- httr2::request(NCBI_BLAST_URL) |>
      httr2::req_url_query(
        CMD         = "Get",
        FORMAT_TYPE = format_type,
        RID         = rid,
        HITLIST_SIZE = "50"
      ) |>
      httr2::req_timeout(30) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    if (httr2::resp_status(resp) == 200) {
      content <- httr2::resp_body_string(resp)
      list(success = TRUE, content = content, rid = rid)
    } else {
      list(success = FALSE, error = paste("Error descargando resultados BLAST:", httr2::resp_status(resp)))
    }
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

# ---- Parsear resultado BLAST en texto a tabla ----
ncbi_blast_parse_text <- function(blast_text) {
  tryCatch({
    lines <- strsplit(blast_text, "\n")[[1]]

    hits <- data.frame(
      Descripcion = character(),
      Identidad   = numeric(),
      Cobertura   = numeric(),
      E_valor     = character(),
      Puntuacion  = numeric(),
      stringsAsFactors = FALSE
    )

    in_hits <- FALSE
    for (ln in lines) {
      if (grepl("^>", ln)) in_hits <- TRUE
      if (in_hits && grepl("Identities\\s*=", ln)) {
        id_match  <- regmatches(ln, regexpr("[0-9]+%", ln))
        evalue_ln <- grep("Expect\\s*=", lines, value = TRUE)[1]
        ev_match  <- if (!is.na(evalue_ln))
          regmatches(evalue_ln, regexpr("[0-9e.+-]+$", evalue_ln)) else "N/A"

        if (length(id_match) > 0) {
          hits <- rbind(hits, data.frame(
            Descripcion = "Ver BLAST completo",
            Identidad   = as.numeric(sub("%", "", id_match)),
            Cobertura   = NA_real_,
            E_valor     = as.character(ev_match),
            Puntuacion  = NA_real_,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    hits
  }, error = function(e) {
    data.frame()
  })
}

# ---- URL directa a resultados en web NCBI ----
ncbi_blast_url <- function(rid) {
  paste0("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID=", rid)
}
