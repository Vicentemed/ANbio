# ============================================================
# CGE/ResFinder REST API — genes de resistencia antimicrobiana
# Servicio publico de la DTU (Technical University of Denmark)
# ============================================================

CGE_RESFINDER_URL <- "https://cge.food.dtu.dk/services/ResFinder/4.6/"

# ---- Enviar secuencia a ResFinder ----
resfinder_submit <- function(fasta_content = NULL, fasta_file = NULL,
                              species = "other", acquired = TRUE,
                              point = FALSE, identity = 90, coverage = 60) {
  tryCatch({
    req <- httr2::request(CGE_RESFINDER_URL) |>
      httr2::req_method("POST") |>
      httr2::req_timeout(60) |>
      httr2::req_error(is_error = \(r) FALSE)

    if (!is.null(fasta_file)) {
      req <- req |> httr2::req_body_multipart(
        inputfile       = curl::form_file(fasta_file),
        inputfiletype   = "assembly",
        acquired        = if (acquired) "true" else "false",
        point           = if (point) "true" else "false",
        species         = species,
        threshold       = as.character(identity),
        min_cov         = as.character(coverage)
      )
    } else if (!is.null(fasta_content)) {
      req <- req |> httr2::req_body_form(
        sequencedata    = fasta_content,
        inputfiletype   = "assembly",
        acquired        = if (acquired) "true" else "false",
        point           = if (point) "true" else "false",
        species         = species,
        threshold       = as.character(identity),
        min_cov         = as.character(coverage)
      )
    } else {
      return(list(success = FALSE, error = "Debe proporcionar fasta_content o fasta_file"))
    }

    resp   <- httr2::req_perform(req)
    status <- httr2::resp_status(resp)
    body   <- httr2::resp_body_string(resp)

    if (status %in% 200:302) {
      parsed <- tryCatch(jsonlite::fromJSON(body, simplifyVector = FALSE), error = \(e) NULL)
      if (!is.null(parsed)) {
        job_id <- parsed$jobid %||% parsed$id %||% NULL
        if (!is.null(job_id))
          return(list(success = TRUE, job_id = job_id, raw = body))
      }
      # Algunos servidores devuelven el resultado directamente
      list(success = TRUE, job_id = NULL, result_direct = body)
    } else {
      list(success = FALSE, error = paste0("ResFinder HTTP ", status))
    }
  }, error = function(e) {
    list(success = FALSE, error = paste("Error ResFinder:", conditionMessage(e)))
  })
}

# ---- Estado / resultado de ResFinder ----
resfinder_status <- function(job_id) {
  url <- paste0(CGE_RESFINDER_URL, "status/", job_id)
  tryCatch({
    resp   <- httr2::request(url) |> httr2::req_timeout(10) |>
      httr2::req_error(is_error = \(r) FALSE) |> httr2::req_perform()
    if (httr2::resp_status(resp) == 200) {
      parsed <- tryCatch(httr2::resp_body_json(resp), error = \(e) NULL)
      status <- if (!is.null(parsed)) parsed$status %||% "running" else "running"
      list(success = TRUE, status = status)
    } else {
      list(success = FALSE, status = "error")
    }
  }, error = function(e) {
    list(success = FALSE, status = "error", error = conditionMessage(e))
  })
}

resfinder_result <- function(job_id) {
  url <- paste0(CGE_RESFINDER_URL, "result/", job_id)
  tryCatch({
    resp <- httr2::request(url) |> httr2::req_timeout(20) |>
      httr2::req_error(is_error = \(r) FALSE) |> httr2::req_perform()
    if (httr2::resp_status(resp) == 200) {
      body   <- httr2::resp_body_string(resp)
      parsed <- tryCatch(jsonlite::fromJSON(body, simplifyVector = FALSE), error = \(e) NULL)
      list(success = TRUE, raw = body, parsed = parsed)
    } else {
      list(success = FALSE, error = "No se pudo descargar resultado de ResFinder")
    }
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

# ---- Parsear resultado de ResFinder a tabla ----
resfinder_to_table <- function(result_parsed) {
  tryCatch({
    genes <- result_parsed$resfinder_results$results
    if (is.null(genes) || length(genes) == 0)
      return(data.frame(Gen = character(), Identidad = numeric(),
                        Cobertura = numeric(), Antibiotico = character(),
                        stringsAsFactors = FALSE))

    do.call(rbind, lapply(genes, function(g) {
      data.frame(
        Gen         = g$resistance_gene %||% g$gene_name %||% "N/A",
        Identidad   = as.numeric(g$identity %||% NA),
        Cobertura   = as.numeric(g$coverage %||% NA),
        Antibiotico = paste(g$phenotypes %||% "N/A", collapse = ", "),
        stringsAsFactors = FALSE
      )
    }))
  }, error = function(e) {
    data.frame()
  })
}

# ---- Especies validas para ResFinder ----
RESFINDER_SPECIES <- c(
  "other",
  "Escherichia coli",
  "Klebsiella pneumoniae",
  "Acinetobacter baumannii",
  "Pseudomonas aeruginosa",
  "Staphylococcus aureus",
  "Enterococcus faecium",
  "Enterococcus faecalis",
  "Salmonella",
  "Campylobacter"
)
