# ============================================================
# Galaxy REST API — upload, ejecucion de herramientas, resultados
# Requiere API key del usuario (generada en usegalaxy.org)
# ============================================================

GALAXY_BASE <- "https://usegalaxy.org"

galaxy_req <- function(endpoint, api_key, method = "GET", body = NULL, multipart = NULL,
                       timeout = 60) {
  url  <- paste0(GALAXY_BASE, endpoint)
  req  <- httr2::request(url) |>
    httr2::req_method(method) |>
    httr2::req_headers("x-api-key" = api_key) |>
    httr2::req_timeout(timeout) |>
    httr2::req_error(is_error = \(r) FALSE)

  if (!is.null(body))      req <- req |> httr2::req_body_json(body)
  if (!is.null(multipart)) req <- req |> httr2::req_body_multipart(!!!multipart)

  httr2::req_perform(req)
}

# ---- Validar API key y obtener info del usuario ----
galaxy_whoami <- function(api_key) {
  tryCatch({
    resp <- galaxy_req("/api/whoami", api_key)
    if (httr2::resp_status(resp) == 200) {
      info <- httr2::resp_body_json(resp)
      list(success = TRUE, username = info$username, email = info$email)
    } else {
      list(success = FALSE, error = "API key invalida o sin conexion a Galaxy")
    }
  }, error = function(e) {
    list(success = FALSE, error = paste("Error Galaxy:", conditionMessage(e)))
  })
}

# ---- Obtener o crear historial activo ----
galaxy_get_or_create_history <- function(api_key, name = "ANbio_LESP") {
  tryCatch({
    resp <- galaxy_req("/api/histories", api_key)
    if (httr2::resp_status(resp) == 200) {
      histories <- httr2::resp_body_json(resp)
      existing  <- Filter(function(h) h$name == name, histories)
      if (length(existing) > 0) {
        return(list(success = TRUE, history_id = existing[[1]]$id))
      }
    }
    # Crear historial nuevo
    resp2 <- galaxy_req("/api/histories", api_key, method = "POST",
                        body = list(name = name))
    if (httr2::resp_status(resp2) == 200) {
      hist <- httr2::resp_body_json(resp2)
      list(success = TRUE, history_id = hist$id)
    } else {
      list(success = FALSE, error = "No se pudo crear historial en Galaxy")
    }
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

# ---- Subir archivo a Galaxy ----
galaxy_upload_file <- function(api_key, history_id, local_path,
                                file_type = "fastqsanger.gz") {
  filename <- basename(local_path)

  tryCatch({
    resp <- galaxy_req(
      endpoint  = "/api/tools",
      api_key   = api_key,
      method    = "POST",
      timeout   = 600,
      multipart = list(
        tool_id    = "upload1",
        history_id = history_id,
        `inputs`   = jsonlite::toJSON(list(
          `files_0|NAME`         = filename,
          `files_0|type`         = "upload_dataset",
          `dbkey`                = "?",
          `file_type`            = file_type
        ), auto_unbox = TRUE),
        `files_0|file_data` = curl::form_file(local_path)
      )
    )

    status <- httr2::resp_status(resp)
    body   <- httr2::resp_body_string(resp)

    if (status == 200) {
      # simplifyVector=FALSE: con TRUE, outputs se convierte en data.frame y
      # parsed$outputs[[1]]$id falla ("$ operator is invalid for atomic vectors")
      parsed     <- jsonlite::fromJSON(body, simplifyVector = FALSE)
      dataset_id <- parsed$outputs[[1]]$id %||% NULL
      list(success = TRUE, dataset_id = dataset_id, filename = filename)
    } else {
      list(success = FALSE, error = paste0("Upload Galaxy fallido (", status, "): ", substr(body, 1, 200)))
    }
  }, error = function(e) {
    list(success = FALSE, error = paste("Error upload Galaxy:", conditionMessage(e)))
  })
}

# ---- Ejecutar herramienta ----
galaxy_run_tool <- function(api_key, history_id, tool_id, inputs_list) {
  body <- list(
    tool_id    = tool_id,
    history_id = history_id,
    inputs     = inputs_list
  )

  tryCatch({
    resp <- galaxy_req("/api/tools", api_key, method = "POST", body = body)
    status <- httr2::resp_status(resp)
    raw    <- httr2::resp_body_string(resp)

    if (status == 200) {
      parsed  <- jsonlite::fromJSON(raw, simplifyVector = FALSE)
      job_ids <- sapply(parsed$jobs, function(j) j$id)
      out_ids <- sapply(parsed$outputs, function(o) o$id)
      list(success = TRUE, job_ids = job_ids, output_ids = out_ids)
    } else {
      list(success = FALSE, error = paste0("HTTP ", status, ": ", substr(raw, 1, 300)))
    }
  }, error = function(e) {
    list(success = FALSE, error = paste("Error run Galaxy:", conditionMessage(e)))
  })
}

# ---- Esperar a que un dataset esté listo (state == "ok") ----
# Después de upload1, el dataset pasa por queued → running → ok.
# galaxy_run_tool lanza ToolInputsNotReadyException (HTTP 400) si se llama antes.
galaxy_wait_dataset <- function(api_key, dataset_id, history_id,
                                 max_wait = 300, poll_secs = 5) {
  deadline <- Sys.time() + max_wait
  repeat {
    state <- tryCatch({
      resp <- galaxy_req(
        paste0("/api/histories/", history_id, "/contents/", dataset_id),
        api_key, timeout = 30)
      if (httr2::resp_status(resp) == 200)
        httr2::resp_body_json(resp)$state %||% "unknown"
      else "unknown"
    }, error = \(e) "unknown")

    if (state == "ok")    return(list(success = TRUE,  state = state))
    if (state == "error") return(list(success = FALSE, state = state,
                                      error = "Dataset en estado error"))
    if (Sys.time() > deadline)
      return(list(success = FALSE, state = state,
                  error = paste0("Timeout esperando dataset (último estado: ", state, ")")))
    Sys.sleep(poll_secs)
  }
}

# ---- Outputs de un job completado ----
# Devuelve lista de dataset IDs producidos por un job
galaxy_job_outputs <- function(api_key, job_id) {
  tryCatch({
    resp <- galaxy_req(paste0("/api/jobs/", job_id, "/outputs"), api_key)
    if (httr2::resp_status(resp) == 200) {
      outs <- httr2::resp_body_json(resp, simplifyVector = FALSE)
      ids  <- vapply(outs, function(o) {
        o$dataset$id %||% o$id %||% ""
      }, character(1))
      nms  <- vapply(outs, function(o) {
        o$label %||% o$dataset$name %||% ""
      }, character(1))
      list(success = TRUE,
           ids   = ids[nchar(ids) > 0],
           names = nms[nchar(ids) > 0])
    } else {
      list(success = FALSE, error = paste0("HTTP ", httr2::resp_status(resp)))
    }
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

# ---- Metadatos de un dataset ----
galaxy_dataset_meta <- function(api_key, dataset_id) {
  tryCatch({
    resp <- galaxy_req(paste0("/api/datasets/", dataset_id), api_key)
    if (httr2::resp_status(resp) == 200)
      list(success = TRUE, meta = httr2::resp_body_json(resp))
    else
      list(success = FALSE)
  }, error = function(e) list(success = FALSE))
}

# ---- Estado del job ----
galaxy_job_status <- function(api_key, job_id) {
  tryCatch({
    resp   <- galaxy_req(paste0("/api/jobs/", job_id), api_key)
    parsed <- httr2::resp_body_json(resp)
    # Estados: new, queued, running, waiting, ok, error, paused, deleted
    list(success = TRUE, status = parsed$state %||% "unknown")
  }, error = function(e) {
    list(success = FALSE, status = "error", error = conditionMessage(e))
  })
}

# ---- Descargar dataset ----
galaxy_download_dataset <- function(api_key, dataset_id, dest_file) {
  tryCatch({
    resp <- galaxy_req(
      paste0("/api/datasets/", dataset_id, "/display"),
      api_key
    )
    if (httr2::resp_status(resp) == 200) {
      writeBin(httr2::resp_body_raw(resp), dest_file)
      list(success = TRUE, file = dest_file)
    } else {
      list(success = FALSE, error = "No se pudo descargar el dataset")
    }
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })
}

# ============================================================
# Constructores de inputs por herramienta Galaxy
# ============================================================

# SPAdes
galaxy_inputs_spades <- function(dataset_id_r1, dataset_id_r2) {
  list(
    `singlePaired|sPaired` = "paired",
    `singlePaired|input_f` = list(src = "hda", id = dataset_id_r1),
    `singlePaired|input_r` = list(src = "hda", id = dataset_id_r2),
    `kmers`                = "auto",
    `cov_cutoff`           = "off"
  )
}

GALAXY_TOOL_SPADES     <- "toolshed.g2.bx.psu.edu/repos/nml/spades/spades/4.2.0+galaxy0"
GALAXY_TOOL_UNICYCLER  <- "toolshed.g2.bx.psu.edu/repos/iuc/unicycler/unicycler/0.5.1+galaxy0"
GALAXY_TOOL_MLST       <- "toolshed.g2.bx.psu.edu/repos/iuc/mlst/mlst/2.22.0"
GALAXY_TOOL_IQTREE     <- "toolshed.g2.bx.psu.edu/repos/iuc/iqtree/iqtree/2.4.0+galaxy1"
GALAXY_TOOL_KRAKEN2    <- "toolshed.g2.bx.psu.edu/repos/iuc/kraken2/kraken2/2.17.1+galaxy0"

# Base de datos Kraken2 por defecto: RefSeq Standard-Full (arquea, bacteria,
# viral, plásmidos, humano, UniVec) — nombres NCBI, ideal para aislados clínicos.
GALAXY_KRAKEN2_DB <- "2022-07-06T094102Z_standard_prebuilt_standard_2022-06-07"

# Kraken2 — clasificación taxonómica (alternativa a BV-BRC TaxonomicClassification)
# Modo single-end: para identificar un aislado basta R1; evita construir
# colecciones pareadas en Galaxy. 'report_output' sale en formato Kraken2
# report (tabulado), que parse_taxonomy() ya sabe leer.
galaxy_inputs_kraken2 <- function(dataset_id_r1, db = GALAXY_KRAKEN2_DB) {
  list(
    `single_paired|single_paired_selector` = "no",
    `single_paired|input_sequences`        = list(src = "hda", id = dataset_id_r1),
    `kraken2_database`                     = db,
    `use_names`                            = "true",
    `report|create_report`                 = "true"
  )
}

# Unicycler
galaxy_inputs_unicycler <- function(dataset_id_r1, dataset_id_r2) {
  list(
    `paired_unpaired|input_type_selector` = "paired",
    `paired_unpaired|fastq_input1` = list(src = "hda", id = dataset_id_r1),
    `paired_unpaired|fastq_input2` = list(src = "hda", id = dataset_id_r2),
    mode = "normal"
  )
}

# MLST
galaxy_inputs_mlst <- function(dataset_id_fasta, scheme = "") {
  list(
    input_file  = list(src = "hda", id = dataset_id_fasta),
    mlst_scheme = scheme
  )
}

# IQ-TREE
galaxy_inputs_iqtree <- function(dataset_id_aln, model = "GTR+G") {
  list(
    s      = list(src = "hda", id = dataset_id_aln),
    m      = model,
    bb     = "1000",
    alrt   = "1000"
  )
}
