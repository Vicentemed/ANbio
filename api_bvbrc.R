# ============================================================
# BV-BRC API — autenticacion, upload, envio y estado de jobs
# ============================================================

BVBRC_AUTH_URL <- "https://user.patricbrc.org/authenticate"
BVBRC_APP_URL  <- "https://p3.theseed.org/services/app_service"
BVBRC_WS_URL   <- "https://p3.theseed.org/services/Workspace"

# Null-coalescing (disponible para todos los modulos)
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# ---- Autenticacion ----
bvbrc_auth <- function(username, password) {
  tryCatch({
    resp <- httr2::request(BVBRC_AUTH_URL) |>
      httr2::req_method("POST") |>
      httr2::req_body_form(username = username, password = password) |>
      httr2::req_timeout(15) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    if (httr2::resp_status(resp) == 200) {
      token    <- trimws(httr2::resp_body_string(resp))
      # El token contiene el username real del workspace: un=USERNAME|tokenid=...|
      ws_user  <- regmatches(token, regexpr("(?<=\\bun=)[^|]+", token, perl = TRUE))
      ws_user  <- if (length(ws_user) > 0 && nchar(ws_user) > 0)
                    utils::URLdecode(ws_user) else username
      list(success = TRUE, token = token, username = ws_user)
    } else {
      list(success = FALSE, error = paste0("Credenciales incorrectas (HTTP ", httr2::resp_status(resp), ")"))
    }
  }, error = function(e) {
    list(success = FALSE, error = paste("Error de red:", conditionMessage(e)))
  })
}

# ---- Crear directorio en workspace (ignora error si ya existe) ----
# Tupla de objeto: [path, type, user_metadata, data, creation_time] — 5 elementos
bvbrc_mkdir <- function(token, ws_dir) {
  rpc <- list(
    version = "1.1", method = "Workspace.create", id = "mkdir1",
    params  = list(list(
      objects   = list(list(ws_dir, "folder", setNames(list(), character(0)), NULL, NULL)),
      overwrite = 0L
    ))
  )
  tryCatch(
    httr2::request(BVBRC_WS_URL) |>
      httr2::req_method("POST") |>
      httr2::req_headers(
        "Authorization" = token,
        "Content-Type"  = "application/json",
        "Origin"        = "https://www.bv-brc.org",
        "Referer"       = "https://www.bv-brc.org/"
      ) |>
      httr2::req_body_json(rpc, auto_unbox = TRUE) |>
      httr2::req_timeout(15) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform(),
    error = function(e) NULL
  )
  invisible(NULL)
}

# ---- Upload de archivo al workspace (protocolo Shock de BV-BRC) ----
# Flujo documentado en WorkspaceClientExt.pm:
#   1. Workspace.create con createUploadNodes=1 → devuelve shock_url en meta[11] (0-based)
#   2. PUT del archivo al shock_url con Authorization: OAuth <token>, campo multipart "upload"
bvbrc_upload_file <- function(token, username, local_path, filename = NULL) {
  if (is.null(filename)) filename <- basename(local_path)
  # Usar el MISMO dominio que ws_home()/ws_out() (@bvbrc) para que el archivo
  # subido y el job que lo consume apunten al mismo workspace. Antes usaba
  # @patricbrc.org como fallback → los uploads iban a otro workspace que los jobs.
  uname   <- if (grepl("@", username, fixed = TRUE)) username else paste0(username, "@bvbrc")
  ws_dir  <- paste0("/", uname, "/home/.fastq_uploads")
  ws_path <- paste0(ws_dir, "/", filename)

  tryCatch({
    # Asegurar directorio
    bvbrc_mkdir(token, ws_dir)

    # Paso 1: crear entrada en workspace y obtener URL de Shock
    rpc_create <- list(
      version = "1.1",
      method  = "Workspace.create",
      id      = paste0("create_", format(Sys.time(), "%Y%m%d%H%M%OS3")),
      params  = list(list(
        objects           = list(list(ws_path, "reads", setNames(list(), character(0)), NULL, NULL)),
        createUploadNodes = 1L,
        overwrite         = 1L
      ))
    )

    resp_c <- httr2::request(BVBRC_WS_URL) |>
      httr2::req_method("POST") |>
      httr2::req_headers("Authorization" = token, "Content-Type" = "application/json") |>
      httr2::req_body_json(rpc_create, auto_unbox = TRUE) |>
      httr2::req_timeout(30) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    st_c  <- httr2::resp_status(resp_c)
    raw_c <- httr2::resp_body_string(resp_c)

    if (st_c != 200) {
      return(list(success = FALSE,
                  error   = paste("Workspace.create fallido:", substr(raw_c, 1, 200)),
                  log_det = paste0("HTTP ", st_c, " | ", substr(raw_c, 1, 300))))
    }

    # shock_url: result[[1]] = lista de objetos; [[1]] = primer objeto;
    # [[12]] = posición R (1-based) que corresponde a Perl [11] (0-based) = shock_node
    parsed_c  <- jsonlite::fromJSON(raw_c, simplifyVector = FALSE)
    obj_meta  <- parsed_c$result[[1]][[1]]
    shock_url <- obj_meta[[12]]

    if (is.null(shock_url) || nchar(trimws(as.character(shock_url))) == 0) {
      return(list(success = FALSE,
                  error   = "Workspace no devolvio URL de Shock en meta[11]",
                  log_det = paste0("result: ", substr(raw_c, 1, 400))))
    }

    # Paso 2: PUT del archivo al nodo Shock
    resp_s <- httr2::request(as.character(shock_url)) |>
      httr2::req_method("PUT") |>
      httr2::req_headers("Authorization" = paste0("OAuth ", token)) |>
      httr2::req_body_multipart(upload = curl::form_file(local_path)) |>
      httr2::req_timeout(600) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    st_s  <- httr2::resp_status(resp_s)
    raw_s <- httr2::resp_body_string(resp_s)
    det   <- paste0("Upload → ", ws_path,
                    "\nShock: ", shock_url,
                    "\nHTTP Shock ", st_s, " | ", substr(raw_s, 1, 200))

    if (st_s == 200) {
      list(success = TRUE, ws_path = ws_path, log_det = det)
    } else {
      list(success = FALSE,
           error   = paste("Shock upload fallido:", substr(raw_s, 1, 200)),
           log_det = det)
    }
  }, error = function(e) {
    list(success = FALSE,
         error   = paste("Error de upload:", conditionMessage(e)),
         log_det = paste("Excepcion:", conditionMessage(e)))
  })
}

# ---- Traducción de nombre interno → app_id real de BV-BRC ----
# Verificado con AppService.enumerate_apps: varios IDs difieren del nombre
# "comercial" de la herramienta. Enviar el ID equivocado hace que el job
# se acepte pero nunca produzca resultados.
BVBRC_APP_IDS <- c(
  Assembly2                   = "GenomeAssembly2",
  Annotation                  = "GenomeAnnotation",
  SeqComparison               = "GenomeComparison",
  PhylogeneticTree            = "CodonTree",     # CodonTree usa genome_ids
  FastqUtils                  = "FastqUtils",
  TaxonomicClassification     = "TaxonomicClassification",
  MetagenomicReadMapping      = "MetagenomicReadMapping",
  ComprehensiveGenomeAnalysis = "ComprehensiveGenomeAnalysis",
  Homology                    = "Homology",
  MSA                         = "MSA"
)
bvbrc_app_id <- function(tool) unname(BVBRC_APP_IDS[tool] %||% tool)

# ---- Envio de job (JSON-RPC 1.1) ----
# start_app(app_id, params, workspace) — el tercer argumento es el path del workspace
bvbrc_submit_job <- function(token, app_name, params, workspace) {
  app_name <- bvbrc_app_id(app_name)
  rpc_id   <- paste0("anbio_", format(Sys.time(), "%Y%m%d%H%M%OS3"))
  rpc_body <- list(
    method  = "AppService.start_app",
    params  = list(app_name, params, workspace),
    id      = rpc_id,
    version = "1.1"
  )

  req_json <- tryCatch(
    jsonlite::toJSON(rpc_body, auto_unbox = TRUE, pretty = FALSE),
    error = function(e) "{}"
  )

  tryCatch({
    resp <- httr2::request(BVBRC_APP_URL) |>   # URL base, sin /start_app/
      httr2::req_method("POST") |>
      httr2::req_headers(
        "Authorization" = token,
        "Content-Type"  = "application/json"
      ) |>
      httr2::req_body_json(rpc_body, auto_unbox = TRUE) |>
      httr2::req_timeout(60) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    status   <- httr2::resp_status(resp)
    raw_body <- httr2::resp_body_string(resp)

    # Construir detalle de log (sin token)
    log_det <- paste0(
      "POST ", BVBRC_APP_URL, "\n",
      "method: AppService.start_app  app: ", app_name, "\n",
      "HTTP ", status, " | resp: ", substr(raw_body, 1, 400)
    )

    if (status == 200) {
      parsed  <- tryCatch(
        jsonlite::fromJSON(raw_body, simplifyVector = FALSE),
        error = \(e) list()
      )
      # JSON-RPC 1.1: result[[1]] es el objeto del job (named list).
      # El UUID del task está en el campo $id; $monitor_url es la URL de seguimiento.
      raw_result  <- parsed$result[[1]]
      fields_info <- if (is.list(raw_result))
                       paste0("campos: ", paste(names(raw_result), collapse = ", "))
                     else
                       paste0("tipo: ", class(raw_result)[1])
      log_det <- paste0(log_det, "\n", fields_info)

      task_id <- if (is.list(raw_result)) {
                   r <- raw_result$id %||% raw_result$task_id %||% NULL
                   if (!is.null(r) && nchar(as.character(r)) > 0) as.character(r) else NULL
                 } else if (!is.null(raw_result) && length(raw_result) == 1) {
                   as.character(raw_result)
                 } else {
                   NULL
                 }

      if (!is.null(task_id) && nchar(task_id) > 0) {
        list(success  = TRUE,
             task_id  = task_id,
             app      = app_name,
             log_det  = log_det)
      } else {
        list(success = FALSE,
             error   = paste0("JSON-RPC OK pero sin task_id. ", fields_info,
                              ". Resp: ", substr(raw_body, 1, 300)),
             log_det = log_det)
      }
    } else {
      list(success = FALSE,
           error   = paste0("HTTP ", status, ": ", substr(raw_body, 1, 300)),
           log_det = log_det)
    }
  }, error = function(e) {
    list(success = FALSE,
         error   = paste("Error de red:", conditionMessage(e)),
         log_det = paste("Excepcion R:", conditionMessage(e)))
  })
}

# ---- ¿El resultado de un job es válido? ----
# Lee el objeto job_result (<output_file>) y revisa success + output_files.
# Devuelve: TRUE (resultado real), FALSE (falló o sin salida), NA (no verificable).
# Sirve para NO marcar como "completado" un job fallido/vacío (así se puede re-correr).
bvbrc_result_ok <- function(token, ws_result_path) {
  g <- tryCatch(bvbrc_get_file(token, ws_result_path),
                error = function(e) list(success = FALSE))
  if (!isTRUE(g$success) || !nchar(trimws(g$content %||% ""))) return(NA)
  j <- tryCatch(jsonlite::fromJSON(g$content, simplifyVector = FALSE),
                error = function(e) NULL)
  if (is.null(j)) return(NA)
  ok_success <- is.null(j$success) || isTRUE(as.integer(j$success) == 1L)
  n_out      <- length(j$output_files %||% list())
  isTRUE(ok_success) && n_out > 0
}

# ---- Estado del job ----
bvbrc_job_status <- function(token, task_id) {
  body <- list(
    method  = "AppService.query_task_summary",
    params  = list(list(task_id)),
    version = "1.1",
    id      = "1"
  )

  tryCatch({
    resp <- httr2::request(BVBRC_APP_URL) |>
      httr2::req_method("POST") |>
      httr2::req_headers(
        "Authorization" = token,
        "Content-Type"  = "application/json"
      ) |>
      httr2::req_body_json(body) |>
      httr2::req_timeout(15) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()

    st  <- httr2::resp_status(resp)
    raw <- httr2::resp_body_string(resp)

    if (st == 200) {
      parsed <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = FALSE), error = \(e) list())
      tasks  <- parsed$result[[1]]
      if (!is.null(tasks) && length(tasks) > 0) {
        t <- tasks[[1]]
        if (!is.list(t)) {
          list(success = TRUE, status = "queued",
               log_det = paste0("query_task_summary HTTP ", st,
                                " | respuesta inesperada (no lista): ", class(t)[1]))
        } else {
          params_p <- tryCatch({
            p <- t$parameters
            if (is.character(p)) jsonlite::fromJSON(p, simplifyVector = FALSE) else (p %||% list())
          }, error = \(e) list())
          out_path <- params_p$output_path %||% NULL
          out_file <- params_p$output_file %||% NULL
          app_name <- t$app %||% t$app_id %||% t$application %||% NULL

          # BV-BRC usa state_code de una letra (Q/R/C/F/E) o cadenas largas
          raw_st    <- t$status %||% t$state_code %||% "queued"
          mapped_st <- switch(as.character(raw_st),
            "C" = "completed", "F" = "failed",  "E" = "error",
            "R" = "running",   "Q" = "queued",
            "in-progress" = "running",
            "Completed"   = "completed",
            "Failed"      = "failed",
            "Running"     = "running",
            "Queued"      = "queued",
            as.character(raw_st)
          )

          list(
            success     = TRUE,
            status      = mapped_st,
            output_path = out_path,
            output_file = out_file,
            app_name    = app_name,
            log_det     = paste0("query_task_summary HTTP ", st,
                                 " | raw=", raw_st, " → ", mapped_st,
                                 " | app=", app_name %||% "?",
                                 " | out=", out_file %||% "?")
          )
        }
      } else {
        list(success = TRUE, status = "queued",
             log_det = paste0("query_task_summary HTTP ", st, " | sin tareas en result"))
      }
    } else {
      list(success = FALSE, status = "error",
           error   = paste0("HTTP ", st, ": ", substr(raw, 1, 200)),
           log_det = paste0("query_task_summary HTTP ", st, " | ", substr(raw, 1, 200)))
    }
  }, error = function(e) {
    list(success = FALSE, status = "error", error = conditionMessage(e),
         log_det = paste("Excepcion status:", conditionMessage(e)))
  })
}

# ---- Listar directorio del workspace (llamada de bajo nivel) ----
# Usamos JSON crudo (req_body_raw) para evitar cualquier problema de serialización
# con req_body_json / jsonlite. El path puede contener @ y / pero no chars que
# necesiten escape JSON, por lo que la interpolación directa es segura.
.bvbrc_ls_raw <- function(token, ws_path) {
  body_str <- paste0(
    '{"version":"1.1","id":"ls1","method":"Workspace.ls",',
    '"params":[{"paths":["', ws_path, '"],"adminmode":false}]}'
  )
  if (isTRUE(getOption("anbio.debug", FALSE))) message("[bvbrc_ls] → ", ws_path)
  tryCatch({
    resp <- httr2::request(BVBRC_WS_URL) |>
      httr2::req_method("POST") |>
      httr2::req_headers(
        "Authorization" = token,
        "Content-Type"  = "application/json",
        "Origin"        = "https://www.bv-brc.org",
        "Referer"       = "https://www.bv-brc.org/"
      ) |>
      httr2::req_body_raw(body_str, type = "application/json") |>
      httr2::req_timeout(20) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()
    st  <- httr2::resp_status(resp)
    raw <- httr2::resp_body_string(resp)
    if (isTRUE(getOption("anbio.debug", FALSE)))
      message("[bvbrc_ls] ← HTTP ", st, " | ", substr(raw, 1, 200))

    # BV-BRC puede devolver errores JSON-RPC en HTTP 200 O HTTP 500
    parsed <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = FALSE),
                       error = \(e) list())
    if (!is.null(parsed$error)) {
      raw_msg <- parsed$error$error %||% parsed$error$message %||%
                 parsed$error$data  %||% substr(raw, 1, 200)
      err_msg <- gsub("_ERROR_", "", as.character(raw_msg))
      return(list(success = FALSE,
                  error   = paste0("HTTP ", st, " — ", trimws(err_msg))))
    }
    if (st != 200)
      return(list(success = FALSE,
                  error   = paste0("HTTP ", st, ": ", substr(raw, 1, 200))))

    dir_map <- parsed$result[[1]] %||% list()

    # BV-BRC puede devolver la clave sin trailing slash o con él
    objs <- dir_map[[ws_path]] %||%
            dir_map[[paste0(ws_path, "/")]] %||%
            (if (length(dir_map) == 1) dir_map[[1]] else list())

    files <- Filter(Negate(is.null), lapply(objs %||% list(), function(obj) {
      if (!is.list(obj) || length(obj) < 1) return(NULL)
      name <- as.character(obj[[1]] %||% "")
      if (!nchar(name)) return(NULL)
      type <- as.character(obj[[2]] %||% "file")
      size <- suppressWarnings(as.numeric(obj[[7]] %||% 0)) %||% 0
      list(name = name, type = type, size = size,
           path = paste0(ws_path, "/", name))
    }))
    list(success = TRUE, files = files, path = ws_path)
  }, error = function(e)
    list(success = FALSE, error = conditionMessage(e)))
}

# ---- Wrapper público: reintenta con dominio alternativo si el primero falla ----
# BV-BRC acepta @patricbrc.org y @bvbrc como alias del mismo espacio.
# Algunos tokens dan permisos solo sobre uno de los dos; probamos ambos.
bvbrc_ls <- function(token, ws_path) {
  r <- .bvbrc_ls_raw(token, ws_path)
  if (isTRUE(r$success)) return(r)
  alt <- if (grepl("@patricbrc.org", ws_path, fixed = TRUE))
           gsub("@patricbrc.org", "@bvbrc",       ws_path, fixed = TRUE)
         else
           gsub("@bvbrc",        "@patricbrc.org", ws_path, fixed = TRUE)
  if (alt == ws_path) return(r)
  r2 <- .bvbrc_ls_raw(token, alt)
  if (isTRUE(r2$success)) r2 else r
}

# ---- Descargar contenido texto de un archivo del workspace ----
bvbrc_get_file <- function(token, ws_path) {
  body_str <- paste0(
    '{"version":"1.1","id":"getf1","method":"Workspace.get",',
    '"params":[{"objects":["', ws_path, '"]}]}'
  )
  tryCatch({
    resp <- httr2::request(BVBRC_WS_URL) |>
      httr2::req_method("POST") |>
      httr2::req_headers(
        "Authorization" = token,
        "Content-Type"  = "application/json",
        "Origin"        = "https://www.bv-brc.org",
        "Referer"       = "https://www.bv-brc.org/"
      ) |>
      httr2::req_body_raw(body_str, type = "application/json") |>
      httr2::req_timeout(60) |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_perform()
    st  <- httr2::resp_status(resp)
    raw <- httr2::resp_body_string(resp)

    parsed  <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = FALSE),
                        error = \(e) list())
    if (!is.null(parsed$error)) {
      raw_msg <- parsed$error$error %||% parsed$error$message %||% substr(raw, 1, 200)
      return(list(success = FALSE,
                  error   = paste0("HTTP ", st, " — ",
                                   trimws(gsub("_ERROR_", "", as.character(raw_msg))))))
    }
    if (st != 200)
      return(list(success = FALSE,
                  error   = paste0("HTTP ", st, ": ", substr(raw, 1, 200))))
    # result[[1]] = lista de pares [metadata_tuple, data_string]
    pair    <- parsed$result[[1]][[1]] %||% NULL
    content <- if (!is.null(pair) && length(pair) >= 2) pair[[2]] %||% "" else ""
    content <- as.character(content %||% "")
    # Archivos grandes: BV-BRC los guarda en Shock y Workspace.get devuelve la
    # URL del nodo (no el contenido). Hay que descargar los bytes desde Shock.
    if (grepl("^https?://.*shock_api/node/", trimws(content))) {
      shock_url <- trimws(content)
      dl <- httr2::request(paste0(shock_url, "?download")) |>
        httr2::req_method("GET") |>
        httr2::req_headers("Authorization" = paste0("OAuth ", token)) |>
        httr2::req_timeout(120) |>
        httr2::req_error(is_error = \(r) FALSE) |>
        httr2::req_perform()
      if (httr2::resp_status(dl) == 200)
        content <- httr2::resp_body_string(dl)
      else
        return(list(success = FALSE,
                    error = paste0("Shock download HTTP ", httr2::resp_status(dl))))
    }
    list(success = TRUE, content = content)
  }, error = function(e)
    list(success = FALSE, error = conditionMessage(e)))
}

# ============================================================
# Constructores de parametros por herramienta
# ============================================================

ws_home <- function(username) {
  uname <- if (grepl("@", username, fixed = TRUE)) username
             else paste0(username, "@bvbrc")
  paste0("/", uname, "/home")
}
# Carpeta de salida para jobs de esta app.
# NO usar ".output" — es una carpeta de sistema solo accesible por el App Service.
ws_out <- function(username) {
  uname <- if (grepl("@", username, fixed = TRUE)) username
             else paste0(username, "@bvbrc")
  paste0("/", uname, "/home/PRUEBAAPP")
}

# Construye campos de librería paired o single según disponibilidad de R2
bvbrc_lib_fields <- function(r1_ws, r2_ws = NULL) {
  has_r2 <- !is.null(r2_ws) && nchar(r2_ws %||% "") > 0 && r2_ws != r1_ws
  if (has_r2) {
    list(paired_end_libs = list(list(
      read1 = r1_ws, read2 = r2_ws, platform = "illumina", interleaved = FALSE
    )))
  } else {
    list(single_end_libs = list(list(
      read = r1_ws, platform = "illumina", interleaved = FALSE
    )))
  }
}

bvbrc_params_fastqutil <- function(username, r1_ws, r2_ws = NULL, sample = "muestra") {
  # 'recipe' es REQUERIDO: sin él el job corre pero no genera salida
  # (output_files vacío). "fastqc" ejecuta el control de calidad tipo FastQC.
  c(bvbrc_lib_fields(r1_ws, r2_ws),
    list(recipe      = list("fastqc"),
         output_path = ws_out(username),
         output_file = paste0("fastqc_", sample)))
}

bvbrc_params_taxonomy <- function(username, r1_ws, r2_ws = NULL, sample = "muestra") {
  # TaxonomicClassification NO acepta 'algorithm'. Requiere 'analysis_type'
  # y 'host_genome'. Para aislados bacterianos: analysis_type="pathogen"
  # con database="bvbrc". (16S requeriría Greengenes/SILVA, no bvbrc → fallaba.)
  c(bvbrc_lib_fields(r1_ws, r2_ws),
    list(analysis_type = "pathogen",
         host_genome   = "no_host",
         database      = "bvbrc",
         output_path   = ws_out(username),
         output_file   = paste0("taxonomy_", sample)))
}

bvbrc_params_assembly <- function(username, r1_ws, r2_ws = NULL,
                                   sample = "muestra", strategy = "auto") {
  c(bvbrc_lib_fields(r1_ws, r2_ws),
    list(recipe      = strategy,
         output_path = ws_out(username),
         output_file = paste0("assembly_", sample)))
}

# GenomeAnnotation: 'domain' y 'code' son REQUERIDOS (enum). 'recipe' es
# string libre; se omite para usar el pipeline por defecto (RASTtk).
bvbrc_params_annotation <- function(username, contigs_ws, organism, sample = "muestra") {
  list(
    contigs         = contigs_ws,
    scientific_name = if (nchar(trimws(organism %||% "")) > 0) organism else "Bacteria",
    taxonomy_id     = 2,
    code            = 11,
    domain          = "Bacteria",
    output_path     = ws_out(username),
    output_file     = paste0("annotation_", sample)
  )
}

# ComprehensiveGenomeAnalysis: 'input_type' y 'domain' son REQUERIDOS.
bvbrc_params_cga <- function(username, r1_ws, r2_ws = NULL, organism, sample = "muestra") {
  c(bvbrc_lib_fields(r1_ws, r2_ws),
    list(input_type      = "reads",
         recipe          = "auto",
         scientific_name = if (nchar(trimws(organism %||% "")) > 0) organism else "Bacteria",
         taxonomy_id     = 2,
         code            = 11,
         domain          = "Bacteria",
         output_path     = ws_out(username),
         output_file     = paste0("cga_", sample)))
}

bvbrc_params_metareads <- function(username, r1_ws, r2_ws = NULL, sample = "muestra") {
  c(bvbrc_lib_fields(r1_ws, r2_ws),
    list(gene_set_name   = "CARD",
         gene_set_type   = "predefined_list",
         output_path     = ws_out(username),
         output_file     = paste0("resistoma_", sample)))
}

bvbrc_params_blast <- function(username, sequence, program = "blastn",
                                db = "PATRIC", sample = "blast") {
  list(
    program     = program,
    db_type     = "fasta",
    database    = db,
    query       = sequence,
    evalue      = "1e-5",
    max_hits    = 50,
    output_path = ws_out(username),
    output_file = paste0("blast_", sample)
  )
}

bvbrc_params_mlst <- function(username, genome_ids, output_name = "mlst") {
  list(
    genome_ids  = as.list(genome_ids),
    output_path = ws_out(username),
    output_file = output_name
  )
}

# CodonTree (app real de PhylogeneticTree): usa 'genome_ids' + bootstraps.
# PhylogeneticTree "puro" exige in_genome_ids Y out_genome_ids (grupo externo),
# que no tenemos; CodonTree es el árbol estándar de genes core de BV-BRC.
bvbrc_params_phylogeny <- function(username, genome_ids, output_name = "filogenia") {
  list(
    genome_ids      = as.list(genome_ids),
    number_of_genes = 100,
    bootstraps      = 100,
    output_path     = ws_out(username),
    output_file     = output_name
  )
}

# GenomeComparison: pese a que enumerate_apps declara 'genome_ids' como string,
# el wrapper Perl (App-GenomeComparison.pl) lo desreferencia como ARRAY.
# Debe enviarse como lista, NO como cadena separada por comas.
bvbrc_params_proteome <- function(username, genome_ids, output_name = "proteoma") {
  list(
    genome_ids             = as.list(as.character(genome_ids)),
    # reference_genome_index es 1-based y NO puede quedar indefinido:
    # sin él el wrapper falla con "Not an ARRAY reference" (línea 523).
    reference_genome_index = 1L,
    output_path            = ws_out(username),
    output_file            = output_name
  )
}

# ---- Mapa app_name por etapa ----
BVBRC_APPS <- list(
  calidad     = "FastqUtils",
  taxonomia   = "TaxonomicClassification",
  ensamblado  = "Assembly2",
  anotacion   = "Annotation",
  resistoma   = list("ComprehensiveGenomeAnalysis", "MetagenomicReadMapping"),
  alineamiento = "Homology",
  mlst        = "CoreGenomeMLST",
  filogenia   = "PhylogeneticTree",
  genomica    = "SeqComparison"
)
