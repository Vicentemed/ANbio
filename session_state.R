# ============================================================
# Persistencia de sesión — guarda y restaura el estado de rv
# Los jobs de BV-BRC/Galaxy siguen corriendo en el servidor
# externo; sólo necesitamos los task_id y metadatos para
# reconectarnos cuando la app se reinicie.
# ============================================================

SESSION_STATE_DIR  <- file.path(path.expand("~"), ".anbio_session")
SESSION_STATE_FILE <- file.path(SESSION_STATE_DIR, "anbio_state.rds")

dir.create(SESSION_STATE_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Guardar estado ----
# Llama desde contexto reactivo; lee rv directamente.
# Los raw file contents (rv$job_results[[x]]$contents) se excluyen
# para mantener el archivo pequeño (sólo metadatos y parsed).
state_save <- function(rv) {
  tryCatch({
    # Versión slim de job_results: sin el campo "contents" (strings grandes)
    jr_slim <- lapply(rv$job_results, function(r) {
      r[setdiff(names(r), "contents")]
    })

    saveRDS(list(
      saved_at     = Sys.time(),
      version      = "1.1",
      jobs         = rv$jobs,
      job_results  = jr_slim,
      calidad      = rv$calidad,
      taxonomia    = rv$taxonomia,
      ensamblado   = rv$ensamblado,
      mlst         = rv$mlst,
      alineamiento = rv$alineamiento,
      ani          = rv$ani,
      resfinder_df = rv$resfinder_df,
      ebi_result   = rv$ebi_result,
      muestras     = rv$muestras,
      fastq_dir    = rv$fastq_dir,
      fastq_files  = rv$fastq_files,
      completados  = rv$completados,
      bvbrc_user   = rv$bvbrc_user
    ), SESSION_STATE_FILE)
  }, error = function(e) warning("state_save: ", e$message))
  invisible(NULL)
}

# ---- Cargar estado ----
# Devuelve NULL si no existe o está corrupto.
# Descarta estados con más de 7 días (jobs expirados).
state_load <- function() {
  if (!file.exists(SESSION_STATE_FILE)) return(NULL)
  tryCatch({
    st <- readRDS(SESSION_STATE_FILE)
    age_days <- as.numeric(difftime(Sys.time(), st$saved_at, units = "days"))
    if (age_days > 7) { unlink(SESSION_STATE_FILE); return(NULL) }
    st
  }, error = function(e) {
    warning("state_load: ", e$message)
    NULL
  })
}

# ---- Borrar estado ----
state_clear <- function() {
  if (file.exists(SESSION_STATE_FILE)) unlink(SESSION_STATE_FILE)
  invisible(NULL)
}

# ============================================================
# Archivo de sesión portable (.anbio)
# Permite cerrar la app y reanudar el análisis después: guarda
# los task_id de los jobs (que siguen corriendo en los servidores)
# junto con los resultados ya importados y los metadatos del proyecto.
# ============================================================

SESSION_DIR <- file.path(path.expand("~"), ".anbio_session", "sesiones")
dir.create(SESSION_DIR, showWarnings = FALSE, recursive = TRUE)

# Construye la lista serializable a partir de rv + metadatos del proyecto
session_payload <- function(rv, proyecto = list()) {
  jr_slim <- lapply(rv$job_results, function(r) r[setdiff(names(r), "contents")])
  list(
    formato      = "anbio-session",
    version      = "2.0",
    saved_at     = Sys.time(),
    proyecto     = proyecto,
    bvbrc_user   = rv$bvbrc_user,
    muestras     = rv$muestras,
    fastq_dir    = rv$fastq_dir,
    fastq_files  = rv$fastq_files,
    jobs         = rv$jobs,
    job_results  = jr_slim,
    completados  = rv$completados,
    auto_pipe    = rv$auto_pipe,
    calidad      = rv$calidad,
    taxonomia    = rv$taxonomia,
    ensamblado   = rv$ensamblado,
    mlst         = rv$mlst,
    alineamiento = rv$alineamiento,
    ani          = rv$ani,
    resfinder_df = rv$resfinder_df,
    amr_sir      = rv$amr_sir,
    anotacion    = rv$anotacion,
    pangenoma    = rv$pangenoma,
    newick       = rv$newick
  )
}

# Guardar a una ruta concreta (.anbio)
session_save_file <- function(rv, path, proyecto = list()) {
  tryCatch({
    saveRDS(session_payload(rv, proyecto), path)
    TRUE
  }, error = function(e) { warning("session_save_file: ", e$message); FALSE })
}

# Cargar desde archivo; valida el formato
session_load_file <- function(path) {
  if (!file.exists(path)) return(NULL)
  st <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(st) || !identical(st$formato, "anbio-session")) return(NULL)
  st
}

# Ruta del autoguardado de la sesión activa
session_autosave_path <- function(nombre = "sesion_activa") {
  file.path(SESSION_DIR, paste0(gsub("[^A-Za-z0-9_-]", "_", nombre), ".anbio"))
}

# Listar sesiones guardadas (más reciente primero)
session_list <- function() {
  fs <- list.files(SESSION_DIR, pattern = "\\.anbio$", full.names = TRUE)
  if (!length(fs)) return(data.frame())
  info <- file.info(fs)
  df <- data.frame(path = fs, nombre = tools::file_path_sans_ext(basename(fs)),
                   fecha = info$mtime, kb = round(info$size / 1024, 1),
                   stringsAsFactors = FALSE)
  df[order(df$fecha, decreasing = TRUE), ]
}

# ---- Resumen legible del estado ----
state_summary <- function(st) {
  n_jobs   <- length(st$jobs %||% list())
  n_active <- sum(vapply(st$jobs %||% list(), function(j)
    isTRUE(j$status %in% c("queued","running","Q","R",
                             "in_progress","in-progress","submitted")),
    logical(1)))
  n_done <- sum(vapply(st$jobs %||% list(), function(j)
    isTRUE(j$status %in% c("completed","FINISHED","ok","done")),
    logical(1)))
  list(n_jobs = n_jobs, n_active = n_active, n_done = n_done,
       has_calidad   = !is.null(st$calidad)    && nrow(st$calidad)    > 0,
       has_taxonomia = !is.null(st$taxonomia)  && nrow(st$taxonomia)  > 0,
       has_ensamblado= !is.null(st$ensamblado) && nrow(st$ensamblado) > 0,
       n_muestras    = length(st$muestras %||% character(0)))
}
