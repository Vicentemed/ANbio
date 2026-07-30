# ============================================================
# SERVER — lógica principal de la app ANbio
# ============================================================

# Carpeta local con los archivos FASTQ del proyecto
DATOS_DIR <- "C:/Escritorio 2026/CARPETAS/Tesis Cristina/ANbio/Datos"

server <- function(input, output, session) {

  rv <- reactiveValues(
    muestras     = character(0),
    fastq_dir    = NULL,
    fastq_files  = character(0),
    staged_files = list(),
    completados  = setNames(rep(FALSE, length(PASOS_IDS)), PASOS_IDS),
    bvbrc_token  = NULL,
    bvbrc_user   = NULL,
    galaxy_key   = NULL,
    galaxy_hist  = NULL,
    jobs         = list(),
    calidad      = data.frame(Muestra=character(),Q30=numeric(),GC=numeric(),
                               Reads_M=numeric(),Duplicados=numeric(),Estado=character(),
                               stringsAsFactors=FALSE),
    taxonomia    = data.frame(Muestra=character(),Especie=character(),
                               Confianza=numeric(),Contaminacion=numeric(),
                               stringsAsFactors=FALSE),
    ensamblado   = data.frame(Muestra=character(),Cobertura=numeric(),
                               Profundidad=numeric(),Contigs=numeric(),N50_kb=numeric(),
                               Tamano_Mb=numeric(),GC=numeric(),stringsAsFactors=FALSE),
    mlst         = data.frame(Muestra=character(),ST=character(),Alelos=character(),
                               Esquema=character(),stringsAsFactors=FALSE),
    alineamiento = data.frame(Par=character(),Herramienta=character(),
                               Identidad=numeric(),Cobertura=numeric(),
                               Gaps=numeric(),SNPs=numeric(),stringsAsFactors=FALSE),
    ani          = data.frame(Muestra1=character(),Muestra2=character(),
                               ANI=numeric(),stringsAsFactors=FALSE),
    reporte_path    = NULL,
    ebi_result      = NULL,
    resfinder_df    = NULL,
    blast_rid       = NULL,
    log_entries          = list(),
    job_results          = list(),   # job_id -> resultado parseado
    viewing_job_id       = NULL,     # job mostrado en modal
    pending_restore      = NULL,     # estado cargado pendiente de confirmar
    backup_import_status = NULL,     # mensaje tras importar JSON backup
    # ── Resultados adicionales para el reporte final ──
    amr_sir     = NULL,   # data.frame  Antibiotico | Prediccion | Confianza (por muestra)
    anotacion   = NULL,   # data.frame  Muestra | CDS | rRNA | tRNA | ...
    pangenoma   = NULL,   # list(core, accesorios, unicos, total)
    newick      = NULL,   # cadena Newick del árbol
    sesion_nombre = "sesion_activa",
    sesion_msg    = NULL,
    reporte_final_path = NULL,
    auto_pipe            = list(     # orquestador del pipeline automático 1→9
      active     = FALSE,
      samples    = character(0),
      organism   = "Unknown",
      ws         = list(),   # muestra -> list(r1=..., r2=...)
      submitted  = list(),   # muestra -> character(): tools ya enviados
      contigs    = list(),   # muestra -> ruta ws de contigs (tras ensamblado)
      genome_ids = list(),   # muestra -> genome_id (tras anotación)
      galaxy_ds  = list(),   # muestra -> dataset_id de Galaxy (para Kraken2)
      global     = character(0),  # etapas globales ya lanzadas (filogenia, etc.)
      started_at = NULL
    )
  )

  # ── Helper de logging reactivo ─────────────────────────────
  # Escribe al archivo Y agrega al rv para la UI
  alog <- function(nivel, msg, det = NULL) {
    entry <- app_log(nivel, msg, det)          # app_log() de R/logger.R
    rv$log_entries <- c(rv$log_entries, list(entry))
    if (length(rv$log_entries) > 400)
      rv$log_entries <- tail(rv$log_entries, 400)
    invisible(entry)
  }

  # ─── Cargar FASTQ desde una carpeta ────────────────────────
  # Detecta pares R1/R2, deriva los nombres de muestra y avisa de
  # descargas incompletas (.part). Devuelve el nº de muestras.
  cargar_fastq_de_carpeta <- function(dir_path, verbose = TRUE) {
    if (is.null(dir_path) || !nzchar(dir_path) || !dir.exists(dir_path)) return(0L)
    all_fastq <- list.files(dir_path, pattern = "\\.f(ast)?q(\\.gz)?$",
                            ignore.case = TRUE, full.names = FALSE)
    part_files <- list.files(dir_path, pattern = "\\.part$", ignore.case = TRUE)
    if (!length(all_fastq)) {
      if (verbose) {
        alog("WARN", "Carpeta sin archivos FASTQ", dir_path)
        showNotification("La carpeta no contiene archivos FASTQ.", type = "warning", duration = 6)
      }
      return(0L)
    }
    nombres <- sub("_R[12].*$", "",
                sub("_L00[0-9].*$", "", all_fastq, ignore.case = TRUE), ignore.case = TRUE)
    nombres <- sub("\\.f(ast)?q(\\.gz)?$", "", nombres, ignore.case = TRUE)
    nombres <- sort(unique(nombres[nchar(nombres) > 0]))
    rv$fastq_dir   <- dir_path
    rv$fastq_files <- all_fastq
    rv$muestras    <- nombres
    if (verbose) {
      alog("OK", paste0("FASTQ cargados desde: ", dir_path),
           paste0(length(all_fastq), " archivo(s) → ", length(nombres),
                  " muestra(s): ", paste(nombres, collapse = ", "),
                  if (length(part_files))
                    paste0("\nINCOMPLETOS (.part): ", paste(part_files, collapse = ", ")) else ""))
      showNotification(
        paste0(length(nombres), " muestra(s) detectada(s): ",
               paste(head(nombres, 5), collapse = ", "),
               if (length(nombres) > 5) "..." else ""),
        type = "message", duration = 6)
      if (length(part_files))
        showNotification(
          paste0("Archivos incompletos (.part): ", paste(part_files, collapse = ", "),
                 " — vuelve a descargarlos."), type = "warning", duration = 12)
    }
    length(nombres)
  }

  # ─── Timers ────────────────────────────────────────────────
  poll_timer    <- reactiveTimer(45000)   # sondeo de estado de jobs
  sidebar_timer <- reactiveTimer(15000)   # refresca elapsed-time en sidebar

  # ─── Sidebar: menuItem dinámico con spinner por etapa ──────
  for (.si in seq_along(PASOS_IDS)) {
    local({
      step <- PASOS_IDS[.si]
      num  <- .si
      lbl  <- PASOS_LABELS[.si]
      icn  <- PASOS_ICONS[.si]
      txt  <- paste0(num, ". ", lbl)

      output[[paste0("sidebar_step_", step)]] <- renderMenu({
        jobs_step <- Filter(function(j) j$step == step, rv$jobs)

        chk <- function(jobs, statuses)
          sum(vapply(jobs, function(j) isTRUE(j$status %in% statuses), logical(1)))

        n_active <- chk(jobs_step, c("queued","running","Q","R",
                                      "in_progress","in-progress","submitted"))
        n_error  <- chk(jobs_step, c("error","failed","ERROR"))
        n_done   <- chk(jobs_step, c("completed","FINISHED","ok","done"))

        if (n_active > 0) {
          menuItem(txt, tabName = step,
            icon       = icon("spinner", class = "fa-spin"),
            badgeLabel = if (n_active > 1) paste0(n_active, " jobs") else "activo",
            badgeColor = "yellow")
        } else if (n_error > 0 && n_done == 0) {
          menuItem(txt, tabName = step, icon = icon(icn),
            badgeLabel = "error", badgeColor = "red")
        } else if (n_done > 0 && n_error > 0) {
          menuItem(txt, tabName = step, icon = icon(icn),
            badgeLabel = paste0(n_done, " ok"), badgeColor = "orange")
        } else if (n_done > 0) {
          menuItem(txt, tabName = step, icon = icon(icn),
            badgeLabel = paste0(n_done, " ok"), badgeColor = "green")
        } else {
          menuItem(txt, tabName = step, icon = icon(icn))
        }
      })
    })
  }

  # ─── Panel lateral de jobs activos / recientes ─────────────
  output$sidebar_active_jobs <- renderUI({
    sidebar_timer()   # refresca elapsed cada 15 s

    jobs_all <- rv$jobs
    active <- Filter(function(j)
      j$status %in% c("queued","running","Q","R",
                       "in_progress","in-progress","submitted"),
      jobs_all)
    recent_done <- Filter(function(j)
      j$status %in% c("completed","FINISHED","ok","done") &&
      as.numeric(difftime(Sys.time(), j$time, units = "mins")) < 5,
      jobs_all)

    items <- c(active, recent_done)
    if (!length(items)) return(NULL)

    step_label <- function(s) {
      idx <- which(PASOS_IDS == s)
      if (length(idx)) PASOS_LABELS[idx] else s
    }

    div(class = "sidebar-jobs-panel",
      div(class = "sj-header",
        tags$i(class = "fa fa-tasks"),
        paste0("En proceso (", length(active), "/", length(items), ")")
      ),
      lapply(items, function(j) {
        is_done  <- j$status %in% c("completed","FINISHED","ok","done")
        is_err   <- j$status %in% c("error","failed","ERROR")
        elapsed  <- as.integer(difftime(Sys.time(), j$time, units = "mins"))
        etxt     <- if (elapsed < 1) "<1m" else paste0(elapsed, "m")
        row_cls  <- paste("sj-row", if (is_done) "sj-done" else if (is_err) "sj-err" else "")

        div(class = row_cls,
          div(class = "sj-icon",
            if (is_done)
              tags$i(class = "fa fa-check-circle")
            else if (is_err)
              tags$i(class = "fa fa-times-circle")
            else
              tags$i(class = "fa fa-spinner fa-spin")
          ),
          div(class = "sj-label",
            title = j$label,
            tags$span(class = "sj-step", step_label(j$step)), " ",
            tags$span(j$tool %||% "")
          ),
          tags$span(class = "sj-time", etxt)
        )
      })
    )
  })

  observe({
    poll_timer()
    jobs <- isolate(rv$jobs)
    if (length(jobs) == 0) return()
    for (i in seq_along(jobs)) {
      j <- jobs[[i]]

      # ── Reintento automático ─────────────────────────────────────────────────
      # Si el job está marcado como completado (FINISHED/completed) pero aún no
      # tiene resultados descargados, intentar la descarga en este ciclo.
      # Cubre: (1) fetch falló la primera vez por red/workspace no listo,
      # (2) app cerrada justo después de recibir FINISHED pero antes de descargar.
      # Se reintenta cada 45 s hasta éxito; luego la condición deja de cumplirse.
      if (j$platform %in% c("bvbrc","galaxy") &&
          j$status %in% c("FINISHED","ok","done","completed") &&
          is.null(rv$job_results[[j$job_id]])) {
        if (j$platform == "bvbrc" &&
            !is.null(j$output_path) && !is.null(j$output_file)) {
          tok_l <- isolate(rv$bvbrc_token)
          if (!is.null(tok_l)) {
            fetch_res <- tryCatch(
              bvbrc_fetch_results(tok_l, j),
              error = function(e) list(success = FALSE, error = conditionMessage(e))
            )
            if (isTRUE(fetch_res$success)) {
              rv$job_results[[j$job_id]] <- fetch_res
              auto_import_results(fetch_res, j)
              rv$jobs[[i]]$status <- "completed"
              alog("OK", paste0("BV-BRC resultados recuperados (reintento): ", j$label),
                   paste0("Archivos: ", length(fetch_res$files)))
              tryCatch(state_save(rv), error = function(e) NULL)
            } else {
              # Límite de reintentos: tras 3 fallos, dejar de re-consultar
              # (evita el spam de Workspace.ls y el bloqueo de la UI cuando
              #  el job falló en BV-BRC o no dejó salida).
              tries <- (rv$jobs[[i]]$fetch_tries %||% 0L) + 1L
              rv$jobs[[i]]$fetch_tries <- tries
              if (tries >= 3L) {
                rv$jobs[[i]]$status <- "sin_resultados"
                alog("WARN", paste0("BV-BRC sin resultados tras 3 intentos: ", j$label),
                     paste0(fetch_res$error %||% "El job pudo fallar en BV-BRC. ",
                            "Revísalo en la web; se puede re-ejecutar."))
                tryCatch(state_save(rv), error = function(e) NULL)
              }
            }
          }
        } else if (j$platform == "galaxy") {
          key_l <- isolate(rv$galaxy_key)
          if (!is.null(key_l)) {
            fetch_res <- tryCatch(
              galaxy_fetch_results(key_l, j),
              error = function(e) list(success = FALSE, error = conditionMessage(e))
            )
            if (isTRUE(fetch_res$success)) {
              rv$job_results[[j$job_id]] <- fetch_res
              auto_import_results(fetch_res, j)
              rv$jobs[[i]]$status <- "completed"
              alog("OK", paste0("Galaxy resultados recuperados (reintento): ", j$label))
              tryCatch(state_save(rv), error = function(e) NULL)
            }
          }
        }
        next  # no hacer status-poll en este ciclo para este job
      }

      if (!j$status %in% c("queued","running","RUNNING","QUEUED","submitted",
                            "in_progress","in-progress","Q","R",
                            "new","waiting")) next   # "new"/"waiting" son estados Galaxy
      poll_res <- tryCatch({
        switch(j$platform,
          bvbrc  = { tok <- isolate(rv$bvbrc_token); if (!is.null(tok)) bvbrc_job_status(tok, j$job_id) else list(status=j$status) },
          galaxy = { key <- isolate(rv$galaxy_key);  if (!is.null(key)) galaxy_job_status(key, j$job_id) else list(status=j$status) },
          ebi    = ebi_status(j$tool, j$job_id),
          ncbi   = ncbi_blast_status(j$job_id),
          list(status = j$status))
      }, error = function(e) list(status = j$status, error = conditionMessage(e)))

      new_st <- poll_res$status %||% j$status

      # Normalizar estados de expiración devueltos por cada plataforma
      # EBI devuelve "NOT_FOUND" (>7 días) / NCBI devuelve "NOT_FOUND" (>36h)
      if (new_st %in% c("NOT_FOUND", "UNKNOWN", "FAILURE")) {
        elapsed_h <- as.numeric(difftime(Sys.time(), j$time, units = "hours"))
        new_st <- if (elapsed_h > 1) "expirado" else "error"
      }

      # Fallback BV-BRC: si query_task_summary sigue devolviendo queued/running
      # pero el job lleva >5 minutos, verificar directamente si el output existe.
      # Esto cubre el caso en que la API no refleja el estado real del job.
      if (j$platform == "bvbrc" &&
          new_st %in% c("queued","running","Q","R","submitted","in_progress") &&
          !is.null(j$output_path) && !is.null(j$output_file)) {
        elapsed_min <- as.numeric(difftime(Sys.time(), j$time, units = "mins"))
        if (elapsed_min > 5) {
          tok_chk <- isolate(rv$bvbrc_token)
          if (!is.null(tok_chk)) {
            out_dir_chk <- paste0(j$output_path, "/", j$output_file)
            ls_chk <- tryCatch(
              bvbrc_ls(tok_chk, out_dir_chk),
              error = function(e) list(success = FALSE))
            if (isTRUE(ls_chk$success) && length(ls_chk$files) > 0) {
              new_st <- "FINISHED"
              alog("OK", paste0("[BV-BRC] Output detectado via workspace (fallback): ", j$label),
                   paste0(length(ls_chk$files), " archivo(s) en ", out_dir_chk))
            }
          }
        }
      }

      # Log solo cuando cambia el estado
      if (new_st != j$status) {
        det_poll <- poll_res$log_det %||% poll_res$error %||% ""
        if (new_st %in% c("FINISHED","ok","completed","done"))
          alog("OK",   paste0("[", toupper(j$platform), "] ", j$label, ": ", j$status, " → ", new_st), det_poll)
        else if (new_st %in% c("error","ERROR","failed"))
          alog("ERROR",paste0("[", toupper(j$platform), "] ", j$label, ": ", j$status, " → ", new_st), det_poll)
        else if (new_st == "expirado")
          alog("WARN", paste0("[", toupper(j$platform), "] ", j$label,
                              ": RID/job_id expirado — no recuperable"), det_poll)
        else
          alog("INFO", paste0("[", toupper(j$platform), "] ", j$label, ": ", j$status, " → ", new_st), det_poll)
      }

      rv$jobs[[i]]$status <- new_st

      # Propagar metadata de workspace desde la consulta de estado (cubre jobs
      # que se guardaron sin output_path/output_file o con tool vacío)
      if (j$platform == "bvbrc") {
        if (is.null(rv$jobs[[i]]$output_path) && !is.null(poll_res$output_path))
          rv$jobs[[i]]$output_path <- poll_res$output_path
        if (is.null(rv$jobs[[i]]$output_file) && !is.null(poll_res$output_file))
          rv$jobs[[i]]$output_file <- poll_res$output_file
        if ((!nchar(rv$jobs[[i]]$tool %||% "")) && !is.null(poll_res$app_name))
          rv$jobs[[i]]$tool <- poll_res$app_name
        # Actualizar j local también para que el bloque de descarga lo use
        j <- rv$jobs[[i]]
      }

      if (new_st %in% c("FINISHED","ok","completed","done")) {
        tryCatch({
          if (j$platform == "ebi") {
            res <- ebi_result(j$tool, j$job_id, "aln-fasta")
            if (res$success) {
              out <- file.path(SESSION_TMPDIR, paste0(j$tool,"_",j$job_id,".fasta"))
              writeLines(res$content, out)
              rv$jobs[[i]]$result_file <- out
              rv$jobs[[i]]$result_text <- res$content
              if (j$step == "alineamiento") rv$ebi_result <- res$content
            }
          } else if (j$platform == "ncbi") {
            res <- ncbi_blast_results(j$job_id)
            if (res$success) { rv$jobs[[i]]$result_text <- res$content; rv$blast_rid <- j$job_id }
          } else if (j$platform == "cge") {
            res <- resfinder_result(j$job_id)
            if (res$success) rv$resfinder_df <- resfinder_to_table(res$parsed)
          } else if (j$platform == "bvbrc") {
            tok_l <- isolate(rv$bvbrc_token)
            if (!is.null(tok_l) && !is.null(j$output_path) && !is.null(j$output_file)) {
              fetch_res <- tryCatch(
                bvbrc_fetch_results(tok_l, j),
                error = function(e) list(success = FALSE, error = conditionMessage(e))
              )
              if (fetch_res$success) {
                rv$job_results[[j$job_id]] <- fetch_res
                auto_import_results(fetch_res, j)
                rv$jobs[[i]]$status <- "completed"
                alog("OK", paste0("BV-BRC resultados importados: ", j$label),
                     paste0("Archivos descargados: ", length(fetch_res$files),
                            "\nDirectorio: ", fetch_res$out_dir))
              } else {
                alog("WARN", paste0("BV-BRC resultados no disponibles aún: ", j$label),
                     fetch_res$error %||% "Error desconocido")
              }
            }
          } else if (j$platform == "galaxy") {
            key_l <- isolate(rv$galaxy_key)
            if (!is.null(key_l)) {
              fetch_res <- tryCatch(
                galaxy_fetch_results(key_l, j),
                error = function(e) list(success = FALSE, error = conditionMessage(e))
              )
              if (fetch_res$success) {
                rv$job_results[[j$job_id]] <- fetch_res
                auto_import_results(fetch_res, j)
                rv$jobs[[i]]$status <- "completed"
                alog("OK", paste0("Galaxy resultados importados: ", j$label),
                     paste0("Datasets: ", length(fetch_res$output_ids)))
              } else {
                alog("WARN", paste0("Galaxy resultados no disponibles: ", j$label),
                     fetch_res$error %||% "Sin outputs")
              }
            }
          }
          if (!j$platform %in% c("bvbrc","galaxy")) rv$jobs[[i]]$status <- "completed"
        }, error = function(e) NULL)
      }
    }
    # Persistir estado al final de cada ciclo de sondeo
    tryCatch(state_save(rv), error = function(e) NULL)
  })

  # ── Auto-login al iniciar sesión ──────────────────────────
  # observe() provee el contexto reactivo necesario para rv$*.
  # isolate() dentro evita crear dependencias reactivas: solo corre una vez.
  observe({
    isolate({
      if (nchar(CRED_BVBRC_USER) > 0 && nchar(CRED_BVBRC_PASS) > 0) {
        alog("API", paste0("BV-BRC: auto-login como ", CRED_BVBRC_USER))
        r <- tryCatch(
          bvbrc_auth(CRED_BVBRC_USER, CRED_BVBRC_PASS),
          error = function(e) list(success = FALSE, error = conditionMessage(e))
        )
        if (r$success) {
          rv$bvbrc_token <- r$token
          rv$bvbrc_user  <- r$username
          alog("OK", paste0("BV-BRC: sesión iniciada como ", r$username),
               paste0("URL: ", BVBRC_AUTH_URL))
        } else {
          alog("ERROR", "BV-BRC: auto-login fallido", r$error)
        }
      }
      if (nchar(CRED_GALAXY_KEY) > 0) {
        alog("API", "Galaxy: auto-conexión con API key guardada")
        r <- tryCatch(
          galaxy_whoami(CRED_GALAXY_KEY),
          error = function(e) list(success = FALSE, error = conditionMessage(e))
        )
        if (r$success) {
          rv$galaxy_key <- CRED_GALAXY_KEY
          h <- tryCatch(
            galaxy_get_or_create_history(CRED_GALAXY_KEY),
            error = function(e) list(success = FALSE)
          )
          if (isTRUE(h$success)) rv$galaxy_hist <- h$history_id
          alog("OK", paste0("Galaxy: auto-conectado como ", r$username),
               if (isTRUE(h$success)) paste0("history_id: ", h$history_id) else "sin historial")
        } else {
          alog("ERROR", "Galaxy: auto-conexión fallida", r$error)
        }
      }

      # ── Cargar FASTQ desde carpeta fija de datos ────────────
      if (nchar(DATOS_DIR) > 0 && dir.exists(DATOS_DIR)) {
        # Solo archivos FASTQ completos (excluye .part, .tmp, etc.)
        all_fastq <- list.files(DATOS_DIR,
          pattern     = "\\.f(ast)?q(\\.gz)?$",
          ignore.case = TRUE, full.names = FALSE)
        part_files <- list.files(DATOS_DIR, pattern = "\\.part$", ignore.case = TRUE)

        if (length(all_fastq) > 0) {
          rv$fastq_dir   <- DATOS_DIR
          rv$fastq_files <- all_fastq

          # Extraer nombres de muestra: quitar _L00X, _R1/_R2, extensión
          nombres <- sub("_R[12].*$", "",
                      sub("_L00[0-9].*$", "", all_fastq, ignore.case = TRUE),
                      ignore.case = TRUE)
          nombres <- sub("\\.f(ast)?q(\\.gz)?$", "", nombres, ignore.case = TRUE)
          nombres <- sort(unique(nombres[nchar(nombres) > 0]))
          rv$muestras <- nombres

          alog("INFO",
            paste0("FASTQ auto-cargados desde: ", DATOS_DIR),
            paste0(length(all_fastq), " archivo(s) → ",
                   length(nombres), " muestra(s): ",
                   paste(nombres, collapse = ", "),
                   if (length(part_files) > 0)
                     paste0("\nARCHIVOS INCOMPLETOS (.part): ",
                            paste(part_files, collapse = ", "))
                   else ""))

          if (length(part_files) > 0) {
            showNotification(
              tagList(icon("exclamation-triangle"),
                paste0(" Archivos incompletos detectados (.part): ",
                       paste(part_files, collapse = ", "),
                       " — vuelve a descargarlos.")),
              type = "warning", duration = 12)
          }
        }
      }

      # ── Verificar sesión guardada en disco ──────────────────
      saved <- state_load()
      if (!is.null(saved)) {
        rv$pending_restore <- saved
        s <- state_summary(saved)
        alog("INFO", "Sesión guardada encontrada en disco",
             paste0("Guardada: ", format(saved$saved_at, "%d/%m/%Y %H:%M"),
                    " | Jobs: ", s$n_jobs,
                    " (", s$n_active, " activos)"))
        showModal(restore_modal(saved))
      } else {
        # Sin sesión guardada: traer del workspace lo que ya esté listo
        auto_sync_import()
      }
    })
  })

  # ── Auto-sincronizar e importar resultados listos ─────────
  # Registra en rv$jobs los resultados válidos del workspace y descarga
  # los que aún no tengan datos. Así, al reconectarse, la app muestra
  # el avance real sin que el usuario tenga que hacer nada.
  auto_sync_import <- function() {
    tok <- isolate(rv$bvbrc_token); user <- isolate(rv$bvbrc_user)
    if (is.null(tok) || is.null(user)) return(invisible(NULL))
    n_new <- tryCatch(bvbrc_sync_registry(tok, user), error = function(e) 0L)

    pend <- Filter(function(j)
      isTRUE(j$platform == "bvbrc") &&
      j$status %in% c("completed","FINISHED","ok","done") &&
      is.null(isolate(rv$job_results)[[j$job_id]]),
      isolate(rv$jobs))
    if (!length(pend)) {
      if (isTRUE(n_new > 0))
        alog("INFO", "[AUTO] Resultados registrados; sin descargas pendientes")
      return(invisible(NULL))
    }

    alog("INFO", paste0("[AUTO] Importando ", length(pend),
                        " resultado(s) listos del workspace..."))
    n_ok <- 0L
    withProgress(message = "Importando resultados disponibles...", value = 0, {
      for (k in seq_along(pend)) {
        j <- pend[[k]]
        incProgress(1 / length(pend), detail = j$output_file %||% j$label)
        fr <- tryCatch(bvbrc_fetch_results(tok, j),
                       error = function(e) list(success = FALSE, error = conditionMessage(e)))
        if (isTRUE(fr$success) &&
            (length(fr$contents %||% list()) > 0)) {
          rv$job_results[[j$job_id]] <- fr
          auto_import_results(fr, j)
          n_ok <- n_ok + 1L
          alog("OK", paste0("[AUTO] Importado: ", j$output_file %||% j$label),
               paste0("tipo: ", fr$parsed$type %||% "?",
                      " | archivos: ", length(fr$contents %||% list())))
        } else {
          alog("WARN", paste0("[AUTO] Sin datos legibles: ", j$output_file %||% j$label),
               fr$error %||% "el job pudo fallar en BV-BRC")
        }
      }
    })
    tryCatch(state_save(rv), error = function(e) NULL)
    if (n_ok > 0)
      showNotification(
        tagList(icon("check-circle"),
                paste0(" ", n_ok, " resultado(s) importados desde BV-BRC.")),
        type = "message", duration = 8)
    invisible(NULL)
  }

  # ─── Modal de restauración de sesión ─────────────────────
  restore_modal <- function(st) {
    s   <- state_summary(st)
    age <- as.integer(difftime(Sys.time(), st$saved_at, units = "mins"))
    age_txt <- if (age < 60) paste0(age, " min atrás")
               else          paste0(round(age / 60, 1), " h atrás")

    modalDialog(
      title     = tagList(icon("history"), " Sesión guardada encontrada"),
      size      = "m",
      easyClose = FALSE,

      div(class = "alert alert-info", style = "padding:8px 12px;",
        icon("clock"),
        strong(paste0(" Guardada el ",
                      format(st$saved_at, "%d/%m/%Y a las %H:%M"),
                      " (", age_txt, ")")),
        if (!is.null(st$bvbrc_user) && nchar(st$bvbrc_user) > 0)
          tags$small(style = "display:block;color:#555;margin-top:2px;",
                     icon("user"), " Usuario BV-BRC: ", st$bvbrc_user)
      ),

      if (s$n_jobs > 0) tagList(
        tags$b("Jobs registrados:"),
        tags$ul(style = "margin:6px 0 4px;",
          tags$li(paste0(s$n_jobs, " jobs en total")),
          if (s$n_active > 0)
            tags$li(icon("spinner", class = "fa-spin", style = "color:#e67e22;"),
              paste0(" ", s$n_active,
                     " activos — se retomarán automáticamente al restaurar")),
          if (s$n_done > 0)
            tags$li(icon("check-circle", style = "color:#27ae60;"),
              paste0(" ", s$n_done, " completados"))
        )
      ),

      if (s$has_calidad || s$has_taxonomia || s$has_ensamblado || s$n_muestras > 0)
        tagList(
          tags$b("Resultados y muestras guardados:"),
          tags$ul(style = "margin:6px 0 4px;",
            if (s$n_muestras > 0)
              tags$li(icon("dna"), paste0(" ", s$n_muestras, " muestra(s) detectadas")),
            if (s$has_calidad)
              tags$li(icon("check", style = "color:#27ae60;"), " Tabla de Calidad"),
            if (s$has_taxonomia)
              tags$li(icon("check", style = "color:#27ae60;"), " Tabla de Taxonomía"),
            if (s$has_ensamblado)
              tags$li(icon("check", style = "color:#27ae60;"), " Tabla de Ensamblado")
          )
        ),

      div(class = "alert alert-warning",
          style = "padding:6px 10px;font-size:12px;margin-top:8px;",
        icon("exclamation-triangle"),
        " Los archivos FASTQ deben volver a cargarse si la app se inició en una nueva ubicación."
      ),

      footer = tagList(
        actionButton("btn_restore_session",
          tagList(icon("history"), " Restaurar sesión"),
          class = "btn-success btn-lg"),
        actionButton("btn_new_session",
          tagList(icon("plus-circle"), " Nueva sesión"),
          class = "btn-default")
      )
    )
  }

  # ─── Modal de credenciales ─────────────────────────────────
  creds_modal <- function() {
    modalDialog(
      title = tagList(icon("key"), " Credenciales de plataformas"),
      size = "m", easyClose = TRUE,
      tabsetPanel(
        tabPanel("BV-BRC", br(),
          div(class="alert alert-info", icon("info-circle"),
              " Cuenta gratuita en ", tags$a("www.bv-brc.org", href="https://www.bv-brc.org", target="_blank"),
              ". Las credenciales NO se guardan en disco."),
          textInput("cred_bvbrc_user", "Usuario (email):",
                    value = if (nchar(CRED_BVBRC_USER) > 0) CRED_BVBRC_USER else "",
                    placeholder = "usuario@correo.com"),
          passwordInput("cred_bvbrc_pass", "Contraseña:"),
          uiOutput("bvbrc_auth_st"), br(),
          actionButton("btn_bvbrc_login", tagList(icon("sign-in-alt")," Iniciar sesión BV-BRC"),
                       class="btn-primary", style="width:100%;")
        ),
        tabPanel("Galaxy", br(),
          div(class="alert alert-info", icon("info-circle"),
              " API key: usegalaxy.org → User → Preferences → Manage API Key. ",
              tags$a("Abrir Galaxy", href="https://usegalaxy.org/user/api_key", target="_blank")),
          textInput("cred_galaxy_key", "API Key:",
                    value = if (nchar(CRED_GALAXY_KEY) > 0) CRED_GALAXY_KEY else "",
                    placeholder = "Pegar aquí..."),
          uiOutput("galaxy_auth_st"), br(),
          actionButton("btn_galaxy_connect", tagList(icon("plug")," Conectar Galaxy"),
                       class="btn-primary", style="width:100%;")
        )
      ),
      footer = modalButton("Cerrar")
    )
  }

  # ─── Restaurar sesión ────────────────────────────────────
  observeEvent(input$btn_restore_session, {
    st <- rv$pending_restore
    if (is.null(st)) { removeModal(); return() }

    rv$jobs         <- st$jobs         %||% list()
    rv$job_results  <- st$job_results  %||% list()
    rv$calidad      <- st$calidad      %||% rv$calidad
    rv$taxonomia    <- st$taxonomia    %||% rv$taxonomia
    rv$ensamblado   <- st$ensamblado   %||% rv$ensamblado
    rv$mlst         <- st$mlst         %||% rv$mlst
    rv$alineamiento <- st$alineamiento %||% rv$alineamiento
    rv$ani          <- st$ani          %||% rv$ani
    rv$resfinder_df <- st$resfinder_df
    rv$ebi_result   <- st$ebi_result
    rv$muestras     <- st$muestras     %||% rv$muestras
    rv$fastq_dir    <- st$fastq_dir
    rv$fastq_files  <- st$fastq_files  %||% rv$fastq_files
    rv$completados  <- st$completados  %||% rv$completados

    s <- state_summary(st)
    rv$pending_restore <- NULL
    removeModal()

    alog("OK", "Sesión restaurada desde disco",
         paste0("Jobs recuperados: ", s$n_jobs,
                " | Activos (retomarán sondeo): ", s$n_active,
                " | Muestras: ", s$n_muestras))
    showNotification(
      tagList(icon("history"),
        paste0(" Sesión restaurada — ", s$n_jobs, " jobs, ",
               s$n_active, " retomarán monitoreo automáticamente.")),
      type = "message", duration = 8)
    # Traer del workspace lo que haya terminado mientras la app estuvo cerrada
    auto_sync_import()
  })

  observeEvent(input$btn_new_session, {
    state_clear()
    rv$pending_restore <- NULL
    removeModal()
    alog("INFO", "Nueva sesión iniciada — estado anterior descartado")
    # Sincronizar registro con el workspace: aunque sea sesión nueva, detecta
    # resultados ya existentes en BV-BRC para no re-subir ni re-analizar.
    tok <- isolate(rv$bvbrc_token); user <- isolate(rv$bvbrc_user)
    n_sync <- tryCatch(bvbrc_sync_registry(tok, user), error = function(e) 0L)
    if (isTRUE(n_sync > 0)) {
      # Marcar como completadas las etapas cuyos resultados ya existen
      for (j in rv$jobs) {
        pid <- TOOL_TO_STEP[j$tool] %||% NULL
        if (!is.null(pid) && pid %in% names(rv$completados)) rv$completados[pid] <- TRUE
      }
      showNotification(
        paste0("Nueva sesión. Se detectaron ", n_sync,
               " resultado(s) ya existentes en BV-BRC — no se re-analizarán."),
        type = "message", duration = 8)
    } else {
      showNotification("Nueva sesión iniciada.", type="message", duration=4)
    }
    # Importar los resultados que ya estén listos en el workspace
    auto_sync_import()
  })

  # ── Botón manual: sincronizar e importar resultados ────────
  observeEvent(input$btn_sync_now, {
    if (is.null(rv$bvbrc_token)) {
      showModal(creds_modal())
      showNotification("Inicia sesión BV-BRC primero.", type = "warning", duration = 5)
      return()
    }
    auto_sync_import()
  }, ignoreInit = TRUE)

  observeEvent(input$btn_credentials, showModal(creds_modal()))
  observeEvent(input$btn_open_creds,  showModal(creds_modal()))

  observeEvent(input$btn_bvbrc_login, {
    req(nchar(input$cred_bvbrc_user)>0, nchar(input$cred_bvbrc_pass)>0)
    alog("API", paste0("BV-BRC: intentando login como ", input$cred_bvbrc_user))
    r <- bvbrc_auth(input$cred_bvbrc_user, input$cred_bvbrc_pass)
    if (r$success) {
      rv$bvbrc_token <- r$token; rv$bvbrc_user <- r$username
      alog("OK", paste0("BV-BRC: sesion iniciada como ", r$username),
           paste0("URL: ", BVBRC_AUTH_URL))
      showNotification(paste("BV-BRC: sesión como", r$username), type="message", duration=5)
    } else {
      alog("ERROR", "BV-BRC: autenticacion fallida", r$error)
      showNotification(paste("BV-BRC:", r$error), type="error", duration=8)
    }
  })

  observeEvent(input$btn_galaxy_connect, {
    req(nchar(input$cred_galaxy_key)>5)
    alog("API", "Galaxy: verificando API key...")
    r <- galaxy_whoami(input$cred_galaxy_key)
    if (r$success) {
      rv$galaxy_key <- input$cred_galaxy_key
      h <- galaxy_get_or_create_history(input$cred_galaxy_key)
      if (h$success) {
        rv$galaxy_hist <- h$history_id
        alog("OK", paste0("Galaxy: conectado como ", r$username),
             paste0("history_id: ", h$history_id))
      } else {
        alog("WARN", paste0("Galaxy: autenticado como ", r$username, " pero no se pudo obtener historial"), h$error)
      }
      showNotification(paste("Galaxy: conectado como", r$username), type="message", duration=5)
    } else {
      alog("ERROR", "Galaxy: error de conexion", r$error)
      showNotification(paste("Galaxy:", r$error), type="error", duration=8)
    }
  })

  output$bvbrc_auth_st  <- renderUI({ if (!is.null(rv$bvbrc_token)) div(class="alert alert-success",style="padding:6px 10px;margin:6px 0;", icon("check-circle"), paste(" Sesión:", rv$bvbrc_user)) else div(class="alert alert-warning",style="padding:6px 10px;margin:6px 0;", icon("exclamation-triangle"), " No conectado") })
  output$galaxy_auth_st <- renderUI({ if (!is.null(rv$galaxy_key))  div(class="alert alert-success",style="padding:6px 10px;margin:6px 0;", icon("check-circle"), " Galaxy conectado") else div(class="alert alert-warning",style="padding:6px 10px;margin:6px 0;", icon("exclamation-triangle"), " Sin API key") })
  output$cred_status_bvbrc  <- renderUI({ if (!is.null(rv$bvbrc_token)) div(class="alert alert-success",style="padding:8px;", icon("check-circle"), strong(" BV-BRC"), br(), tags$small(rv$bvbrc_user)) else div(class="alert alert-danger",style="padding:8px;", icon("times-circle"), strong(" BV-BRC"), br(), tags$small("Sin sesión")) })
  output$cred_status_galaxy <- renderUI({ if (!is.null(rv$galaxy_key))  div(class="alert alert-success",style="padding:8px;", icon("check-circle"), strong(" Galaxy"), br(), tags$small("API key activa")) else div(class="alert alert-danger",style="padding:8px;", icon("times-circle"), strong(" Galaxy"), br(), tags$small("Sin API key")) })

  # ─── Archivos FASTQ ────────────────────────────────────────
  volumes <- c("Escritorio 2026"="C:/Escritorio 2026", "Documentos"=path.expand("~/Documents"), "C:"="C:/", "D:"="D:/")
  volumes <- volumes[dir.exists(volumes)]
  shinyDirChoose(input, "carpeta_fastq", roots=volumes, session=session, restrictions=system.file(package="base"))

  observeEvent(input$fastq_upload, {
    files <- input$fastq_upload
    if (is.null(files)) return()
    staged <- list()
    for (i in seq_len(nrow(files))) {
      fname <- files$name[i]; src <- files$datapath[i]; dest <- file.path(SESSION_TMPDIR, fname)
      file.copy(src, dest, overwrite=TRUE); staged[[fname]] <- dest
    }
    rv$staged_files <- c(rv$staged_files, staged)
    nombres <- sort(unique(sub("_R[12].*$","",sub("_L00[0-9].*$","",names(staged),ignore.case=TRUE),ignore.case=TRUE)))
    nombres <- sub("\\.f(ast)?q(\\.gz)?$","",nombres,ignore.case=TRUE)
    rv$muestras <- sort(unique(c(rv$muestras, nombres[nchar(nombres)>0])))
    rv$fastq_files <- unique(c(rv$fastq_files, files$name))
    alog("INFO",
         paste0(nrow(files), " archivo(s) cargado(s) → ", length(rv$muestras), " muestra(s)"),
         paste(files$name, collapse="\n"))
    showNotification(paste(nrow(files),"archivo(s) subido(s)."), type="message", duration=4)
  })

  carpeta_path <- reactive({
    req(input$carpeta_fastq); if (is.integer(input$carpeta_fastq)) return(NULL)
    parseDirPath(volumes, input$carpeta_fastq)
  })
  observeEvent(carpeta_path(), {
    ruta <- carpeta_path(); if (is.null(ruta)||ruta=="") return()
    archivos <- list.files(ruta, pattern="\\.f(ast)?q(\\.gz)?$", ignore.case=TRUE)
    rv$fastq_dir   <- ruta
    rv$fastq_files <- unique(c(rv$fastq_files, archivos))
    # Misma logica que fileInput: quitar _L00X, _R1/_R2 y extension
    nombres <- sort(unique(
      sub("\\.f(ast)?q(\\.gz)?$", "",
        sub("_R[12].*$", "",
          sub("_L00[0-9]_.*$|_L00[0-9]\\..*$|_L00[0-9]$", "",
              archivos, ignore.case = TRUE),
          ignore.case = TRUE),
        ignore.case = TRUE)
    ))
    nombres <- nombres[nchar(nombres) > 0]
    rv$muestras <- sort(unique(c(rv$muestras, nombres)))
    alog("INFO", paste0(length(archivos), " FASTQ en carpeta → ", length(nombres), " muestra(s)"),
         paste(archivos, collapse="\n"))
    showNotification(paste(length(archivos),"archivos FASTQ en carpeta."), type="message", duration=4)
  })

  output$carpeta_seleccionada <- renderUI({
    ruta <- carpeta_path(); if (is.null(ruta)||ruta=="") return(NULL)
    div(style="background:#eaf4fb;border:1px solid #aed6f1;border-radius:4px;padding:7px 10px;margin-top:8px;font-size:12px;", icon("folder",style="color:#2e86c1;"), " ", tags$code(ruta))
  })

  output$muestras_detectadas <- renderUI({
    if (length(rv$muestras)==0) return(div(style="color:#888;font-size:13px;", icon("arrow-circle-up"), " Sube archivos FASTQ para detectar muestras."))
    tagList(div(style="margin-bottom:6px;", tags$b(icon("check-circle",style="color:#27ae60;"), paste0(" ", length(rv$fastq_files)," archivo(s) — ", length(rv$muestras)," muestra(s):"))), div(lapply(rv$muestras, function(m) span(class="muestra-chip", icon("dna"), " ", m))))
  })

  output$bs_status <- renderUI({
    req(input$bs_cpf, input$bs_q30)
    ok <- !is.na(input$bs_cpf)&&!is.na(input$bs_q30)&&input$bs_cpf>=80&&input$bs_q30>=75
    if (ok) div(class="alert alert-success",style="margin:0;",icon("check-circle")," Corrida: APROBADA") else div(class="alert alert-danger",style="margin:0;",icon("times-circle")," Corrida: REVISAR")
  })

  output$inicio_validacion <- renderUI({
    if (length(rv$muestras)==0) div(class="alert alert-warning",style="display:inline-block;", icon("exclamation-triangle"), " Sube archivos FASTQ antes de continuar.")
    else div(class="alert alert-success",style="display:inline-block;", icon("check-circle"), paste0(" ", length(rv$muestras)," muestra(s): ", paste(rv$muestras,collapse=", ")))
  })

  observeEvent(input$btn_iniciar_analisis, {
    if (length(rv$muestras)==0) showNotification("Sube archivos FASTQ primero.", type="warning", duration=4)
    else updateTabItems(session, "sidebar", selected="calidad")
  }, ignoreInit=TRUE)

  # ─── Helpers de ejecución ──────────────────────────────────
  get_sample_files <- function(muestra) {
    # ── 1. Buscar en archivos staged via fileInput ─────────────
    all_f  <- rv$staged_files
    r1_key <- names(all_f)[grepl(paste0(muestra, ".*_R1"), names(all_f), ignore.case = TRUE)]
    r2_key <- names(all_f)[grepl(paste0(muestra, ".*_R2"), names(all_f), ignore.case = TRUE)]
    if (length(r1_key) == 0)
      r1_key <- names(all_f)[grepl(muestra, names(all_f), ignore.case = TRUE)][1]
    if (length(r2_key) == 0) r2_key <- r1_key

    r1_path <- if (length(r1_key) > 0 && !is.na(r1_key[1])) all_f[[r1_key[1]]] else NULL
    # R2: solo si se encontró una clave distinta de R1
    r2_path <- if (length(r2_key) > 0 && !is.na(r2_key[1]) &&
                   !identical(r2_key[1], r1_key[1]))
                 all_f[[r2_key[1]]] else NULL

    # ── 2. Fallback: carpeta seleccionada con el explorador ────
    if (is.null(r1_path) && !is.null(rv$fastq_dir) && dir.exists(rv$fastq_dir)) {
      dir_full <- list.files(rv$fastq_dir,
                             pattern     = "\\.f(ast)?q(\\.gz)?$",
                             ignore.case = TRUE, full.names = TRUE)
      bnm <- basename(dir_full)

      i1 <- which(grepl(paste0(muestra, ".*_R1"), bnm, ignore.case = TRUE))
      i2 <- which(grepl(paste0(muestra, ".*_R2"), bnm, ignore.case = TRUE))
      # sin _R1/_R2 explícito: cualquier archivo que contenga la muestra va a R1
      if (length(i1) == 0) i1 <- which(grepl(muestra, bnm, ignore.case = TRUE))
      # R2 se deja NULL si no existe archivo distinto al R1

      r1_path <- if (length(i1) > 0) dir_full[i1[1]] else NULL
      r2_path <- if (length(i2) > 0) dir_full[i2[1]] else NULL

      alog("DEBUG", paste0("Archivos en carpeta para '", muestra, "'"),
           paste0("R1: ", r1_path %||% "ninguno", "\nR2: ", r2_path %||% "ninguno"))
    }

    list(r1 = r1_path, r2 = r2_path)
  }

  add_job <- function(job_id, platform, tool, step, label,
                       output_path = NULL, output_file = NULL,
                       muestra = NULL, output_ids = NULL) {
    rv$jobs <- c(rv$jobs, list(list(
      job_id = as.character(job_id), platform = platform, tool = tool,
      step = step, label = label, status = "queued",
      time = Sys.time(), result_file = NULL, result_text = NULL,
      output_path = output_path, output_file = output_file,
      muestra = muestra, output_ids = output_ids
    )))
    # Persistir inmediatamente para que el job_id sobreviva reinicios
    tryCatch(state_save(rv), error = function(e) NULL)
  }

  # ─── Helpers: detectar jobs BV-BRC ya completados ─────────
  bvbrc_expected_out <- function(tool_id, muestra) {
    switch(tool_id,
      FastqUtils              = paste0("fastqc_",     muestra),
      TaxonomicClassification = paste0("taxonomy_",   muestra),
      GenomeDistance          = paste0("taxonomy_",   muestra),
      Assembly2               = paste0("assembly_",   muestra),
      MetagenomicReadMapping  = paste0("resistoma_",  muestra),
      Annotation              = paste0("annotation_", muestra),
      ComprehensiveGenomeAnalysis = paste0("cga_",    muestra),
      CoreGenomeMLST          = paste0("mlst",        muestra),
      PhylogeneticTree        = paste0("filogenia",   muestra),
      NULL
    )
  }

  bvbrc_get_done_job <- function(jobs, tool_id, muestra) {
    exp_out <- bvbrc_expected_out(tool_id, muestra)
    if (is.null(exp_out)) return(NULL)
    matches <- Filter(function(j)
      isTRUE(j$platform == "bvbrc") &&
      isTRUE(j$tool == tool_id) &&
      isTRUE((j$output_file %||% "") == exp_out) &&
      j$status %in% c("FINISHED", "completed", "ok", "done", "C"),
      jobs)
    if (length(matches)) matches[[1]] else NULL
  }

  # ── Sincronizar registro desde el workspace ────────────────
  # Lista ws_out(user) (PRUEBAAPP) y registra CADA carpeta de resultado
  # existente (fastqc_*, taxonomy_*, assembly_*, …) como un job completado,
  # AUNQUE su import haya fallado. Así el dedup evita re-subir/re-analizar
  # muestras que ya se procesaron, incluso en una sesión nueva.
  # Devuelve el nº de entradas nuevas registradas.
  bvbrc_sync_registry <- function(tok, user, verbose = TRUE) {
    if (is.null(tok) || is.null(user)) return(invisible(0L))
    ls <- tryCatch(bvbrc_ls(tok, ws_out(user)), error = function(e) list(success = FALSE))
    if (!isTRUE(ls$success)) {
      if (verbose)
        alog("WARN", "[SYNC] No se pudo listar el workspace para sincronizar",
             ls$error %||% "Workspace.ls falló")
      return(invisible(0L))
    }
    n_added <- 0L; n_skip_fail <- 0L
    for (f in ls$files %||% list()) {
      # Solo objetos job_result (no las carpetas ocultas .name ni folders sueltos)
      if (!identical(f$type %||% "", "job_result")) next
      tool <- ws_folder_tool(f$name %||% "")
      if (is.null(tool)) next
      prefix  <- WS_FOLDER_MAP[[tool]] %||% "^[^_]+_"
      muestra <- sub(prefix, "", f$name)
      job_id  <- paste0("ws_", f$name)
      ya <- any(vapply(rv$jobs, function(j)
        identical(j$job_id, job_id) ||
        (isTRUE(j$tool == tool) && isTRUE((j$output_file %||% "") == f$name)),
        logical(1)))
      if (ya) next
      # No registrar como completado un job fallido o sin salida → re-ejecutable
      valido <- bvbrc_result_ok(tok, f$path %||% paste0(ws_out(user), "/", f$name))
      if (identical(valido, FALSE)) {
        n_skip_fail <- n_skip_fail + 1L
        next
      }
      rv$jobs <- c(rv$jobs, list(list(
        job_id = job_id, platform = "bvbrc", tool = tool,
        step = TOOL_TO_STEP[tool] %||% "inicio",
        label = paste0(tool, " — ", muestra), status = "completed",
        time = Sys.time(),
        output_path = ws_out(user), output_file = f$name,
        muestra = muestra, output_ids = NULL
      )))
      n_added <- n_added + 1L
    }
    if (n_added > 0) {
      tryCatch(state_save(rv), error = function(e) NULL)
      if (verbose)
        alog("OK", paste0("[SYNC] ", n_added,
                          " resultado(s) válidos registrados desde el workspace"),
             paste0("Carpeta: ", ws_out(user),
                    " — el dedup evitará re-analizar estas muestras.",
                    if (n_skip_fail > 0)
                      paste0("\n", n_skip_fail,
                             " carpeta(s) con job fallido/sin salida se ignoraron ",
                             "(se pueden re-ejecutar).") else ""))
    } else if (n_skip_fail > 0 && verbose) {
      alog("INFO", paste0("[SYNC] ", n_skip_fail,
                          " resultado(s) fallidos/vacíos en el workspace — re-ejecutables"))
    }
    invisible(n_added)
  }

  # ── Refrescar TODAS las etapas desde los servidores ────────
  # Se llama al reanudar una sesión: sincroniza el registro con el
  # workspace, consulta el estado de los jobs pendientes e importa
  # los resultados que ya estén disponibles. Devuelve un resumen.
  refresh_all_stages <- function(verbose = TRUE) {
    tok  <- isolate(rv$bvbrc_token); user <- isolate(rv$bvbrc_user)
    key  <- isolate(rv$galaxy_key)
    n_new <- 0L; n_imp <- 0L; n_err <- 0L

    # 1) Descubrir resultados existentes en el workspace de BV-BRC
    if (!is.null(tok) && !is.null(user))
      n_new <- tryCatch(bvbrc_sync_registry(tok, user, verbose = verbose),
                        error = function(e) 0L)

    # 2) Actualizar estado e importar resultados de cada job pendiente
    jobs <- isolate(rv$jobs)
    for (i in seq_along(jobs)) {
      j <- jobs[[i]]
      if (!is.null(rv$job_results[[j$job_id]])) next   # ya importado

      # Actualizar estado si sigue en marcha
      if (j$status %in% c("queued","running","Q","R","submitted","in_progress")) {
        st <- tryCatch(switch(j$platform,
          bvbrc  = if (!is.null(tok)) bvbrc_job_status(tok, j$job_id)$status else NULL,
          galaxy = if (!is.null(key)) galaxy_job_status(key, j$job_id)$status else NULL,
          NULL), error = function(e) NULL)
        if (!is.null(st)) rv$jobs[[i]]$status <- st
      }

      # Intentar importar (el fetch decide si hay algo que traer)
      fr <- tryCatch(switch(j$platform,
        bvbrc  = if (!is.null(tok)) bvbrc_fetch_results(tok, j) else NULL,
        galaxy = if (!is.null(key)) galaxy_fetch_results(key, j) else NULL,
        NULL), error = function(e) list(success = FALSE, error = conditionMessage(e)))

      if (isTRUE(fr$success)) {
        rv$job_results[[j$job_id]] <- fr
        auto_import_results(fr, j)
        rv$jobs[[i]]$status <- "completed"
        n_imp <- n_imp + 1L
      } else if (!is.null(fr)) n_err <- n_err + 1L
    }

    if (verbose)
      alog("OK", "[SESIÓN] Etapas actualizadas desde los servidores",
           paste0("Nuevos resultados detectados: ", n_new,
                  " | Importados ahora: ", n_imp,
                  " | Aún sin resultado: ", n_err))
    tryCatch(state_save(rv), error = function(e) NULL)
    list(nuevos = n_new, importados = n_imp, pendientes = n_err)
  }

  # ─── Ejecutar job por herramienta ─────────────────────────
  run_exec <- function(step_id, tool_spec, muestra) {
    if (is.null(tool_spec)||is.null(muestra)||nchar(muestra)==0) {
      showNotification("Selecciona muestra y herramienta.", type="warning"); return()
    }
    parts <- strsplit(tool_spec,":")[[1]]; platform <- parts[1]; tool_id <- parts[2]
    sfiles <- get_sample_files(muestra)
    alog("API",
         paste0("[", toupper(platform), "] Ejecutando: ", tool_id, " | muestra: ", muestra, " | etapa: ", step_id),
         paste0("Archivos encontrados — R1: ", sfiles$r1 %||% "ninguno", " | R2: ", sfiles$r2 %||% "ninguno"))

    # EBI JDispatcher ─────────────────────────────────────────
    if (platform == "ebi") {
      ffile <- sfiles$r1
      if (is.null(ffile)||!file.exists(ffile)) {
        showNotification("Archivo no encontrado para esta muestra.", type="warning"); return()
      }
      raw <- readLines(ffile, n=100, warn=FALSE)
      fasta <- if (any(grepl("^>",raw))) paste(head(raw,20),collapse="\n") else paste0(">",muestra,"\n",paste(head(raw,5),collapse=""))
      res <- switch(tool_id, muscle=ebi_muscle(fasta), clustalo=ebi_clustalo(fasta), mafft=ebi_mafft(fasta), list(success=FALSE,error="Herramienta desconocida"))
      if (res$success) {
        add_job(res$job_id,"ebi",tool_id,step_id,paste0(toupper(tool_id)," — ",muestra))
        alog("OK", paste0("EBI job enviado: ", toupper(tool_id)), paste0("job_id: ", res$job_id, "\nmuestra: ", muestra))
        showNotification(paste("EBI job enviado. ID:",res$job_id,"(verificar en ~2 min)"), type="message", duration=6)
      } else {
        alog("ERROR", paste0("EBI error — ", tool_id), res$error)
        showNotification(paste("Error EBI:",res$error), type="error", duration=8)
      }

    # NCBI BLAST ──────────────────────────────────────────────
    } else if (platform == "ncbi") {
      ffile <- sfiles$r1
      if (is.null(ffile)||!file.exists(ffile)) { showNotification("Archivo no encontrado.", type="warning"); return() }
      raw <- readLines(ffile, n=50, warn=FALSE)
      fasta <- paste(head(raw,10),collapse="\n")
      res <- ncbi_blast_submit(fasta, program=tool_id)
      if (res$success) {
        add_job(res$rid,"ncbi",tool_id,step_id,paste0("BLAST — ",muestra))
        alog("OK", paste0("NCBI BLAST enviado (", tool_id, ")"), paste0("RID: ", res$rid, "\nTiempo estimado: ", res$estimated_seconds, " s"))
        showNotification(paste("BLAST enviado. RID:",res$rid,"— Est.:",res$estimated_seconds,"s"), type="message", duration=7)
      } else {
        alog("ERROR", "NCBI BLAST error", res$error)
        showNotification(paste("Error NCBI:",res$error), type="error", duration=8)
      }

    # ResFinder / CGE ─────────────────────────────────────────
    } else if (platform == "cge") {
      ffile <- sfiles$r1
      if (is.null(ffile)||!file.exists(ffile)) { showNotification("Archivo no encontrado.", type="warning"); return() }
      res <- resfinder_submit(fasta_file=ffile, species="other")
      if (res$success) {
        if (!is.null(res$job_id)) {
          add_job(res$job_id,"cge","resfinder",step_id,paste0("ResFinder — ",muestra))
          alog("OK", "ResFinder job enviado (CGE/DTU)", paste0("job_id: ", res$job_id))
          showNotification(paste("ResFinder enviado. ID:",res$job_id), type="message", duration=5)
        } else if (!is.null(res$result_direct)) {
          parsed <- tryCatch(jsonlite::fromJSON(res$result_direct,simplifyVector=FALSE), error=function(e)NULL)
          if (!is.null(parsed)) rv$resfinder_df <- resfinder_to_table(parsed)
          alog("OK", "ResFinder: resultado directo recibido", paste0("Genes AMR: ", if(!is.null(rv$resfinder_df)) nrow(rv$resfinder_df) else "?"))
          showNotification("ResFinder: resultados recibidos.", type="message", duration=5)
        }
      } else {
        alog("ERROR", "ResFinder error (CGE/DTU)", res$error)
        showNotification(paste("Error ResFinder:",res$error), type="error", duration=8)
      }

    # BV-BRC ──────────────────────────────────────────────────
    } else if (platform == "bvbrc") {
      tok <- rv$bvbrc_token; user <- rv$bvbrc_user
      if (is.null(tok)) { showModal(creds_modal()); showNotification("Inicia sesión BV-BRC primero.", type="warning", duration=5); return() }

      # Sincronizar registro con el workspace para no re-analizar lo ya hecho
      tryCatch(bvbrc_sync_registry(tok, user), error = function(e) NULL)

      # ── Verificar si ya existe un job completado para esta muestra ──
      ej_done <- bvbrc_get_done_job(rv$jobs, tool_id, muestra)
      if (!is.null(ej_done)) {
        if (!is.null(rv$job_results[[ej_done$job_id]])) {
          showNotification(
            paste0(tool_id, " — ", muestra, ": resultados ya importados. Ver pestaña de resultados."),
            type = "message", duration = 7)
          return()
        }
        showNotification(
          paste0(tool_id, " — ", muestra, ": job ya completado. Importando resultados existentes..."),
          type = "message", duration = 5)
        fetch_res <- tryCatch(
          bvbrc_fetch_results(tok, ej_done),
          error = function(e) list(success = FALSE, error = conditionMessage(e))
        )
        ji <- which(vapply(rv$jobs, function(j) j$job_id == ej_done$job_id, logical(1)))
        if (isTRUE(fetch_res$success)) {
          rv$job_results[[ej_done$job_id]] <- fetch_res
          auto_import_results(fetch_res, ej_done)
          if (length(ji)) rv$jobs[[ji[1]]]$status <- "completed"
          alog("OK", paste0("Re-import: ", ej_done$label),
               paste0(length(fetch_res$files), " archivo(s)"))
          showNotification(paste0("Resultados importados para ", muestra, "."),
                           type = "message", duration = 5)
        } else {
          alog("WARN", paste0("Re-import fallido: ", ej_done$label),
               fetch_res$error %||% "Error desconocido")
          showNotification(paste0("No se pudieron importar resultados: ",
                                  fetch_res$error %||% "Error desconocido"),
                           type = "warning", duration = 8)
        }
        return()
      }

      showModal(modalDialog(
        title    = tagList(icon("cloud-upload-alt"), " BV-BRC — Procesando..."),
        tags$p(icon("spinner", class = "fa-spin"), " Subiendo archivos y enviando job a BV-BRC..."),
        tags$small(style = "color:#888;", "Este proceso puede tardar 1-3 minutos. No cierre la ventana."),
        footer   = NULL,
        easyClose = FALSE
      ))

      # Crear carpeta de salida si no existe (ws_out usa /home/ANbio_output)
      tryCatch(bvbrc_mkdir(tok, ws_out(user)), error = function(e) NULL)

      r1 <- sfiles$r1; r2 <- sfiles$r2
      ws_r1 <- paste0(ws_home(user), "/.fastq_uploads/", muestra, "_R1.fastq.gz")
      ws_r2 <- NULL  # NULL = single-end por defecto hasta confirmar upload R2
      if (!is.null(r1)&&file.exists(r1)) {
        alog("API", paste0("BV-BRC upload R1: ", basename(r1)))
        up1 <- bvbrc_upload_file(tok,user,r1)
        if (!up1$success) {
          removeModal()
          alog("ERROR", "BV-BRC upload R1 fallido", up1$log_det %||% up1$error)
          showNotification(paste("Upload error:",up1$error),type="error",duration=8); return()
        }
        alog("OK", "BV-BRC upload R1 exitoso", up1$log_det)
        ws_r1 <- up1$ws_path
        if (!is.null(r2)&&file.exists(r2)&&r2!=r1) {
          alog("API", paste0("BV-BRC upload R2: ", basename(r2)))
          up2 <- bvbrc_upload_file(tok,user,r2)
          if (up2$success) {
            ws_r2 <- up2$ws_path
            alog("OK", "BV-BRC upload R2 exitoso", up2$log_det)
          } else {
            alog("WARN", "BV-BRC upload R2 fallido — modo single-end", up2$log_det %||% up2$error)
          }
        } else if (is.null(r2)) {
          alog("INFO", paste0("BV-BRC R2 no disponible — modo single-end: ", muestra))
        }
      }
      organism <- input$proy_organismo %||% "Unknown"
      params <- switch(tool_id,
        FastqUtils                  = bvbrc_params_fastqutil(user,ws_r1,ws_r2,muestra),
        TaxonomicClassification     = bvbrc_params_taxonomy(user,ws_r1,ws_r2,muestra),
        GenomeDistance              = bvbrc_params_taxonomy(user,ws_r1,ws_r2,muestra),
        Assembly2                   = bvbrc_params_assembly(user,ws_r1,ws_r2,muestra),
        ComprehensiveGenomeAnalysis = bvbrc_params_cga(user,ws_r1,ws_r2,organism,muestra),
        Annotation                  = bvbrc_params_annotation(user,ws_r1,organism,muestra),
        MetagenomicReadMapping      = bvbrc_params_metareads(user,ws_r1,ws_r2,muestra),
        CoreGenomeMLST              = bvbrc_params_mlst(user,list(muestra),"mlst"),
        PhylogeneticTree            = bvbrc_params_phylogeny(user,list(muestra),"filogenia"),
        SeqComparison               = bvbrc_params_proteome(user,muestra,muestra,"comparacion"),
        Homology                    = bvbrc_params_blast(user,muestra,"blastn",muestra),
        list()
      )
      alog("API", paste0("BV-BRC submit: AppService.start_app [", tool_id, "]"),
           paste0("JSON-RPC → ", BVBRC_APP_URL, "\nparams keys: ", paste(names(params), collapse=", ")))
      res <- bvbrc_submit_job(tok, tool_id, params, ws_home(user))
      removeModal()
      if (res$success) {
        add_job(res$task_id, "bvbrc", tool_id, step_id,
                paste0(tool_id, " — ", muestra),
                output_path = params$output_path,
                output_file = params$output_file,
                muestra     = muestra)
        alog("OK", paste0("BV-BRC job aceptado: ", tool_id), paste0("task_id: ", res$task_id, "\n", res$log_det %||% ""))
        showNotification(paste("BV-BRC job enviado. ID:",res$task_id), type="message", duration=5)
      } else {
        alog("ERROR", paste0("BV-BRC submit fallido: ", tool_id), paste0(res$error, "\n", res$log_det %||% ""))
        showNotification(paste("Error BV-BRC:",res$error), type="error", duration=8)
      }

    # Galaxy ──────────────────────────────────────────────────
    } else if (platform == "galaxy") {
      key <- rv$galaxy_key; hist <- rv$galaxy_hist
      if (is.null(key)) { showModal(creds_modal()); showNotification("Configura API key Galaxy.", type="warning", duration=5); return() }
      r1 <- sfiles$r1; r2 <- sfiles$r2
      if (is.null(r1)||!file.exists(r1)) {
        alog("ERROR", paste0("Galaxy: archivo R1 no encontrado para muestra ", muestra))
        showNotification("Archivo no encontrado.", type="warning"); return()
      }

      showModal(modalDialog(
        title    = tagList(icon("cloud-upload-alt"), " Galaxy — Subiendo archivos..."),
        tags$p(icon("spinner", class="fa-spin"), " Subiendo archivos FASTQ a Galaxy..."),
        tags$small(style="color:#888;",
          paste0("Archivo: ", basename(r1), " — Este proceso puede tardar varios minutos para archivos grandes.")),
        footer    = NULL,
        easyClose = FALSE
      ))

      alog("API", paste0("Galaxy upload R1: ", basename(r1)), paste0("history_id: ", hist))
      up1 <- galaxy_upload_file(key, hist, r1, "fastqsanger.gz")
      if (!up1$success) {
        removeModal()
        alog("ERROR", paste0("Galaxy upload R1 fallido — ", muestra), up1$error)
        showNotification(paste("Galaxy upload error:", up1$error), type="error", duration=8)
        return()
      }
      alog("OK", paste0("Galaxy upload R1 exitoso — ", basename(r1)), paste0("dataset_id: ", up1$dataset_id))

      up2 <- up1
      if (!is.null(r2) && file.exists(r2) && r2 != r1) {
        alog("API", paste0("Galaxy upload R2: ", basename(r2)))
        up2 <- galaxy_upload_file(key, hist, r2, "fastqsanger.gz")
        if (!up2$success) {
          removeModal()
          alog("ERROR", paste0("Galaxy upload R2 fallido — ", muestra), up2$error)
          showNotification(paste("Galaxy upload R2 error:", up2$error), type="error", duration=8)
          return()
        }
        alog("OK", paste0("Galaxy upload R2 exitoso — ", basename(r2)), paste0("dataset_id: ", up2$dataset_id))
      }

      # Esperar a que los datasets estén en state="ok" antes de lanzar la herramienta.
      # Sin esta espera Galaxy devuelve HTTP 400 ToolInputsNotReadyException.
      showModal(modalDialog(
        title    = tagList(icon("clock"), " Galaxy — Preparando datasets..."),
        tags$p(icon("spinner", class="fa-spin"),
               " Esperando que Galaxy procese los archivos subidos..."),
        tags$small(style="color:#888;", "Esto puede tomar 1-2 minutos."),
        footer = NULL, easyClose = FALSE
      ))
      alog("INFO", "Galaxy: esperando que datasets estén listos (state=ok)")
      w1 <- galaxy_wait_dataset(key, up1$dataset_id, hist, max_wait = 300)
      if (!w1$success) {
        removeModal()
        alog("ERROR", paste0("Galaxy dataset R1 no listo — ", muestra), w1$error)
        showNotification(paste("Galaxy dataset error:", w1$error), type="error", duration=8)
        return()
      }
      if (up2$dataset_id != up1$dataset_id) {
        w2 <- galaxy_wait_dataset(key, up2$dataset_id, hist, max_wait = 300)
        if (!w2$success) {
          removeModal()
          alog("ERROR", paste0("Galaxy dataset R2 no listo — ", muestra), w2$error)
          showNotification(paste("Galaxy dataset R2 error:", w2$error), type="error", duration=8)
          return()
        }
      }
      alog("OK", "Galaxy: datasets listos — lanzando herramienta")
      removeModal()

      inputs <- switch(tool_id,
        spades    = galaxy_inputs_spades(up1$dataset_id, up2$dataset_id),
        unicycler = galaxy_inputs_unicycler(up1$dataset_id, up2$dataset_id),
        mlst      = galaxy_inputs_mlst(up1$dataset_id),
        iqtree    = galaxy_inputs_iqtree(up1$dataset_id),
        kraken2   = galaxy_inputs_kraken2(up1$dataset_id),
        list()
      )
      gid <- switch(tool_id,
        spades    = GALAXY_TOOL_SPADES,
        unicycler = GALAXY_TOOL_UNICYCLER,
        mlst      = GALAXY_TOOL_MLST,
        iqtree    = GALAXY_TOOL_IQTREE,
        kraken2   = GALAXY_TOOL_KRAKEN2,
        tool_id)
      alog("API", paste0("Galaxy run_tool: ", gid),
           paste0("tool_id interno: ", tool_id, "\nhistory_id: ", hist,
                  "\ndataset R1: ", up1$dataset_id, "\ndataset R2: ", up2$dataset_id))
      res <- galaxy_run_tool(key, hist, gid, inputs)
      if (res$success && length(res$job_ids) > 0) {
        add_job(res$job_ids[1], "galaxy", tool_id, step_id,
                paste0(toupper(tool_id), " (Galaxy) — ", muestra),
                output_ids = res$output_ids)
        alog("OK", paste0("Galaxy job enviado: ", toupper(tool_id)),
             paste0("job_id: ", res$job_ids[1], "\noutput_ids: ",
                    paste(res$output_ids, collapse=", ")))
        showNotification(paste("Galaxy job enviado:", res$job_ids[1]), type="message", duration=5)
      } else {
        alog("ERROR", paste0("Galaxy run_tool fallido: ", tool_id),
             res$error %||% "Sin respuesta del servidor")
        showNotification(paste("Error Galaxy:", res$error %||% "Sin respuesta"), type="error", duration=8)
      }
    }
  }

  # ─── Botones de ejecución por paso ───────────────────────
  for (.step in PASOS_IDS) {
    local({
      step <- .step
      observeEvent(input[[paste0("exec_btn_",step)]], {
        run_exec(step, input[[paste0("exec_tool_",step)]], input[[paste0("exec_muestra_",step)]])
      }, ignoreInit=TRUE)
    })
  }

  # ─── "Todas las muestras" por paso ────────────────────────
  for (.step in PASOS_IDS) {
    local({
      step <- .step
      observeEvent(input[[paste0("exec_btn_all_",step)]], {
        tool_spec <- input[[paste0("exec_tool_",step)]]
        muestras  <- isolate(rv$muestras)
        if (!length(muestras)) {
          showNotification("No hay muestras cargadas.", type="warning", duration=4); return()
        }
        run_exec_all(step, tool_spec, muestras)
      }, ignoreInit=TRUE)
    })
  }

  # ─── Helper: ejecutar un tool para TODAS las muestras ─────
  run_exec_all <- function(step_id, tool_spec, muestras_list) {
    if (is.null(tool_spec) || !length(muestras_list)) return()
    parts    <- strsplit(tool_spec, ":")[[1]]
    platform <- parts[1]; tool_id <- parts[2]
    n_total  <- length(muestras_list)
    n_ok     <- 0L; n_err <- 0L

    alog("INFO", paste0("[TODAS] Iniciando: ", tool_id,
                        " para ", n_total, " muestra(s)"),
         paste(muestras_list, collapse = ", "))

    # ── BV-BRC ───────────────────────────────────────────────
    if (platform == "bvbrc") {
      tok  <- rv$bvbrc_token; user <- rv$bvbrc_user
      if (is.null(tok)) { showModal(creds_modal()); return() }
      organism <- input$proy_organismo %||% "Unknown"

      # Sincronizar registro con el workspace (evita re-analizar lo ya hecho)
      tryCatch(bvbrc_sync_registry(tok, user), error = function(e) NULL)

      # Crear carpeta de salida una sola vez antes del loop
      tryCatch(bvbrc_mkdir(tok, ws_out(user)), error = function(e) NULL)

      withProgress(
        message = paste0(tool_id, " — todas las muestras"),
        value   = 0, {

        for (i in seq_along(muestras_list)) {
          muestra <- muestras_list[i]
          incProgress(0, detail = paste0(i,"/",n_total," — subiendo: ", muestra))

          # ── Skip si ya existe un job completado ──────────────
          if (!is.null(bvbrc_get_done_job(rv$jobs, tool_id, muestra))) {
            alog("INFO", paste0("[TODAS] Skip (ya completado): ", tool_id, " — ", muestra))
            n_ok <- n_ok + 1L
            incProgress(1/n_total, detail = paste0("Ya completado: ", muestra))
            next
          }

          # Upload R1 + R2 una vez por muestra
          sfiles <- get_sample_files(muestra)
          r1 <- sfiles$r1; r2 <- sfiles$r2
          ws_r1 <- paste0(ws_home(user), "/.fastq_uploads/", muestra, "_R1.fastq.gz")
          ws_r2 <- NULL  # NULL = single-end por defecto

          if (!is.null(r1) && file.exists(r1)) {
            up1 <- tryCatch(bvbrc_upload_file(tok, user, r1),
                            error = function(e) list(success=FALSE, error=e$message))
            if (isTRUE(up1$success)) {
              ws_r1 <- up1$ws_path
              alog("OK", paste0("[TODAS] R1 subido: ", muestra), up1$log_det)
            } else {
              alog("WARN", paste0("[TODAS] R1 upload fallido: ", muestra),
                   up1$error %||% "sin detalle")
            }
          }
          if (!is.null(r2) && file.exists(r2) && !identical(r2, r1)) {
            up2 <- tryCatch(bvbrc_upload_file(tok, user, r2),
                            error = function(e) list(success=FALSE, error=e$message))
            if (isTRUE(up2$success)) {
              ws_r2 <- up2$ws_path
              alog("OK", paste0("[TODAS] R2 subido: ", muestra), up2$log_det)
            } else {
              alog("WARN", paste0("[TODAS] R2 upload fallido — modo single-end: ", muestra),
                   up2$error %||% "sin detalle")
            }
          } else if (!is.null(r2)) {
            alog("DEBUG", paste0("[TODAS] R2 no disponible — modo single-end: ", muestra))
          }

          incProgress(0.5/n_total,
            detail = paste0(i,"/",n_total," — enviando job: ", muestra))

          params <- tryCatch(switch(tool_id,
            FastqUtils                  = bvbrc_params_fastqutil(user,ws_r1,ws_r2,muestra),
            TaxonomicClassification     = bvbrc_params_taxonomy(user,ws_r1,ws_r2,muestra),
            GenomeDistance              = bvbrc_params_taxonomy(user,ws_r1,ws_r2,muestra),
            Assembly2                   = bvbrc_params_assembly(user,ws_r1,ws_r2,muestra),
            MetagenomicReadMapping      = bvbrc_params_metareads(user,ws_r1,ws_r2,muestra),
            ComprehensiveGenomeAnalysis = bvbrc_params_cga(user,ws_r1,ws_r2,organism,muestra),
            Annotation                  = bvbrc_params_annotation(user,ws_r1,organism,muestra),
            list()
          ), error = function(e) list())

          res <- tryCatch(
            bvbrc_submit_job(tok, tool_id, params, ws_home(user)),
            error = function(e) list(success=FALSE, error=e$message))

          if (isTRUE(res$success)) {
            add_job(res$task_id, "bvbrc", tool_id, step_id,
                    paste0(tool_id, " — ", muestra),
                    output_path = params$output_path,
                    output_file = params$output_file,
                    muestra     = muestra)
            n_ok <- n_ok + 1L
            alog("OK", paste0("[TODAS] Job enviado: ", tool_id, " — ", muestra),
                 paste0("task_id: ", res$task_id))
          } else {
            n_err <- n_err + 1L
            alog("ERROR", paste0("[TODAS] Submit fallido: ", tool_id, " — ", muestra),
                 res$error %||% "sin detalle")
          }
          incProgress(0.5/n_total)
        }
      })   # fin withProgress

    # ── Galaxy — requiere datasets en estado "ok" por muestra ─
    } else if (platform == "galaxy") {
      key  <- rv$galaxy_key; hist <- rv$galaxy_hist
      if (is.null(key)) { showModal(creds_modal()); return() }

      withProgress(
        message = paste0("Galaxy: ", toupper(tool_id), " — todas las muestras"),
        value   = 0, {

        for (i in seq_along(muestras_list)) {
          muestra <- muestras_list[i]
          incProgress(0, detail = paste0(i,"/",n_total," — subiendo: ", muestra))

          sfiles <- get_sample_files(muestra)
          r1 <- sfiles$r1; r2 <- sfiles$r2
          if (is.null(r1) || !file.exists(r1)) {
            alog("WARN", paste0("[TODAS] Galaxy: sin archivo para ", muestra))
            n_err <- n_err + 1L; incProgress(1/n_total); next
          }

          up1 <- tryCatch(galaxy_upload_file(key, hist, r1, "fastqsanger.gz"),
                          error = function(e) list(success=FALSE, error=e$message))
          if (!isTRUE(up1$success)) {
            alog("ERROR", paste0("[TODAS] Galaxy upload R1 fallido: ", muestra), up1$error)
            n_err <- n_err + 1L; incProgress(1/n_total); next
          }

          up2 <- up1
          if (!is.null(r2) && file.exists(r2) && !identical(r2, r1)) {
            up2 <- tryCatch(galaxy_upload_file(key, hist, r2, "fastqsanger.gz"),
                            error = function(e) list(success=FALSE, error=e$message))
            if (!isTRUE(up2$success)) up2 <- up1
          }

          incProgress(0.4/n_total,
            detail = paste0(i,"/",n_total," — esperando dataset: ", muestra))

          w1 <- galaxy_wait_dataset(key, up1$dataset_id, hist, max_wait=300)
          if (!w1$success) {
            alog("ERROR", paste0("[TODAS] Galaxy dataset no listo: ", muestra), w1$error)
            n_err <- n_err + 1L; incProgress(0.6/n_total); next
          }
          if (up2$dataset_id != up1$dataset_id) {
            w2 <- galaxy_wait_dataset(key, up2$dataset_id, hist, max_wait=300)
            if (!w2$success) up2 <- up1
          }

          incProgress(0.3/n_total,
            detail = paste0(i,"/",n_total," — lanzando: ", muestra))

          inputs <- switch(tool_id,
            spades    = galaxy_inputs_spades(up1$dataset_id, up2$dataset_id),
            unicycler = galaxy_inputs_unicycler(up1$dataset_id, up2$dataset_id),
            mlst      = galaxy_inputs_mlst(up1$dataset_id),
            iqtree    = galaxy_inputs_iqtree(up1$dataset_id),
            kraken2   = galaxy_inputs_kraken2(up1$dataset_id),
            list())
          gid <- switch(tool_id,
            spades    = GALAXY_TOOL_SPADES,
            unicycler = GALAXY_TOOL_UNICYCLER,
            mlst      = GALAXY_TOOL_MLST,
            iqtree    = GALAXY_TOOL_IQTREE,
            kraken2   = GALAXY_TOOL_KRAKEN2,
            tool_id)

          res <- tryCatch(galaxy_run_tool(key, hist, gid, inputs),
                          error = function(e) list(success=FALSE, error=e$message))

          if (isTRUE(res$success) && length(res$job_ids) > 0) {
            add_job(res$job_ids[1], "galaxy", tool_id, step_id,
                    paste0(toupper(tool_id), " (Galaxy) — ", muestra),
                    output_ids = res$output_ids, muestra = muestra)
            n_ok <- n_ok + 1L
            alog("OK", paste0("[TODAS] Galaxy job enviado: ", tool_id, " — ", muestra),
                 paste0("job_id: ", res$job_ids[1]))
          } else {
            n_err <- n_err + 1L
            alog("ERROR", paste0("[TODAS] Galaxy run fallido: ", tool_id, " — ", muestra),
                 res$error %||% "sin detalle")
          }
          incProgress(0.3/n_total)
        }
      })   # fin withProgress

    # ── EBI JDispatcher ──────────────────────────────────────
    } else if (platform == "ebi") {
      withProgress(
        message = paste0("EBI: ", toupper(tool_id), " — todas las muestras"),
        value   = 0, {

        for (i in seq_along(muestras_list)) {
          muestra <- muestras_list[i]
          incProgress(1/n_total, detail = paste0(i,"/",n_total,": ", muestra))
          sfiles <- get_sample_files(muestra)
          ffile  <- sfiles$r1
          if (is.null(ffile) || !file.exists(ffile)) { n_err <- n_err+1L; next }
          raw   <- readLines(ffile, n=100, warn=FALSE)
          fasta <- if (any(grepl("^>",raw)))
            paste(head(raw,20), collapse="\n")
          else
            paste0(">", muestra, "\n", paste(head(raw,5), collapse=""))
          res <- switch(tool_id,
            muscle  = ebi_muscle(fasta),
            clustalo= ebi_clustalo(fasta),
            mafft   = ebi_mafft(fasta),
            list(success=FALSE, error="Herramienta EBI desconocida"))
          if (isTRUE(res$success)) {
            add_job(res$job_id, "ebi", tool_id, step_id,
                    paste0(toupper(tool_id), " — ", muestra))
            n_ok <- n_ok + 1L
            alog("OK", paste0("[TODAS] EBI job enviado: ", toupper(tool_id), " — ", muestra),
                 paste0("job_id: ", res$job_id))
          } else {
            n_err <- n_err + 1L
            alog("ERROR", paste0("[TODAS] EBI fallido: ", tool_id, " — ", muestra),
                 res$error %||% "sin detalle")
          }
        }
      })

    # ── NCBI BLAST ───────────────────────────────────────────
    } else if (platform == "ncbi") {
      withProgress(
        message = paste0("NCBI BLAST — todas las muestras"),
        value   = 0, {

        for (i in seq_along(muestras_list)) {
          muestra <- muestras_list[i]
          incProgress(1/n_total, detail = paste0(i,"/",n_total,": ", muestra))
          sfiles <- get_sample_files(muestra)
          ffile  <- sfiles$r1
          if (is.null(ffile) || !file.exists(ffile)) { n_err <- n_err+1L; next }
          raw   <- readLines(ffile, n=50, warn=FALSE)
          fasta <- paste(head(raw,10), collapse="\n")
          res   <- ncbi_blast_submit(fasta, program=tool_id)
          if (isTRUE(res$success)) {
            add_job(res$rid, "ncbi", tool_id, step_id,
                    paste0("BLAST — ", muestra))
            n_ok <- n_ok + 1L
            alog("OK", paste0("[TODAS] NCBI BLAST enviado: ", muestra),
                 paste0("RID: ", res$rid))
          } else {
            n_err <- n_err + 1L
            alog("ERROR", paste0("[TODAS] NCBI BLAST fallido: ", muestra),
                 res$error %||% "sin detalle")
          }
        }
      })

    # ── CGE / ResFinder ──────────────────────────────────────
    } else if (platform == "cge") {
      withProgress(
        message = "ResFinder (CGE) — todas las muestras",
        value   = 0, {

        for (i in seq_along(muestras_list)) {
          muestra <- muestras_list[i]
          incProgress(1/n_total, detail = paste0(i,"/",n_total,": ", muestra))
          sfiles <- get_sample_files(muestra)
          ffile  <- sfiles$r1
          if (is.null(ffile) || !file.exists(ffile)) { n_err <- n_err+1L; next }
          res <- resfinder_submit(fasta_file=ffile, species="other")
          if (isTRUE(res$success) && !is.null(res$job_id)) {
            add_job(res$job_id, "cge", "resfinder", step_id,
                    paste0("ResFinder — ", muestra))
            n_ok <- n_ok + 1L
            alog("OK", paste0("[TODAS] ResFinder enviado: ", muestra),
                 paste0("job_id: ", res$job_id))
          } else if (isTRUE(res$success) && !is.null(res$result_direct)) {
            parsed <- tryCatch(
              jsonlite::fromJSON(res$result_direct, simplifyVector=FALSE),
              error=function(e) NULL)
            if (!is.null(parsed)) rv$resfinder_df <- resfinder_to_table(parsed)
            n_ok <- n_ok + 1L
          } else {
            n_err <- n_err + 1L
            alog("ERROR", paste0("[TODAS] ResFinder fallido: ", muestra),
                 res$error %||% "sin detalle")
          }
        }
      })
    }

    # ── Notificación final ────────────────────────────────────
    msg_all <- paste0(n_ok, "/", n_total, " jobs enviados")
    if (n_err > 0) msg_all <- paste0(msg_all, " (", n_err, " error(es) — ver Log)")
    alog(if(n_err==0)"OK" else "WARN",
         paste0("[TODAS] ", tool_id, " completado"),
         paste0("OK: ", n_ok, " | ERR: ", n_err, " | Total: ", n_total))
    showNotification(
      tagList(icon("layer-group"), " ", msg_all),
      type     = if (n_err == 0) "message" else "warning",
      duration = 10)
  }

  # ─── Renderizar panel de ejecución por paso ───────────────
  for (.step in PASOS_IDS) {
    local({
      step <- .step
      output[[paste0("exec_panel_",step)]] <- renderUI({
        muestras_disp <- rv$muestras
        tools_list    <- EXEC_TOOLS[[step]]
        jobs_step     <- Filter(function(j) j$step==step, rv$jobs)

        tagList(
          hr(),
          div(style="background:#fffbea;border:1px solid #e8c84a;border-radius:6px;padding:12px 14px;",
            div(style="font-weight:bold;font-size:13px;color:#856404;margin-bottom:8px;",
                icon("bolt"), " Ejecución directa via API"),
            if (length(muestras_disp)==0)
              div(class="alert alert-warning",style="padding:6px 10px;font-size:12px;margin:0;",
                  icon("exclamation-triangle"), " Sube archivos FASTQ en Inicio primero.")
            else
              tagList(
                fluidRow(
                  column(4, selectInput(paste0("exec_muestra_",step), "Muestra:", choices=muestras_disp, width="100%")),
                  column(5, selectInput(paste0("exec_tool_",step),    "Herramienta:", choices=tools_list, width="100%")),
                  column(3, br(), actionButton(paste0("exec_btn_",step), tagList(icon("play")," Enviar"),
                                               class="btn-warning btn-sm", style="width:100%;margin-top:5px;"))
                ),
                actionButton(paste0("exec_btn_all_",step),
                  tagList(icon("layer-group"),
                          " Enviar para todas las muestras (",
                          length(muestras_disp), ")"),
                  class = "btn-info btn-sm",
                  style = "width:100%;margin-top:5px;")
              ),
            if (length(jobs_step)>0)
              div(style="margin-top:8px;",
                lapply(rev(jobs_step), function(j) {
                  cls <- if (j$status %in% c("completed","FINISHED","ok")) "job-card job-done"
                         else if (j$status %in% c("error","ERROR","failed")) "job-card job-error"
                         else "job-card"
                  ic  <- if (j$status %in% c("completed","FINISHED","ok"))    icon("check-circle", style="color:#27ae60;")
                         else if (j$status %in% c("error","ERROR","failed"))  icon("times-circle", style="color:#e74c3c;")
                         else if (j$status == "queued")                        icon("clock", style="color:#888;")
                         else                                                   icon("spinner", class="fa-spin", style="color:#f39c12;")
                  div(class=cls,
                    ic, " ", tags$b(j$label), " — ", em(j$status),
                    tags$small(style="float:right;color:#666;", format(j$time,"%H:%M:%S")),
                    if (!is.null(j$result_file))
                      div(style="margin-top:4px;",
                        downloadLink(paste0("dl_res_",gsub("[^a-zA-Z0-9]","",j$job_id)),
                                     tagList(icon("download")," Descargar resultado"))),
                    if (j$status %in% c("completed","FINISHED","ok","done") &&
                        j$platform == "bvbrc") {
                      has_res <- !is.null(rv$job_results[[j$job_id]])
                      div(style="margin-top:5px;",
                        tags$button(
                          class = paste("btn btn-xs", if(has_res) "btn-success" else "btn-info"),
                          onclick = paste0("Shiny.setInputValue('view_result_job','",
                                          j$job_id, "',{priority:'event'})"),
                          HTML(paste0(
                            '<i class="fa fa-', if(has_res) 'chart-bar' else 'sync-alt', '"></i> ',
                            if(has_res) 'Ver resultados' else 'Buscar resultados'
                          ))
                        )
                      )
                    }
                  )
                })
              )
          )
        )
      })
    })
  }

  # ─── Resultados EBI en alineamiento ───────────────────────
  output$ebi_result_display <- renderUI({
    res <- rv$ebi_result; if (is.null(res)) return(NULL)
    div(class="alert alert-success",style="margin-top:10px;",
      icon("check-circle"), strong(" Alineamiento completado (EBI)"), br(),
      tags$pre(style="font-size:11px;max-height:200px;overflow-y:auto;background:#f8f9fa;border-radius:4px;padding:8px;", substr(res,1,1000)),
      downloadButton("dl_ebi_aln", "Descargar alineamiento (.fasta)", class="btn-success btn-sm")
    )
  })
  output$dl_ebi_aln <- downloadHandler(
    filename = function() paste0("alineamiento_",format(Sys.Date(),"%Y%m%d"),".fasta"),
    content  = function(file) writeLines(rv$ebi_result%||%"", file)
  )

  # ─── Resultados ResFinder en resistoma ────────────────────
  output$resfinder_result_display <- renderUI({
    df <- rv$resfinder_df; if (is.null(df)||nrow(df)==0) return(NULL)
    tagList(
      div(class="alert alert-danger",style="padding:6px 10px;", icon("shield-alt"),
          strong(paste0(" ResFinder: ",nrow(df)," genes AMR detectados"))),
      renderDT(datatable(df, rownames=FALSE, options=list(dom='t',pageLength=10,scrollX=TRUE),
                         colnames=c("Gen AMR","Identidad (%)","Cobertura (%)","Antibiótico")))
    )
  })

  # ─── Barra de progreso ────────────────────────────────────
  output$barra_progreso_global <- renderUI({
    tab <- req(input$sidebar); if (tab %in% c("inicio","reporte")) return(NULL)
    compl <- rv$completados; n_done <- sum(compl,na.rm=TRUE); pct <- round(n_done/length(PASOS_IDS)*100)
    step_items <- lapply(seq_along(PASOS_IDS), function(i) {
      pid <- PASOS_IDS[i]
      estado <- if(isTRUE(compl[pid]))"done" else if(tab==pid)"active" else "pending"
      bg  <- switch(estado, done="background:#27ae60;color:white;", active="background:#2e86c1;color:white;box-shadow:0 0 0 4px #aed6f1;", "background:#e9ecef;color:#868e96;")
      ico <- if(estado=="done")"✓" else as.character(i)
      con <- if(i<length(PASOS_IDS)) div(class="step-connector",style=paste0("background:",if(isTRUE(compl[pid]))"#27ae60" else "#e9ecef",";")) else NULL
      tagList(div(class="step-dot", onclick=paste0("Shiny.setInputValue('nav_click','",pid,"',{priority:'event'})"),
        div(class=paste("step-circle",estado),style=bg,ico),
        div(class=paste("step-label",estado),PASOS_LABELS[i])), con)
    })
    div(id="barra_progreso_global",
      div(style="display:flex;align-items:center;overflow-x:auto;padding-bottom:4px;", step_items),
      div(style="margin-top:6px;",
        div(style="background:#e9ecef;border-radius:4px;height:6px;",
          div(style=paste0("background:#27ae60;height:6px;border-radius:4px;width:",pct,"%;transition:width 0.4s;"))),
        div(style="font-size:11px;color:#666;margin-top:3px;text-align:right;", paste0(n_done," de ",length(PASOS_IDS)," — ",pct,"%"))
      ))
  })
  output$sidebar_progreso <- renderUI({
    compl <- rv$completados; n <- sum(compl,na.rm=TRUE); tot <- length(compl); pct <- round(n/tot*100)
    div(div(style="font-size:11px;color:#aaa;margin-bottom:3px;", paste0("Progreso: ",n,"/",tot," (",pct,"%)")),
        div(class="sidebar-prog-bar", div(class="sidebar-prog-fill",style=paste0("width:",pct,"%;"))))
  })
  observeEvent(input$nav_click, { req(input$nav_click); updateTabItems(session,"sidebar",selected=input$nav_click) })

  for (i in seq_along(PASOS_IDS)) {
    local({
      pid     <- PASOS_IDS[i]
      prev_id <- if(i==1)"inicio" else PASOS_IDS[i-1]
      next_id <- if(i==length(PASOS_IDS))"reporte" else PASOS_IDS[i+1]
      observeEvent(input[[paste0("btn_sig_",pid)]], { rv$completados[pid]<-TRUE; updateCheckboxInput(session,paste0("chk_",pid),value=TRUE); updateTabItems(session,"sidebar",selected=next_id) }, ignoreInit=TRUE)
      observeEvent(input[[paste0("btn_ant_",pid)]], { updateTabItems(session,"sidebar",selected=prev_id) }, ignoreInit=TRUE)
      observeEvent(input[[paste0("chk_",pid)]],     { rv$completados[pid]<-isTRUE(input[[paste0("chk_",pid)]]) }, ignoreInit=TRUE)
    })
  }

  # ─── Tablas y gráficas ────────────────────────────────────
  # ─── Estado vacío: sin datos dummy ────────────────────────
  # Placeholder para gráficas/tablas hasta que existan resultados REALES.
  empty_plot <- function(msg = "Sin datos aún — ejecuta esta etapa para ver resultados reales.") {
    plotly::plot_ly(type = "scatter", mode = "markers", x = numeric(0), y = numeric(0)) |>
      plotly::add_annotations(text = msg, x = 0.5, y = 0.5,
                              xref = "paper", yref = "paper", showarrow = FALSE,
                              font = list(size = 14, color = "#7f8c8d")) |>
      plotly::layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
                     plot_bgcolor = "#fafafa", paper_bgcolor = "#fafafa")
  }
  empty_dt <- function(msg = "Sin datos aún — los resultados aparecerán aquí al completarse la etapa.") {
    DT::datatable(data.frame(` ` = msg, check.names = FALSE),
                  rownames = FALSE, options = list(dom = 't', ordering = FALSE))
  }

  observeEvent(input$add_calidad, {
    req(input$cq_muestra)
    rv$calidad <- rbind(rv$calidad, data.frame(Muestra=input$cq_muestra,Q30=input$cq_q30,GC=input$cq_gc,Reads_M=input$cq_reads,Duplicados=input$cq_dup,Estado=ifelse(!is.na(input$cq_q30)&&input$cq_q30>=75,"APROBADO","REVISAR"),stringsAsFactors=FALSE))
    showNotification(paste("Muestra",input$cq_muestra,"agregada."),type="message",duration=3)
  })
  output$tabla_calidad <- renderDT({
    if(nrow(rv$calidad)==0) return(empty_dt())
    datatable(rv$calidad,editable=TRUE,rownames=FALSE,options=list(dom='t',pageLength=10,scrollX=TRUE),colnames=c("Muestra","% Q30","% GC","Reads (M)","% Dup.","Estado"))%>%formatStyle("Estado",backgroundColor=styleEqual(c("APROBADO","REVISAR"),c("#d5f5e3","#fadbd8")),fontWeight="bold") })
  output$plot_calidad <- renderPlotly({
    if(nrow(rv$calidad)==0) return(empty_plot())
    df <- rv$calidad
    p <- ggplot(df,aes(x=Muestra,y=Q30,fill=Estado))+geom_col(width=0.6)+geom_hline(yintercept=75,linetype="dashed",color="red",linewidth=1)+annotate("text",x=0.7,y=76.5,label="Min 75%",color="red",size=3.5)+scale_fill_manual(values=c("APROBADO"="#27AE60","REVISAR"="#E74C3C"))+scale_y_continuous(limits=c(0,100))+labs(title="% Q30 por muestra",y="% Q30",x="")+theme_minimal(base_size=13)
    ggplotly(p)%>%layout(legend=list(x=0.01,y=0.99))
  })
  output$plot_gc <- renderPlotly({
    if(nrow(rv$calidad)==0||all(is.na(rv$calidad$GC))) return(empty_plot())
    df <- rv$calidad
    p <- ggplot(df,aes(x=Muestra,y=GC,fill=Muestra))+geom_col(width=0.6,show.legend=FALSE)+geom_hline(yintercept=c(38,40),linetype="dotted",color="#888")+scale_y_continuous(limits=c(30,50))+scale_fill_brewer(palette="Set2")+labs(title="GC (%)",y="%GC",x="")+theme_minimal(base_size=13)
    ggplotly(p)
  })

  observeEvent(input$add_taxonomia, {
    req(input$tax_muestra)
    rv$taxonomia <- rbind(rv$taxonomia, data.frame(Muestra=input$tax_muestra,Especie=input$tax_especie,Confianza=input$tax_conf,Contaminacion=input$tax_contam,stringsAsFactors=FALSE))
    showNotification(paste("Muestra",input$tax_muestra,"agregada."),type="message",duration=3)
  })
  output$tabla_taxonomia <- renderDT({
    if(nrow(rv$taxonomia)==0) return(empty_dt())
    datatable(rv$taxonomia,editable=TRUE,rownames=FALSE,options=list(dom='t',pageLength=10)) })
  output$plot_taxonomia_pie <- renderPlotly({
    if(nrow(rv$taxonomia)==0||all(is.na(rv$taxonomia$Especie))) return(empty_plot())
    df <- rv$taxonomia%>%count(Especie)
    plot_ly(df,labels=~Especie,values=~n,type="pie",marker=list(colors=c("#2ECC71","#3498DB","#E74C3C","#F39C12","#9B59B6")))%>%layout(title="Distribución taxonómica")
  })
  output$plot_contaminacion <- renderPlotly({
    if(nrow(rv$taxonomia)==0||all(is.na(rv$taxonomia$Contaminacion))) return(empty_plot())
    df <- rv$taxonomia
    p <- ggplot(df,aes(x=Muestra,y=Contaminacion,fill=Contaminacion>5))+geom_col(width=0.6)+geom_hline(yintercept=5,linetype="dashed",color="red",linewidth=1)+scale_fill_manual(values=c("FALSE"="#3498DB","TRUE"="#E74C3C"),name=">5%")+labs(title="% Contaminación",y="%",x="")+theme_minimal(base_size=13)
    ggplotly(p)
  })

  observeEvent(input$add_ensamblado, {
    req(input$ens_muestra)
    rv$ensamblado <- rbind(rv$ensamblado, data.frame(Muestra=input$ens_muestra,Cobertura=input$ens_cobert,Profundidad=input$ens_prof,Contigs=input$ens_contigs,N50_kb=input$ens_n50,Tamano_Mb=input$ens_tam,GC=input$ens_gc,stringsAsFactors=FALSE))
    showNotification(paste("Muestra",input$ens_muestra,"agregada."),type="message",duration=3)
  })
  output$tabla_ensamblado <- renderDT({
    if(nrow(rv$ensamblado)==0) return(empty_dt())
    datatable(rv$ensamblado,editable=TRUE,rownames=FALSE,options=list(dom='t',pageLength=10,scrollX=TRUE)) })
  output$plot_ensamblado <- renderPlotly({
    if(nrow(rv$ensamblado)==0) return(empty_plot())
    df <- rv$ensamblado
    p1 <- plot_ly(df,x=~Muestra,y=~Cobertura,type="bar",name="Cobertura (%)",marker=list(color="#2ECC71"))
    p1 <- add_trace(p1,x=df$Muestra,y=rep(95,nrow(df)),type="scatter",mode="lines",inherit=FALSE,name="Min 95%",line=list(color="red",dash="dash"))
    p2 <- plot_ly(df,x=~Muestra,y=~N50_kb,type="bar",name="N50 (kb)",marker=list(color="#3498DB"))
    subplot(p1,p2,nrows=1,shareX=TRUE,titleY=TRUE)%>%layout(title="Ensamblado",showlegend=TRUE)
  })
  output$eval_ensamblado <- renderUI({
    if(nrow(rv$ensamblado)==0)return(div(class="alert alert-info","Ingrese resultados."))
    df <- rv$ensamblado; ok <- sum(!is.na(df$Cobertura)&df$Cobertura>=95&!is.na(df$N50_kb)&df$N50_kb>=50&!is.na(df$Contigs)&df$Contigs<200,na.rm=TRUE)
    div(div(class="metric-value",paste0(ok,"/",nrow(df))),div(class="metric-label","muestras aprobadas"),br(),
      if(ok==nrow(df))div(class="alert alert-success",icon("check-circle")," Todos cumplen criterios.") else div(class="alert alert-warning",icon("exclamation-triangle"),paste(nrow(df)-ok,"requieren revisión.")))
  })

  output$plot_anotacion <- renderPlotly({
    txt <- input$anot_cog
    df  <- tryCatch({
      if(!is.null(txt)&&nchar(trimws(txt))>10){
        ln <- strsplit(trimws(txt),"\n")[[1]]; ln <- ln[nchar(trimws(ln))>0]
        pt <- strsplit(ln,":")
        tmp <- data.frame(Cat=trimws(sapply(pt,`[[`,1)),Genes=suppressWarnings(as.numeric(trimws(gsub("[^0-9]","",sapply(pt,function(x)ifelse(length(x)>1,x[2],"0")))))),stringsAsFactors=FALSE)
        tmp[!is.na(tmp$Genes)&tmp$Genes>0,]
      } else NULL
    },error=function(e)NULL)
    if(is.null(df)||nrow(df)==0) return(empty_plot("Sin datos de anotación aún — completa la etapa 4."))
    p <- ggplot(df,aes(x=reorder(Cat,Genes),y=Genes,fill=Genes))+geom_col()+coord_flip()+scale_fill_gradient(low="#AED6F1",high="#1A5276")+labs(title="COG",x="",y="Genes")+theme_minimal(base_size=12)+theme(legend.position="none")
    ggplotly(p)
  })
  output$resumen_anotacion <- renderUI({
    div(fluidRow(column(6,div(class="metric-value",ifelse(is.na(input$anot_cds),"—",format(input$anot_cds,big.mark=","))),div(class="metric-label","CDS")),column(6,div(class="metric-value",ifelse(is.na(input$anot_funcion),"—",paste0(input$anot_funcion,"%"))),div(class="metric-label","Con función"))),fluidRow(column(6,div(class="metric-value",ifelse(is.na(input$anot_rrna),"—",input$anot_rrna)),div(class="metric-label","ARNr")),column(6,div(class="metric-value",ifelse(is.na(input$anot_trna),"—",input$anot_trna)),div(class="metric-label","ARNt"))))
  })

  output$tabla_resistencia <- renderDT({
    df <- rv$resfinder_df
    if(is.null(df)||!is.data.frame(df)||nrow(df)==0)
      return(empty_dt("Sin genes de resistencia aún — ejecuta el resistoma (MetagenomicReadMapping / ResFinder)."))
    datatable(df,editable=TRUE,rownames=FALSE,options=list(dom='tp',pageLength=15,scrollX=TRUE))
  })
  output$plot_resistoma <- renderPlotly({
    df <- rv$resfinder_df
    if(is.null(df)||!is.data.frame(df)||nrow(df)==0||!("Antibiotico"%in%names(df)))
      return(empty_plot("Sin datos de resistoma aún — completa la etapa 5."))
    ab <- df$Antibiotico; ab[is.na(ab)|nchar(trimws(ab))==0] <- "Sin clasificar"
    tab <- as.data.frame(table(Clase=ab), stringsAsFactors=FALSE)
    if(nrow(tab)==0) return(empty_plot("Sin datos de resistoma aún — completa la etapa 5."))
    p <- ggplot(tab,aes(x=reorder(Clase,Freq),y=Freq,fill=Freq))+geom_col()+coord_flip()+
      scale_fill_gradient(low="#F5B7B1",high="#C0392B")+
      labs(title="Genes de resistencia por clase de antibiótico",x="",y="N° de genes")+
      theme_minimal(base_size=12)+theme(legend.position="none")
    ggplotly(p)
  })
  output$clasif_badge <- renderUI({
    cl <- input$res_clasif; if(is.null(cl)||length(cl)==0)return(div(class="alert alert-light",style="font-size:12px;","Seleccione clasificación"))
    col <- if("PDR"%in%cl)"danger" else if("XDR"%in%cl)"danger" else "warning"
    div(class=paste0("alert alert-",col),style="margin:0;",icon("exclamation-triangle"),strong(paste(cl,collapse=" / ")))
  })

  observeEvent(input$add_ali, {
    req(input$ali_par)
    rv$alineamiento <- rbind(rv$alineamiento,data.frame(Par=input$ali_par,Herramienta=input$ali_tool,Identidad=input$ali_ident,Cobertura=input$ali_cobert,Gaps=input$ali_gaps,SNPs=input$ali_snps,stringsAsFactors=FALSE))
    showNotification("Comparación agregada.",type="message",duration=3)
  })
  output$tabla_alineamiento <- renderDT({
    if(nrow(rv$alineamiento)==0) return(empty_dt())
    datatable(rv$alineamiento,editable=TRUE,rownames=FALSE,options=list(dom='t',pageLength=10,scrollX=TRUE)) })
  output$plot_alineamiento  <- renderPlotly({
    df <- rv$alineamiento
    if(nrow(df)==0||all(is.na(df$Identidad))) return(empty_plot("Sin alineamientos aún — completa la etapa 6."))
    p <- ggplot(df,aes(x=reorder(Par,Identidad),y=Identidad,fill=Identidad))+geom_col()+coord_flip()+
      scale_fill_gradient(low="#AED6F1",high="#1A5276")+
      labs(title="% Identidad por comparación",x="",y="% Identidad")+
      theme_minimal(base_size=12)+theme(legend.position="none")
    ggplotly(p)
  })

  observeEvent(input$add_mlst, {
    req(input$mlst_muestra)
    rv$mlst <- rbind(rv$mlst,data.frame(Muestra=input$mlst_muestra,ST=input$mlst_st,Alelos=input$mlst_alelos,Esquema=input$mlst_esquema,stringsAsFactors=FALSE))
    showNotification(paste("MLST",input$mlst_muestra,"agregado."),type="message",duration=3)
  })
  output$tabla_mlst <- renderDT({
    if(nrow(rv$mlst)==0) return(empty_dt())
    datatable(rv$mlst,editable=TRUE,rownames=FALSE,options=list(dom='t',pageLength=10))
  })
  output$plot_mlst <- renderPlotly({
    if(nrow(rv$mlst)==0) return(empty_plot("Sin datos de MLST aún — completa la etapa 7."))
    df <- as.data.frame(table(ST=rv$mlst$ST),stringsAsFactors=FALSE)
    p <- ggplot(df,aes(x=reorder(ST,-Freq),y=Freq,fill=ST))+geom_col(width=0.5,show.legend=FALSE)+geom_text(aes(label=Freq),vjust=-0.4,size=5,fontface="bold")+scale_fill_brewer(palette="Set1")+labs(title="Sequence Types",x="ST",y="n")+theme_minimal(base_size=14)
    ggplotly(p)
  })

  # Mapa genome_id → nombre de muestra (BV-BRC etiqueta todos los
  # genomas propios con el mismo nombre científico)
  mapa_genomas <- reactive({
    gids <- isolate(rv$auto_pipe)$genome_ids %||% list()
    if (!length(gids)) return(NULL)
    setNames(names(gids), unlist(gids, use.names = FALSE))
  })

  output$plot_filogenia <- renderPlotly({
    nwk <- rv$newick %||% input$fil_newick
    if (is.null(nwk) || nchar(trimws(nwk)) == 0)
      return(empty_plot("Sin árbol filogenético aún — completa la etapa 8."))
    p <- tryCatch(
      plot_newick(nwk, titulo = NULL,
                  resaltar  = c(rv$muestras, unlist(mapa_genomas(), use.names = FALSE)),
                  renombrar = mapa_genomas()),
      error = function(e) NULL)
    if (is.null(p))
      return(empty_plot("No se pudo interpretar el Newick; revísalo en el campo de abajo."))
    ggplotly(p) %>% layout(showlegend = FALSE)
  })

  observeEvent(input$add_ani, {
    req(input$ani_m1,input$ani_m2)
    rv$ani <- rbind(rv$ani,data.frame(Muestra1=input$ani_m1,Muestra2=input$ani_m2,ANI=input$ani_val,stringsAsFactors=FALSE))
    showNotification("ANI agregado.",type="message",duration=3)
  })
  output$tabla_ani <- renderDT({
    if(nrow(rv$ani)==0) return(empty_dt())
    datatable(rv$ani,editable=TRUE,rownames=FALSE,options=list(dom='t',pageLength=10))
  })
  output$plot_ani <- renderPlotly({
    df <- rv$ani
    if(nrow(df)==0||all(is.na(df$ANI))) return(empty_plot("Sin datos de ANI aún — completa la etapa 9."))
    # Matriz simétrica a partir de los pares reales ingresados/importados
    m <- sort(unique(c(df$Muestra1, df$Muestra2)))
    mat <- matrix(NA_real_, length(m), length(m), dimnames=list(m,m))
    diag(mat) <- 100
    for(i in seq_len(nrow(df))){
      a<-df$Muestra1[i]; b<-df$Muestra2[i]; v<-df$ANI[i]
      if(a%in%m && b%in%m){ mat[a,b]<-v; mat[b,a]<-v }
    }
    plot_ly(z=mat,x=m,y=m,type="heatmap",colorscale=list(c(0,"#FEF9C3"),c(0.5,"#F39C12"),c(1,"#C0392B")),zmin=min(96,min(mat,na.rm=TRUE)),zmax=100,text=round(mat,2),texttemplate="%{text}",colorbar=list(title="ANI (%)"))%>%layout(title="ANI")
  })
  output$plot_pangenoma <- renderPlotly({
    if(is.na(input$gc_core%||%NA)&&is.na(input$gc_acc%||%NA)&&is.na(input$gc_uniq%||%NA))
      return(empty_plot("Sin datos de pangenoma aún — ingresa/importa core, accesorios y únicos."))
    core<-ifelse(is.na(input$gc_core),0,input$gc_core); acc<-ifelse(is.na(input$gc_acc),0,input$gc_acc); uniq<-ifelse(is.na(input$gc_uniq),0,input$gc_uniq)
    plot_ly(labels=c("Core","Accesorios","Únicos"),values=c(core,acc,uniq),type="pie",hole=0.42,marker=list(colors=c("#2ECC71","#3498DB","#E74C3C")))%>%layout(title="Pangenoma")
  })

  # ─── Reporte ──────────────────────────────────────────────
  output$preview_reporte <- renderUI({
    div(h4(icon("file-code")," ",input$rep_titulo),hr(),fluidRow(column(6,p(strong("Proyecto:"),input$proy_nombre)),column(6,p(strong("Organismo:"),input$proy_organismo))),fluidRow(column(6,p(strong("Analista:"),input$rep_autor)),column(6,p(strong("Fecha:"),format(input$rep_fecha,"%d/%m/%Y")))),if(nchar(trimws(input$rep_resumen))>0)div(hr(),h5("Resumen:"),p(input$rep_resumen)),hr(),h5("Secciones:"),tags$ol(lapply(input$rep_secciones,function(s)tags$li(paste0(herramientas_info[[s]]$numero,". ",herramientas_info[[s]]$titulo)))))
  })
  observeEvent(input$generar_reporte, {
    showModal(modalDialog(title="Generando reporte...",tags$p(icon("spinner",class="fa-spin")," Espere..."),footer=NULL))
    tryCatch({
      params_rep <- list(titulo=input$rep_titulo,autor=input$rep_autor,fecha=format(input$rep_fecha,"%d de %B de %Y"),proyecto=input$proy_nombre,organismo=input$proy_organismo,institucion=input$proy_institucion,analista=input$proy_analista,resumen=input$rep_resumen,secciones=input$rep_secciones,muestras=if(length(rv$muestras)>0)rv$muestras else "(sin muestras cargadas)",calidad=rv$calidad,calidad_notas=input$calidad_notas,taxonomia=rv$taxonomia,tax_notas=input$tax_notas,ensamblado=rv$ensamblado,ens_notas=input$ens_notas,anot_cds=input$anot_cds,anot_funcion=input$anot_funcion,anot_cog=input$anot_cog,anot_notas=input$anot_notas,res_genes=input$res_genes,res_mutaciones=input$res_mutaciones,res_clasif=input$res_clasif,alineamiento=rv$alineamiento,ali_notas=input$ali_notas,mlst=rv$mlst,mlst_esquema=input$mlst_esquema,mlst_notas=input$mlst_notas,fil_metodo=input$fil_metodo,fil_modelo=input$fil_modelo,fil_bootstrap=input$fil_bootstrap,fil_newick=input$fil_newick,fil_notas=input$fil_notas,ani=rv$ani,gc_core=input$gc_core,gc_acc=input$gc_acc,gc_uniq=input$gc_uniq,gc_notas=input$gc_notas)
      fmt <- input$rep_formato
      ext <- switch(fmt,"html_document"=".html","pdf_document"=".pdf","word_document"=".docx")
      outfile <- file.path(tempdir(),paste0("reporte_bioinformatico_",format(Sys.Date(),"%Y%m%d"),ext))
      rmarkdown::render(input="reporte_template.Rmd",output_format=fmt,output_file=outfile,params=list(datos=params_rep),envir=new.env(parent=globalenv()),quiet=TRUE)
      rv$reporte_path <- outfile; removeModal()
      showNotification("Reporte generado exitosamente.",type="message",duration=5)
    }, error=function(e){ removeModal(); showNotification(paste("Error:",e$message),type="error",duration=10) })
  })
  output$reporte_descarga <- renderUI({ req(rv$reporte_path); if(file.exists(rv$reporte_path)) downloadButton("dl_reporte","Descargar Reporte",class="btn-primary btn-lg",style="width:100%;margin-top:8px;") })
  output$dl_reporte <- downloadHandler(
    filename = function(){ ext<-switch(input$rep_formato,"html_document"=".html","pdf_document"=".pdf","word_document"=".docx"); paste0("reporte_bioinformatico_",format(Sys.Date(),"%Y%m%d"),ext) },
    content  = function(file){ req(rv$reporte_path); file.copy(rv$reporte_path,file) }
  )

  # ============================================================
  # ASISTENTE DE INICIO (carpeta → credenciales → iniciar)
  # ============================================================
  shinyFiles::shinyDirChoose(input, "wiz_carpeta",
    roots = c(Datos = DATOS_DIR, Inicio = path.expand("~"),
              `C:` = "C:/"),
    session = session)

  observeEvent(input$wiz_carpeta, {
    rutas <- c(Datos = DATOS_DIR, Inicio = path.expand("~"), `C:` = "C:/")
    ruta <- tryCatch(shinyFiles::parseDirPath(rutas, input$wiz_carpeta),
                     error = function(e) character(0))
    if (!length(ruta) || !nzchar(ruta)) return()
    cargar_fastq_de_carpeta(as.character(ruta))
  }, ignoreInit = TRUE)

  output$wiz_carpeta_info <- renderUI({
    n <- length(rv$muestras)
    if (n == 0)
      return(div(class = "alert alert-warning", style = "padding:6px 10px;font-size:12px;margin:0;",
                 icon("exclamation-triangle"), " Sin carpeta seleccionada."))
    div(class = "alert alert-success", style = "padding:6px 10px;font-size:12px;margin:0;",
      icon("check-circle"), strong(paste0(" ", n, " muestra(s)")), br(),
      tags$small(paste(head(rv$muestras, 6), collapse = ", "),
                 if (n > 6) paste0(" +", n - 6, " más")),
      if (!is.null(rv$fastq_dir))
        tags$small(style = "display:block;color:#555;margin-top:3px;",
                   icon("folder"), " ", basename(rv$fastq_dir)))
  })

  observeEvent(input$wiz_creds, showModal(creds_modal()), ignoreInit = TRUE)

  output$wiz_creds_info <- renderUI({
    ok_b <- !is.null(rv$bvbrc_token); ok_g <- !is.null(rv$galaxy_key)
    est <- function(ok, txt) div(style = "font-size:12px;margin:2px 0;",
      icon(if (ok) "check-circle" else "times-circle",
           style = paste0("color:", if (ok) "#27ae60" else "#e74c3c", ";")),
      " ", txt, if (ok) "" else " — falta")
    div(class = paste("alert", if (ok_b && ok_g) "alert-success" else "alert-warning"),
        style = "padding:6px 10px;margin:0;",
      est(ok_b, paste0("BV-BRC", if (ok_b) paste0(" (", rv$bvbrc_user, ")") else "")),
      est(ok_g, "Galaxy"))
  })

  output$wiz_estado <- renderUI({
    ap <- rv$auto_pipe
    if (isTRUE(ap$active))
      return(div(class = "alert alert-info", style = "padding:6px 10px;font-size:12px;margin:0;",
                 icon("spinner", class = "fa-spin"), " Análisis en curso..."))
    n_ok <- sum(rv$completados, na.rm = TRUE)
    if (n_ok > 0)
      return(div(class = "alert alert-success", style = "padding:6px 10px;font-size:12px;margin:0;",
                 icon("check"), paste0(" ", n_ok, "/", length(PASOS_IDS), " etapas completadas")))
    div(style = "font-size:12px;color:#888;", "Listo para iniciar.")
  })

  # Lanzar el pipeline completo desde el asistente
  observeEvent(input$wiz_iniciar, {
    if (length(rv$muestras) == 0) {
      showNotification("Paso 1: selecciona la carpeta con los FASTQ.", type = "warning", duration = 6)
      return()
    }
    if (is.null(rv$bvbrc_token)) {
      showModal(creds_modal())
      showNotification("Paso 2: configura las credenciales de BV-BRC.", type = "warning", duration = 6)
      return()
    }
    if (is.null(rv$galaxy_key))
      showNotification(paste("Sin Galaxy: la taxonomía (Kraken2) y el MLST se omitirán.",
                             "Configura la API key para el análisis completo."),
                       type = "warning", duration = 9)
    updateTextInput(session, "proy_organismo", value = input$wiz_organismo %||% "")
    rv$sesion_nombre <- input$sesion_nombre %||% "sesion_activa"
    alog("INFO", "[ASISTENTE] Iniciando análisis completo",
         paste0("Muestras: ", paste(rv$muestras, collapse = ", "),
                "\nOrganismo: ", input$wiz_organismo %||% ""))
    # Dispara el mismo flujo que el botón del pipeline automático
    pipe_arrancar(input$wiz_organismo %||% "Acinetobacter baumannii")
  }, ignoreInit = TRUE)

  # ============================================================
  # SESIÓN DE TRABAJO (.anbio)
  # ============================================================
  proyecto_meta <- function() list(
    nombre      = isolate(input$proy_nombre)      %||% "",
    organismo   = isolate(input$proy_organismo)   %||% isolate(input$wiz_organismo) %||% "",
    analista    = isolate(input$proy_analista)    %||% "",
    institucion = isolate(input$proy_institucion) %||% "LESP Aguascalientes",
    fecha       = as.character(isolate(input$proy_fecha) %||% Sys.Date()),
    leyenda     = isolate(input$rep_leyenda)      %||% LEYENDA_REPORTE_DEFAULT
  )

  guardar_sesion <- function(nombre = NULL, verbose = TRUE) {
    nm <- nombre %||% isolate(rv$sesion_nombre) %||% "sesion_activa"
    p  <- session_autosave_path(nm)
    ok <- session_save_file(rv, p, proyecto_meta())
    if (ok && verbose)
      alog("OK", paste0("[SESIÓN] Guardada: ", nm),
           paste0(p, "\nJobs: ", length(rv$jobs),
                  " | Muestras: ", length(rv$muestras)))
    ok
  }

  observeEvent(input$btn_sesion_guardar, {
    rv$sesion_nombre <- input$sesion_nombre %||% "sesion_activa"
    if (guardar_sesion(rv$sesion_nombre)) {
      rv$sesion_msg <- paste0("Sesión '", rv$sesion_nombre, "' guardada — ",
                              format(Sys.time(), "%H:%M:%S"))
      showNotification("Sesión guardada.", type = "message", duration = 4)
    } else showNotification("No se pudo guardar la sesión.", type = "error", duration = 6)
  }, ignoreInit = TRUE)

  output$sesion_lista_ui <- renderUI({
    invalidateLater(10000, session)
    df <- session_list()
    if (!nrow(df))
      return(div(style = "font-size:12px;color:#888;padding:6px 0;", "Sin sesiones guardadas."))
    etiquetas <- paste0(df$nombre, "  (", format(df$fecha, "%d/%m %H:%M"), ")")
    selectInput("sesion_sel", NULL, choices = setNames(df$path, etiquetas), width = "100%")
  })

  # Aplica un payload de sesión a rv
  aplicar_sesion <- function(st) {
    rv$muestras     <- st$muestras     %||% rv$muestras
    rv$fastq_dir    <- st$fastq_dir
    rv$fastq_files  <- st$fastq_files  %||% rv$fastq_files
    rv$jobs         <- st$jobs         %||% list()
    rv$job_results  <- st$job_results  %||% list()
    rv$completados  <- st$completados  %||% rv$completados
    rv$auto_pipe    <- st$auto_pipe    %||% rv$auto_pipe
    rv$calidad      <- st$calidad      %||% rv$calidad
    rv$taxonomia    <- st$taxonomia    %||% rv$taxonomia
    rv$ensamblado   <- st$ensamblado   %||% rv$ensamblado
    rv$mlst         <- st$mlst         %||% rv$mlst
    rv$alineamiento <- st$alineamiento %||% rv$alineamiento
    rv$ani          <- st$ani          %||% rv$ani
    rv$resfinder_df <- st$resfinder_df
    rv$amr_sir      <- st$amr_sir
    rv$anotacion    <- st$anotacion
    rv$pangenoma    <- st$pangenoma
    rv$newick       <- st$newick
    pr <- st$proyecto %||% list()
    if (nchar(pr$nombre    %||% "")) updateTextInput(session, "proy_nombre",    value = pr$nombre)
    if (nchar(pr$organismo %||% "")) {
      updateTextInput(session, "proy_organismo", value = pr$organismo)
      updateTextInput(session, "wiz_organismo",  value = pr$organismo)
    }
    if (nchar(pr$analista  %||% "")) updateTextInput(session, "proy_analista",  value = pr$analista)
    if (nchar(pr$leyenda   %||% "")) updateTextAreaInput(session, "rep_leyenda", value = pr$leyenda)
    if (!is.null(st$newick)) updateTextAreaInput(session, "fil_newick", value = st$newick)
    invisible(TRUE)
  }

  # Cargar sesión guardada + refrescar avance real
  cargar_y_reanudar <- function(path) {
    st <- session_load_file(path)
    if (is.null(st)) {
      showNotification("Archivo de sesión inválido.", type = "error", duration = 6); return()
    }
    aplicar_sesion(st)
    rv$sesion_nombre <- tools::file_path_sans_ext(basename(path))
    alog("OK", paste0("[SESIÓN] Cargada: ", basename(path)),
         paste0("Guardada: ", format(st$saved_at, "%d/%m/%Y %H:%M"),
                " | Jobs: ", length(st$jobs %||% list())))
    showModal(modalDialog(
      title = tagList(icon("sync", class = "fa-spin"), " Reanudando análisis..."),
      tags$p("Consultando el avance real en BV-BRC y Galaxy..."),
      footer = NULL, easyClose = FALSE))
    r <- tryCatch(refresh_all_stages(), error = function(e) list(nuevos=0, importados=0, pendientes=0))
    removeModal()
    rv$sesion_msg <- paste0("Sesión reanudada — ", r$importados,
                            " resultado(s) nuevos importados, ",
                            r$pendientes, " aún en proceso.")
    showNotification(rv$sesion_msg, type = "message", duration = 10)
  }

  observeEvent(input$btn_sesion_cargar, {
    req(input$sesion_sel); cargar_y_reanudar(input$sesion_sel)
  }, ignoreInit = TRUE)

  observeEvent(input$up_sesion, {
    req(input$up_sesion$datapath); cargar_y_reanudar(input$up_sesion$datapath)
  }, ignoreInit = TRUE)

  output$dl_sesion <- downloadHandler(
    filename = function() paste0(isolate(rv$sesion_nombre) %||% "sesion", "_",
                                 format(Sys.Date(), "%Y%m%d"), ".anbio"),
    content  = function(file) session_save_file(rv, file, proyecto_meta())
  )

  output$sesion_estado_ui <- renderUI({
    if (is.null(rv$sesion_msg)) return(NULL)
    div(class = "alert alert-success", style = "padding:7px 11px;font-size:12px;margin-top:8px;",
        icon("check-circle"), " ", rv$sesion_msg)
  })

  # ── AUTOGUARDADO al cerrar la app ──
  session$onSessionEnded(function() {
    tryCatch({
      nm <- isolate(rv$sesion_nombre) %||% "sesion_activa"
      session_save_file(rv, session_autosave_path(nm), isolate(proyecto_meta()))
      state_save(rv)
    }, error = function(e) NULL)
  })

  # ============================================================
  # REPORTE FINAL COMPLETO (HTML, 9 etapas)
  # ============================================================
  # Estado real de cada etapa según los datos importados
  etapas_estado <- reactive({
    list(
      "1 · Control de calidad"      = nrow(rv$calidad)    > 0,
      "2 · Clasificación taxonómica"= nrow(rv$taxonomia)  > 0,
      "3 · Ensamblado"              = nrow(rv$ensamblado) > 0,
      "4 · Anotación"               = !is.null(rv$anotacion) && nrow(rv$anotacion) > 0,
      "5 · Resistoma"               = (!is.null(rv$resfinder_df) && nrow(rv$resfinder_df) > 0) ||
                                      (!is.null(rv$amr_sir) && nrow(rv$amr_sir) > 0),
      "6 · Alineamiento"            = nrow(rv$alineamiento) > 0 || !is.null(rv$pangenoma),
      "7 · MLST"                    = nrow(rv$mlst) > 0,
      "8 · Filogenia"               = !is.null(rv$newick) && nchar(trimws(rv$newick %||% "")) > 0,
      "9 · Genómica comparativa"    = !is.null(rv$pangenoma)
    )
  })

  output$etapas_check_ui <- renderUI({
    et <- etapas_estado()
    n_ok <- sum(vapply(et, isTRUE, logical(1)))
    tagList(
      div(style = "display:flex;flex-wrap:wrap;gap:6px;margin:8px 0;",
        lapply(names(et), function(nm) {
          ok <- isTRUE(et[[nm]])
          div(style = paste0("padding:5px 10px;border-radius:14px;font-size:12px;",
                             "background:", if (ok) "#d5f5e3" else "#f4f6f6", ";",
                             "color:", if (ok) "#1e8449" else "#909497", ";",
                             "border:1px solid ", if (ok) "#27ae60" else "#d5d8dc", ";"),
            icon(if (ok) "check-circle" else "circle-notch"), " ", nm)
        })),
      div(style = "font-size:13px;font-weight:bold;color:#1a5276;",
          paste0(n_ok, " de ", length(et), " etapas con resultados")),
      if (n_ok < length(et))
        div(style = "font-size:12px;color:#7d6608;background:#fef9e7;padding:6px 10px;border-radius:4px;margin-top:5px;",
            icon("info-circle"),
            " Puedes generar el reporte con las etapas disponibles; las pendientes se marcarán como tales.")
    )
  })

  observeEvent(input$btn_reporte_final, {
    showModal(modalDialog(title = "Generando reporte final...",
      tags$p(icon("spinner", class = "fa-spin"), " Compilando resultados y gráficas..."),
      footer = NULL, easyClose = FALSE))
    tryCatch({
      dat <- list(
        titulo      = paste0("Reporte de Análisis Bioinformático",
                             if (nchar(input$proy_nombre %||% "")) paste0(" — ", input$proy_nombre) else ""),
        proyecto    = input$proy_nombre      %||% "",
        organismo   = input$proy_organismo   %||% input$wiz_organismo %||% "",
        analista    = input$proy_analista    %||% "",
        institucion = input$proy_institucion %||% "LESP Aguascalientes",
        fecha       = format(Sys.Date(), "%d/%m/%Y"),
        muestras    = rv$muestras,
        etapas      = etapas_estado(),
        calidad     = if (nrow(rv$calidad)      > 0) rv$calidad      else NULL,
        taxonomia   = if (nrow(rv$taxonomia)    > 0) rv$taxonomia    else NULL,
        ensamblado  = if (nrow(rv$ensamblado)   > 0) rv$ensamblado   else NULL,
        anotacion   = rv$anotacion,
        amr_sir     = rv$amr_sir,
        resfinder_df= rv$resfinder_df,
        alineamiento= if (nrow(rv$alineamiento) > 0) rv$alineamiento else NULL,
        mlst        = if (nrow(rv$mlst)         > 0) rv$mlst         else NULL,
        ani         = if (nrow(rv$ani)          > 0) rv$ani          else NULL,
        newick      = rv$newick %||% input$fil_newick,
        renombrar_arbol = mapa_genomas(),
        resaltar_arbol  = c(rv$muestras, unlist(mapa_genomas(), use.names = FALSE)),
        pangenoma   = rv$pangenoma,
        notas       = input$rep_notas_final %||% "",
        leyenda     = input$rep_leyenda %||% LEYENDA_REPORTE_DEFAULT
      )
      out <- file.path(tempdir(),
                       paste0("ANbio_reporte_final_", format(Sys.Date(), "%Y%m%d"), ".html"))
      rmarkdown::render("reporte_final.Rmd", output_file = out,
                        params = list(d = dat), envir = new.env(parent = globalenv()),
                        quiet = TRUE)
      rv$reporte_final_path <- out
      removeModal()
      alog("OK", "Reporte final generado", out)
      showNotification("Reporte final generado. Descárgalo abajo.", type = "message", duration = 8)
    }, error = function(e) {
      removeModal()
      alog("ERROR", "Fallo al generar el reporte final", conditionMessage(e))
      showNotification(paste("Error al generar el reporte:", conditionMessage(e)),
                       type = "error", duration = 12)
    })
  }, ignoreInit = TRUE)

  output$reporte_final_dl <- renderUI({
    req(rv$reporte_final_path)
    if (!file.exists(rv$reporte_final_path)) return(NULL)
    downloadButton("dl_reporte_final", "Descargar reporte HTML",
                   class = "btn-primary", style = "width:100%;")
  })

  output$dl_reporte_final <- downloadHandler(
    filename = function() paste0("ANbio_reporte_final_", format(Sys.Date(), "%Y%m%d"), ".html"),
    content  = function(file) file.copy(rv$reporte_final_path, file, overwrite = TRUE)
  )

  # ─── Log / Consola ────────────────────────────────────────
  LOG_COLORS <- list(
    API   = list(bg="#5b2c6f", text="#e8daef", badge="#8e44ad"),
    OK    = list(bg="#1e4620", text="#a9dfbf", badge="#27ae60"),
    INFO  = list(bg="#1a3a5c", text="#aed6f1", badge="#2980b9"),
    WARN  = list(bg="#4a3000", text="#f9e79f", badge="#e67e22"),
    ERROR = list(bg="#4a1010", text="#f1948a", badge="#e74c3c"),
    DEBUG = list(bg="#2c3e50", text="#aab7b8", badge="#7f8c8d")
  )

  output$log_console <- renderUI({
    entries <- rv$log_entries
    if (length(entries)==0)
      return(div(style="color:#7f8c8d;font-size:12px;padding:8px;", "[ sin entradas de log ]"))

    lapply(rev(entries), function(e) {
      col <- LOG_COLORS[[e$nivel]] %||% LOG_COLORS$INFO
      has_det <- nchar(trimws(e$detalle)) > 0
      div(style=paste0("border-left:3px solid ",col$badge,";background:",col$bg,
                       ";padding:5px 8px 5px 10px;margin-bottom:2px;border-radius:2px;"),
        tags$span(style="color:#566573;font-size:10px;", paste0("[", e$ts, "]")), " ",
        tags$span(style=paste0("background:",col$badge,";color:white;",
                               "padding:1px 5px;border-radius:2px;font-size:10px;font-weight:bold;"), e$nivel), " ",
        tags$span(style=paste0("color:",col$text,";"), e$mensaje),
        if (has_det)
          tags$pre(style=paste0("color:#aab7b8;font-size:10px;margin:3px 0 0 0;",
                                "white-space:pre-wrap;word-break:break-all;",
                                "background:rgba(0,0,0,.2);border-radius:2px;padding:4px 6px;"),
                   e$detalle)
      )
    })
  })

  output$log_stats_ui <- renderUI({
    n <- length(rv$log_entries)
    n_err  <- sum(vapply(rv$log_entries, function(e) isTRUE(e$nivel == "ERROR"), logical(1)))
    n_warn <- sum(vapply(rv$log_entries, function(e) isTRUE(e$nivel == "WARN"),  logical(1)))
    tagList(
      tags$span(style="font-size:12px;color:#888;", paste0(n, " entradas")),
      if (n_err  > 0) tags$span(style="color:#e74c3c;margin-left:6px;font-weight:bold;", paste0(n_err, " error(es)")),
      if (n_warn > 0) tags$span(style="color:#e67e22;margin-left:6px;", paste0(n_warn, " aviso(s)"))
    )
  })

  # ─── Diagnóstico end-to-end BV-BRC (self-test de la cadena API) ──
  # Prueba en vivo: token → mkdir → upload → Workspace.ls → Workspace.get
  # Cada paso se registra en el log y el resultado se muestra como checklist.
  observeEvent(input$btn_diag_bvbrc, {
    tok <- rv$bvbrc_token; user <- rv$bvbrc_user
    if (is.null(tok)) {
      showModal(creds_modal())
      showNotification("Inicia sesión BV-BRC primero.", type = "warning", duration = 5)
      return()
    }

    showModal(modalDialog(
      title = tagList(icon("stethoscope"), " Diagnóstico BV-BRC en curso..."),
      tags$p(icon("spinner", class = "fa-spin"),
             " Probando conexión, carga y descarga de archivos..."),
      tags$small(style = "color:#888;", "Sube un archivo de prueba minúsculo y lo vuelve a leer."),
      footer = NULL, easyClose = FALSE
    ))
    alog("INFO", "[DIAG] Iniciando diagnóstico end-to-end BV-BRC",
         paste0("Usuario: ", user))

    steps <- list()
    add_step <- function(name, ok, detail = "") {
      steps[[length(steps) + 1]] <<- list(name = name, ok = isTRUE(ok),
                                          detail = as.character(detail %||% ""))
    }

    # 1. Token / sesión
    dom <- if (grepl("@", user, fixed = TRUE))
             sub("^[^@]+", "", user)
           else "(sin dominio en el token — se usará @bvbrc)"
    add_step("Sesión / token activo", TRUE,
             paste0("Usuario: ", user, "  |  Dominio workspace: ", dom))
    alog("OK", "[DIAG] 1/6 Token presente", paste0("Usuario: ", user, " | Dominio: ", dom))

    # 2. Crear carpeta de salida
    mk <- tryCatch({ bvbrc_mkdir(tok, ws_out(user)); TRUE }, error = function(e) FALSE)
    add_step("Crear carpeta de resultados (Workspace.create)", mk, ws_out(user))
    alog(if (mk) "OK" else "WARN", "[DIAG] 2/6 mkdir carpeta de salida", ws_out(user))

    # 3. Subir archivo de prueba (protocolo Shock)
    tf <- file.path(tempdir(), "anbio_diag_test.txt")
    writeLines(paste0("ANbio diagnostico ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), tf)
    up <- tryCatch(bvbrc_upload_file(tok, user, tf, "anbio_diag_test.txt"),
                   error = function(e) list(success = FALSE, error = conditionMessage(e)))
    add_step("Subir archivo (Shock upload)", isTRUE(up$success),
             up$log_det %||% up$error %||% "")
    alog(if (isTRUE(up$success)) "OK" else "ERROR", "[DIAG] 3/6 Upload de prueba",
         up$log_det %||% up$error %||% "sin detalle")

    # 4. Listar el workspace (Workspace.ls) — el punto históricamente frágil
    up_dir <- paste0(ws_home(user), "/.fastq_uploads")
    ls1 <- tryCatch(bvbrc_ls(tok, up_dir),
                    error = function(e) list(success = FALSE, error = conditionMessage(e)))
    found <- isTRUE(ls1$success) &&
      any(vapply(ls1$files %||% list(),
                 function(f) grepl("anbio_diag_test", f$name %||% ""), logical(1)))
    add_step("Listar workspace (Workspace.ls)", isTRUE(ls1$success),
             if (isTRUE(ls1$success))
               paste0(length(ls1$files), " archivo(s) en ", up_dir,
                      "  |  archivo de prueba hallado: ", if (found) "SÍ" else "no")
             else paste0("ERROR: ", ls1$error %||% "desconocido"))
    alog(if (isTRUE(ls1$success)) "OK" else "ERROR", "[DIAG] 4/6 Workspace.ls",
         if (isTRUE(ls1$success)) paste0(length(ls1$files), " archivo(s); test: ", found)
         else ls1$error %||% "")

    # 5. Descargar el archivo (Workspace.get)
    test_path <- up$ws_path %||% paste0(up_dir, "/anbio_diag_test.txt")
    g <- tryCatch(bvbrc_get_file(tok, test_path),
                  error = function(e) list(success = FALSE, error = conditionMessage(e)))
    g_ok <- isTRUE(g$success) && nchar(trimws(g$content %||% "")) > 0
    add_step("Descargar archivo (Workspace.get)", g_ok,
             if (g_ok) paste0("Contenido leído: ", substr(trimws(g$content), 1, 80))
             else paste0("ERROR: ", g$error %||% "sin contenido"))
    alog(if (g_ok) "OK" else "ERROR", "[DIAG] 5/6 Workspace.get",
         if (g_ok) "archivo recuperado correctamente" else g$error %||% "")

    # 6. Listar la carpeta de resultados (donde caen los outputs de los jobs)
    ls2 <- tryCatch(bvbrc_ls(tok, ws_out(user)),
                    error = function(e) list(success = FALSE, error = conditionMessage(e)))
    add_step("Listar carpeta de resultados", isTRUE(ls2$success),
             if (isTRUE(ls2$success))
               paste0(length(ls2$files), " item(s) en ", ws_out(user))
             else paste0("ERROR: ", ls2$error %||% "desconocido"))
    alog(if (isTRUE(ls2$success)) "OK" else "WARN", "[DIAG] 6/6 ls carpeta de resultados",
         if (isTRUE(ls2$success)) paste0(length(ls2$files), " item(s)") else ls2$error %||% "")

    rv$diag_result <- list(steps = steps, ts = format(Sys.time(), "%H:%M:%S"))
    removeModal()
    n_ok <- sum(vapply(steps, function(s) isTRUE(s$ok), logical(1)))
    alog(if (n_ok == length(steps)) "OK" else "WARN",
         paste0("[DIAG] Diagnóstico finalizado: ", n_ok, "/", length(steps), " pasos OK"))
    showNotification(
      paste0("Diagnóstico BV-BRC: ", n_ok, "/", length(steps),
             " pasos OK. Revisa el resumen y el log."),
      type = if (n_ok == length(steps)) "message" else "warning", duration = 8)
  }, ignoreInit = TRUE)

  output$diag_bvbrc_result <- renderUI({
    dr <- rv$diag_result
    if (is.null(dr)) return(NULL)
    n_ok <- sum(vapply(dr$steps, function(s) isTRUE(s$ok), logical(1)))
    n    <- length(dr$steps)
    head_cls <- if (n_ok == n) "#1e4620" else "#4a3000"
    div(style = paste0("background:", head_cls, ";border-radius:6px;padding:10px 14px;",
                       "margin-bottom:12px;border:1px solid #34495e;"),
      div(style = "color:#ecf0f1;font-weight:bold;font-size:13px;margin-bottom:6px;",
        icon("stethoscope"),
        paste0(" Diagnóstico BV-BRC — ", n_ok, "/", n, " pasos OK  (", dr$ts, ")")),
      lapply(dr$steps, function(s) {
        col <- if (isTRUE(s$ok)) "#2ecc71" else "#e74c3c"
        ic  <- if (isTRUE(s$ok)) "check-circle" else "times-circle"
        div(style = "margin:3px 0;font-size:12px;color:#dfe6e9;",
          tags$span(style = paste0("color:", col, ";font-weight:bold;"),
                    icon(ic), " ", s$name),
          if (nchar(trimws(s$detail)) > 0)
            tags$pre(style = paste0("color:#aab7b8;font-size:10px;margin:2px 0 0 22px;",
                                    "white-space:pre-wrap;word-break:break-all;",
                                    "background:rgba(0,0,0,.25);border-radius:2px;padding:3px 6px;"),
                     s$detail)
        )
      })
    )
  })

  # ─── Helpers de resultados BV-BRC ────────────────────────

  # Fallback cuando Workspace.ls falla: intenta descargar archivos por rutas conocidas
  bvbrc_fetch_by_known_paths <- function(tok, job, out_dir) {
    m  <- job$muestra %||% ""
    # Rutas conocidas por herramienta (sufijos comunes de BV-BRC)
    known <- switch(job$tool,
      FastqUtils = list(
        paste0(out_dir, "/multiqc_data/multiqc_general_stats.txt"),
        paste0(out_dir, "/multiqc_data/multiqc_fastqc.txt"),
        paste0(out_dir, "/multiqc_data/multiqc_data.json")
      ),
      TaxonomicClassification = list(
        paste0(out_dir, "/", m, ".report"),
        paste0(out_dir, "/TaxClassificationReport.txt"),
        paste0(out_dir, "/", m, "_classification.txt")
      ),
      Assembly2 = list(
        paste0(out_dir, "/quast_results/report.tsv"),
        paste0(out_dir, "/report.txt"),
        paste0(out_dir, "/assembly_info.txt")
      ),
      MetagenomicReadMapping = list(
        paste0(out_dir, "/", m, ".CARD.json"),
        paste0(out_dir, "/", m, ".tsv"),
        paste0(out_dir, "/CARD_kmer_data.json")
      ),
      list()
    )
    contents <- list()
    for (p in known) {
      g <- tryCatch(bvbrc_get_file(tok, p), error = function(e) list(success=FALSE))
      if (isTRUE(g$success) && nchar(trimws(g$content %||% "")) > 0)
        contents[[basename(p)]] <- g$content
    }
    if (!length(contents))
      return(list(success = FALSE,
                  error = paste0("Workspace.ls falló y Workspace.get no encontró archivos en ", out_dir)))
    pick <- function(pats) {
      idx <- grep(paste(pats, collapse="|"), names(contents), ignore.case=TRUE)
      if (length(idx)) contents[[idx[1]]] else (contents[[1]] %||% "")
    }
    parsed <- switch(job$tool,
      FastqUtils              = list(type="calidad",   data=parse_fastqutils(contents)),
      TaxonomicClassification = list(type="taxonomia", data=parse_taxonomy(pick(c(".report","classif")))),
      Assembly2               = list(type="ensamblado",data=parse_assembly(pick(c("report","quast")))),
      MetagenomicReadMapping  = list(type="resistoma", data=parse_resistome(pick(c("CARD","amr","tsv")))),
      list(type="raw", data=NULL)
    )
    list(success=TRUE, tool=job$tool, muestra=m,
         out_dir=out_dir, files=list(), contents=contents, parsed=parsed)
  }

  bvbrc_fetch_results <- function(tok, job) {
    out_dir <- paste0(job$output_path, "/", job$output_file)

    patterns <- switch(job$tool,
      FastqUtils              = c("fastqc","general_stats",".html",".txt",".tsv"),
      TaxonomicClassification = c(".report","summary","classification",".tsv",".txt"),
      Assembly2               = c("contigs\\.fasta$","report","quast","assembly_info",".txt",".tsv"),
      Annotation              = c("\\.features\\.txt$","genome_quality_details","amr-sir",
                                  "genome_stat","summary"),
      MetagenomicReadMapping  = c("kma\\.res$","CARD","amr","resistance",".tsv",".txt"),
      SeqComparison           = c("genome_comparison\\.txt$"),
      c(".txt",".tsv",".json")
    )
    want     <- function(nm) any(vapply(patterns, function(p)
                  grepl(p, nm, ignore.case = TRUE), logical(1)))
    skip_big <- function(nm) grepl("\\.(fastq|fq|bam|sam|bai|cram|gz)$", nm, ignore.case = TRUE)

    contents <- list()

    # 1) Fuente fiable: el objeto job_result trae la lista 'output_files' con las
    #    rutas reales (en la carpeta oculta .<output_file>). Leerlas y descargar.
    jr <- tryCatch(bvbrc_get_file(tok, out_dir), error = function(e) list(success = FALSE))
    if (isTRUE(jr$success) && nchar(trimws(jr$content %||% "")) > 0) {
      jj <- tryCatch(jsonlite::fromJSON(jr$content, simplifyVector = FALSE),
                     error = function(e) NULL)
      of <- jj$output_files %||% list()
      for (item in of) {
        p  <- as.character(if (is.list(item)) item[[1]] else item)
        nm <- basename(p)
        if (skip_big(nm) || !want(nm)) next
        g <- tryCatch(bvbrc_get_file(tok, p), error = function(e) list(success = FALSE))
        if (isTRUE(g$success) && nchar(trimws(g$content %||% "")) > 0)
          contents[[nm]] <- g$content
      }
    }

    # 2) Fallback: listar la carpeta oculta .<output_file> (y la visible) + subcarpetas
    if (!length(contents)) {
      hidden    <- paste0(job$output_path, "/.", job$output_file)
      all_files <- list()
      for (dir_try in c(hidden, out_dir)) {
        lsx <- tryCatch(bvbrc_ls(tok, dir_try), error = function(e) list(success = FALSE))
        if (isTRUE(lsx$success)) {
          all_files <- c(all_files, lsx$files)
          for (f in lsx$files) {
            if ((f$type %||% "") %in% c("folder","directory","job_result")) {
              ls2 <- tryCatch(bvbrc_ls(tok, f$path), error = function(e) list(success = FALSE))
              if (isTRUE(ls2$success)) all_files <- c(all_files, ls2$files)
            }
          }
        }
      }
      for (f in all_files) {
        if ((f$type %||% "") %in% c("folder","directory")) next
        if (!is.na(f$size %||% NA) && f$size > 5e6)        next
        if (skip_big(f$name) || !want(f$name))             next
        g <- tryCatch(bvbrc_get_file(tok, f$path), error = function(e) list(success = FALSE))
        if (isTRUE(g$success) && nchar(trimws(g$content %||% "")) > 0)
          contents[[f$name]] <- g$content
      }
    }

    # 3) Último recurso: rutas de archivo conocidas
    if (!length(contents)) {
      kp <- bvbrc_fetch_by_known_paths(tok, job, out_dir)
      if (isTRUE(kp$success)) return(kp)
    }

    pick_content <- function(patterns_pick) {
      if (!length(contents)) return("")   # evita "subíndice fuera de los límites"
      idx <- grep(paste(patterns_pick, collapse="|"), names(contents), ignore.case=TRUE)
      if (length(idx)) contents[[idx[1]]] else contents[[1]]
    }

    parsed <- switch(job$tool,
      FastqUtils = list(type="calidad",
        data = parse_fastqutils(contents)),
      TaxonomicClassification = list(type="taxonomia",
        data = parse_taxonomy(pick_content(c(".report","summary","classif")))),
      Assembly2 = list(type="ensamblado",
        data = parse_assembly(pick_content(c("contigs\\.fasta$","report","quast")))),
      Annotation = list(type="anotacion",
        data = parse_annotation(contents)),
      MetagenomicReadMapping = list(type="resistoma",
        data = parse_resistome(pick_content(c("kma\\.res$","CARD","amr","resistance")))),
      CoreGenomeMLST = list(type="mlst_galaxy",
        data = parse_mlst_tsv(pick_content(c("mlst","allele","ST",".tsv",".txt")))),
      PhylogeneticTree = list(type="newick",
        data = parse_newick(pick_content(c(".nwk",".tree","newick","phylo")))),
      SeqComparison = list(type="genomica",
        data = parse_genome_comparison(pick_content(c("genome_comparison\\.txt$")))),
      list(type="raw", data=NULL)
    )

    list(success=TRUE, tool=job$tool, muestra=job$muestra %||% "",
         out_dir=out_dir, files=as.list(names(contents)), contents=contents, parsed=parsed)
  }

  # ─── Helpers de resultados Galaxy ───────────────────────
  galaxy_fetch_results <- function(key, job) {
    # 1. Obtener IDs de output (guardados al enviar o consultados ahora)
    ids <- job$output_ids
    if (is.null(ids) || !length(ids)) {
      outs <- galaxy_job_outputs(key, job$job_id)
      if (!outs$success) return(list(success = FALSE, error = outs$error))
      ids <- outs$ids
    }
    if (!length(ids)) return(list(success = FALSE, error = "Job sin outputs"))

    # 2. Descargar cada dataset (solo si está ok y < 10 MB)
    contents <- list()
    for (did in ids) {
      m <- galaxy_dataset_meta(key, did)
      if (!m$success || m$meta$state != "ok") next
      fname <- m$meta$name %||% did
      fsize <- m$meta$file_size %||% 0
      if (fsize > 10e6) {
        contents[[fname]] <- paste0("[", round(fsize/1e6,1), " MB — disponible en Galaxy]")
        next
      }
      dl <- tryCatch(
        galaxy_req(paste0("/api/datasets/", did, "/display"), key),
        error = \(e) NULL)
      if (!is.null(dl) && httr2::resp_status(dl) == 200)
        contents[[fname]] <- httr2::resp_body_string(dl)
    }

    # 3. Parsear según herramienta
    parsed <- switch(job$tool,
      spades = , unicycler = {
        fk <- names(contents)[grepl("scaffold|contig|assembly|fasta",
                                    names(contents), ignore.case = TRUE)][1]
        txt <- if (!is.na(fk %||% NA)) contents[[fk]] else (contents[[1]] %||% "")
        if (nchar(trimws(txt)) > 0 && !grepl("MB —", txt)) {
          n_ctg  <- length(gregexpr(">", txt)[[1]])
          bp     <- nchar(gsub(">\\S+[^\n]*\n|[^ACGTNacgtn]", "", txt))
          list(type = "ensamblado", data = list(
            contigs  = n_ctg,
            n50_kb   = NA,
            total_mb = round(bp / 1e6, 2),
            gc       = NA, largest = NA))
        } else list(type = "raw", data = NULL)
      },
      mlst = {
        tk <- names(contents)[grepl("\\.tsv$|mlst|results",
                                    names(contents), ignore.case = TRUE)][1]
        if (!is.na(tk %||% NA))
          list(type = "mlst_galaxy", data = contents[[tk]])
        else list(type = "raw", data = NULL)
      },
      iqtree = {
        nk <- names(contents)[grepl("\\.treefile|\\.nwk|tree",
                                    names(contents), ignore.case = TRUE)][1]
        if (!is.na(nk %||% NA))
          list(type = "newick", data = contents[[nk]])
        else list(type = "raw", data = NULL)
      },
      kraken2 = {
        # El dataset "Report" trae el formato kraken2 report (tabulado),
        # que parse_taxonomy() interpreta igual que el de BV-BRC.
        rk <- names(contents)[grepl("report", names(contents), ignore.case = TRUE)][1]
        txt <- if (!is.na(rk %||% NA)) contents[[rk]] else (contents[[1]] %||% "")
        if (nchar(trimws(txt)) > 0 && !grepl("MB —", txt))
          list(type = "taxonomia", data = parse_taxonomy(txt))
        else list(type = "raw", data = NULL)
      },
      list(type = "raw", data = NULL)
    )

    list(success = TRUE, tool = job$tool, muestra = job$muestra %||% "",
         output_ids = ids, contents = contents, parsed = parsed)
  }

  auto_import_results <- function(res, job) {
    parsed  <- res$parsed
    muestra <- res$muestra %||% job$muestra %||% ""

    if (parsed$type == "calidad" && is.list(parsed$data) && length(parsed$data) > 0) {
      for (sdata in parsed$data) {
        m <- sdata$sample %||% muestra
        if (!nchar(m %||% "")) next
        estado <- if (!is.na(sdata$q30 %||% NA) && (sdata$q30 %||% 0) >= 75) "APROBADO" else "REVISAR"
        row <- data.frame(Muestra=m, Q30=sdata$q30 %||% NA, GC=sdata$gc %||% NA,
          Reads_M=sdata$reads_m %||% NA, Duplicados=sdata$dup %||% NA,
          Estado=estado, stringsAsFactors=FALSE)
        rv$calidad <- rbind(rv$calidad[rv$calidad$Muestra != m, ], row)
      }
    } else if (parsed$type == "taxonomia" && !is.null(parsed$data)) {
      d <- parsed$data; m <- muestra
      if (!nchar(m %||% "")) return(invisible(NULL))
      # Confianza: se prioriza el % a nivel de GÉNERO. Kraken2 deja la mayoría
      # de lecturas en el género cuando la especie no es discriminable
      # (complejo ACB), por lo que el % de especie subestima la identificación.
      conf <- d$pct_genero %||% d$confianza %||% NA
      row <- data.frame(Muestra=m, Especie=d$top_especie %||% "",
        Confianza=conf, Contaminacion=d$contaminacion %||% NA,
        stringsAsFactors=FALSE)
      rv$taxonomia <- rbind(rv$taxonomia[rv$taxonomia$Muestra != m, ], row)
    } else if (parsed$type == "ensamblado" && !is.null(parsed$data)) {
      d <- parsed$data; m <- muestra
      if (!nchar(m %||% "")) return(invisible(NULL))
      row <- data.frame(Muestra=m, Cobertura=NA,
        Profundidad=d$profundidad %||% NA,
        Contigs=d$contigs %||% NA, N50_kb=d$n50_kb %||% NA,
        Tamano_Mb=d$total_mb %||% NA, GC=d$gc %||% NA, stringsAsFactors=FALSE)
      rv$ensamblado <- rbind(rv$ensamblado[rv$ensamblado$Muestra != m, ], row)
    } else if (parsed$type == "anotacion" && !is.null(parsed$data)) {
      d <- parsed$data
      if (!is.na(d$cds %||% NA)) {
        updateNumericInput(session, "anot_cds",     value = d$cds)
        updateNumericInput(session, "anot_rrna",    value = d$rrna %||% NA)
        updateNumericInput(session, "anot_trna",    value = d$trna %||% NA)
        updateNumericInput(session, "anot_funcion", value = d$funcion_pct %||% NA)
      }
      # Acumular por muestra para el reporte final
      if (nchar(muestra %||% "")) {
        fila <- data.frame(Muestra = muestra, CDS = d$cds %||% NA,
          rRNA = d$rrna %||% NA, tRNA = d$trna %||% NA,
          Hipoteticos = d$hypo %||% NA, Pct_funcion = d$funcion_pct %||% NA,
          Completitud = d$completitud %||% NA,
          Contaminacion = d$contaminacion %||% NA, stringsAsFactors = FALSE)
        prev <- rv$anotacion
        rv$anotacion <- if (is.null(prev)) fila
                        else rbind(prev[prev$Muestra != muestra, , drop = FALSE], fila)
      }
      # Antibiograma predicho (viene en el mismo job de anotación)
      if (!is.null(res$contents)) {
        sir_txt <- res$contents[["amr-sir.txt"]] %||% NULL
        if (!is.null(sir_txt)) {
          s <- tryCatch(parse_amr_sir(sir_txt), error = function(e) NULL)
          if (!is.null(s) && nrow(s) > 0) {
            s$Muestra <- muestra
            prev <- rv$amr_sir
            rv$amr_sir <- if (is.null(prev)) s
                          else rbind(prev[prev$Muestra != muestra, , drop = FALSE], s)
          }
        }
      }
    } else if (parsed$type == "resistoma" && !is.null(parsed$data) &&
               is.data.frame(parsed$data) && nrow(parsed$data) > 0) {
      rv$resfinder_df <- parsed$data

    } else if (parsed$type == "mlst_galaxy" && !is.null(parsed$data)) {
      # parsed$data puede ser data.frame (BV-BRC) o texto TSV (Galaxy)
      df <- if (is.data.frame(parsed$data)) {
        parsed$data
      } else {
        parse_mlst_tsv(as.character(parsed$data))
      }
      if (!is.null(df) && nrow(df) > 0) {
        for (i in seq_len(nrow(df))) {
          m <- df$Muestra[i] %||% muestra
          rv$mlst <- rbind(rv$mlst[rv$mlst$Muestra != m, ], df[i, , drop = FALSE])
        }
      }

    } else if (parsed$type == "newick" && !is.null(parsed$data) &&
               nchar(trimws(parsed$data)) > 0) {
      updateTextAreaInput(session, "fil_newick", value = parsed$data)
      rv$newick <- parsed$data

    } else if (parsed$type == "genomica" && !is.null(parsed$data)) {
      d <- parsed$data
      updateNumericInput(session, "gc_core", value = d$core       %||% NA)
      updateNumericInput(session, "gc_acc",  value = d$accesorios %||% NA)
      updateNumericInput(session, "gc_uniq", value = d$unicos     %||% NA)
      rv$pangenoma <- d
    }

    # Marcar el paso correspondiente como completado en la barra de progreso
    step_of_type <- c(calidad = "calidad", taxonomia = "taxonomia",
                      ensamblado = "ensamblado", anotacion = "anotacion",
                      resistoma = "resistoma", mlst_galaxy = "mlst",
                      newick = "filogenia")
    pid_done <- step_of_type[[parsed$type %||% ""]]
    if (!is.null(pid_done) && pid_done %in% names(rv$completados)) {
      rv$completados[pid_done] <- TRUE
      updateCheckboxInput(session, paste0("chk_", pid_done), value = TRUE)
    }

    tryCatch(state_save(rv), error = function(e) NULL)
    invisible(NULL)
  }

  # ─── Modal de resultados ──────────────────────────────────
  observeEvent(input$view_result_job, {
    req(input$view_result_job)
    jid <- input$view_result_job
    rv$viewing_job_id <- jid

    job_list <- Filter(function(j) j$job_id == jid, rv$jobs)
    if (!length(job_list)) return()
    job <- job_list[[1]]
    res <- rv$job_results[[jid]]

    showModal(modalDialog(
      title     = tagList(icon("chart-bar"), " Resultados — ", job$label),
      size      = "l",
      easyClose = TRUE,
      div(style="font-size:11px;color:#888;margin-bottom:6px;",
        icon("folder-open"), " ",
        if (!is.null(job$output_path) && !is.null(job$output_file))
          code(paste0(job$output_path, "/", job$output_file))
        else
          "Directorio de salida no disponible"
      ),
      if (is.null(res)) {
        div(class="alert alert-info",
          icon("info-circle"),
          " Los resultados aún se están descargando o el job acaba de completar. ",
          tags$button(
            class   = "btn btn-sm btn-primary",
            style   = "margin-left:8px;",
            onclick = paste0("Shiny.setInputValue('fetch_result_job','",
                             jid, "',{priority:'event'})"),
            HTML('<i class="fa fa-download"></i> Descargar ahora')
          )
        )
      } else {
        tagList(
          if (length(res$files) > 0)
            tags$details(
              style = "margin-bottom:8px;",
              tags$summary(style="cursor:pointer;color:#2980b9;font-size:12px;",
                icon("folder"), paste0(" Archivos en workspace (", length(res$files), ")")),
              tags$ul(style="font-size:11px;max-height:90px;overflow-y:auto;margin:4px 0 0;",
                lapply(head(res$files, 25), function(f)
                  tags$li(code(f$name),
                    tags$small(style="color:#999;",
                      if (!is.na(f$size %||% NA) && f$size > 0)
                        paste0(" — ", round(f$size/1024), " KB")
                      else "")))
              )
            ),
          hr(style="margin:8px 0;"),
          uiOutput("res_modal_content")
        )
      },
      footer = tagList(
        if (!is.null(res))
          tags$a(class = "btn btn-default btn-sm",
            href    = paste0("https://bv-brc.org/workspace", res$out_dir),
            target  = "_blank",
            icon("external-link-alt"), " Ver en BV-BRC"),
        modalButton("Cerrar")
      )
    ))
  })

  # Descarga manual de resultados
  observeEvent(input$fetch_result_job, {
    req(input$fetch_result_job)
    jid <- input$fetch_result_job
    tok <- rv$bvbrc_token
    job_list <- Filter(function(j) j$job_id == jid, rv$jobs)
    if (!length(job_list) || is.null(tok)) return()
    job <- job_list[[1]]
    if (is.null(job$output_path) || is.null(job$output_file)) {
      showNotification("output_path no disponible.", type="warning"); return()
    }
    showModal(modalDialog(
      title  = "Descargando resultados...",
      tags$p(icon("spinner", class="fa-spin"), " Listando y descargando archivos del workspace..."),
      footer = NULL, easyClose = FALSE
    ))
    fetch_res <- tryCatch(bvbrc_fetch_results(tok, job),
      error = function(e) list(success=FALSE, error=conditionMessage(e)))
    removeModal()
    if (fetch_res$success) {
      rv$job_results[[jid]] <- fetch_res
      auto_import_results(fetch_res, job)
      rv$viewing_job_id <- jid
      alog("OK", paste0("Resultados descargados manualmente: ", job$label),
           paste0("Archivos: ", length(fetch_res$files)))
      showNotification("Resultados descargados e importados a las tablas.", type="message", duration=5)
      # Reabrir modal con contenido
      showModal(modalDialog(
        title     = tagList(icon("chart-bar"), " Resultados — ", job$label),
        size      = "l", easyClose = TRUE,
        if (length(fetch_res$files) > 0)
          tags$details(
            style="margin-bottom:8px;",
            tags$summary(style="cursor:pointer;color:#2980b9;font-size:12px;",
              icon("folder"), paste0(" Archivos (", length(fetch_res$files), ")")),
            tags$ul(style="font-size:11px;max-height:90px;overflow-y:auto;",
              lapply(head(fetch_res$files,25), function(f)
                tags$li(code(f$name))))
          ),
        hr(style="margin:8px 0;"),
        uiOutput("res_modal_content"),
        footer = tagList(
          tags$a(class="btn btn-default btn-sm",
            href=paste0("https://bv-brc.org/workspace", fetch_res$out_dir),
            target="_blank", icon("external-link-alt"), " Ver en BV-BRC"),
          modalButton("Cerrar"))
      ))
    } else {
      alog("WARN", paste0("Fetch manual fallido: ", job$label), fetch_res$error)
      showNotification(paste("Error:", fetch_res$error), type="error", duration=8)
    }
  })

  # Importación manual desde modal
  observeEvent(input$import_job_result, {
    req(input$import_job_result)
    info <- input$import_job_result
    jid  <- info$job_id
    res  <- rv$job_results[[jid]]
    if (is.null(res)) return()
    job_list <- Filter(function(j) j$job_id == jid, rv$jobs)
    job <- if (length(job_list)) job_list[[1]] else list(job_id=jid, muestra=res$muestra%||%"")
    auto_import_results(res, job)
    alog("INFO", paste0("Datos reimportados: ", res$parsed$type), jid)
    showNotification(paste0("Datos importados a tabla de ", res$parsed$type), type="message", duration=4)
  })

  # Contenido reactivo del modal de resultados
  output$res_modal_content <- renderUI({
    jid <- rv$viewing_job_id; req(jid)
    res <- rv$job_results[[jid]]
    if (is.null(res)) return(div(class="alert alert-warning", "Sin resultados disponibles."))

    parsed  <- res$parsed
    muestra <- res$muestra %||% ""

    tbl_html <- function(df, max_rows = 15) {
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
      df <- head(df, max_rows)
      df[] <- lapply(df, function(x) if (is.numeric(x)) round(x, 2) else x)
      HTML(paste0(
        '<div style="overflow-x:auto;max-height:200px;overflow-y:auto;margin-top:6px;">',
        '<table class="table table-condensed table-striped" style="font-size:12px;margin:0;">',
        '<thead><tr>',
        paste0('<th style="white-space:nowrap;background:#ecf0f1;">',
               names(df), '</th>', collapse=''),
        '</tr></thead><tbody>',
        paste0(apply(df, 1, function(row)
          paste0('<tr>', paste0('<td>', row, '</td>', collapse=''), '</tr>')),
          collapse=''),
        '</tbody></table></div>'
      ))
    }

    import_btn <- function(type, label, cls = "btn-primary") {
      tags$button(
        class   = paste("btn btn-sm", cls),
        style   = "margin-top:10px;",
        onclick = paste0(
          "Shiny.setInputValue('import_job_result',",
          '{\"job_id\":\"', jid, '\",\"type\":\"', type, '\"}',
          ",{priority:'event'})"),
        HTML(paste0('<i class="fa fa-table"></i> ', label))
      )
    }

    switch(parsed$type,

      calidad = {
        rows <- parsed$data
        if (is.null(rows) || !length(rows))
          return(div(class="alert alert-warning",
            icon("exclamation-triangle"),
            " No se encontraron estadísticas de calidad. ",
            "Los archivos FastQC/MultiQC pueden estar aún procesándose."))
        df <- do.call(rbind, lapply(rows, function(r)
          data.frame(Muestra=r$sample%||%"", `Q30 (%)`=r$q30%||%NA,
            `GC (%)`=r$gc%||%NA, `Reads (M)`=r$reads_m%||%NA,
            `Dup (%)`=r$dup%||%NA, check.names=FALSE, stringsAsFactors=FALSE)))
        tagList(
          div(class="alert alert-success", style="padding:6px 10px;",
            icon("check-circle"), strong(" FastQC / MultiQC — Estadísticas por muestra")),
          tbl_html(df),
          import_btn("calidad", "Importar a tabla de Calidad y actualizar gráficas")
        )
      },

      taxonomia = {
        d <- parsed$data
        if (is.null(d))
          return(div(class="alert alert-warning",
            "No se pudo parsear el reporte taxonómico (Kraken2)."))
        top5_df <- if (length(d$top5) > 0)
          data.frame(
            Especie    = sapply(d$top5, `[[`, "name"),
            `% Reads`  = sapply(d$top5, `[[`, "pct"),
            check.names = FALSE, stringsAsFactors = FALSE)
        else NULL
        tagList(
          div(class="alert alert-info", style="padding:8px 12px;",
            icon("leaf"), strong(" Clasificación taxonómica (Kraken2)"), br(),
            tags$table(style="margin:6px 0 0;font-size:13px;width:100%;",
              tags$tr(
                tags$td(style="width:160px;", tags$b("Especie principal:")),
                tags$td(strong(d$top_especie))),
              tags$tr(
                tags$td(tags$b("Confianza:")),
                tags$td(paste0(d$confianza, " %"))),
              tags$tr(
                tags$td(tags$b("% No clasificado:")),
                tags$td(style=if(!is.na(d$contaminacion%||%NA)&&d$contaminacion>5)
                                "color:#e74c3c;font-weight:bold;" else "",
                  paste0(d$contaminacion, " %")))
            )
          ),
          if (!is.null(top5_df)) tbl_html(top5_df),
          import_btn("taxonomia", "Importar a tabla Taxonómica y actualizar gráficas")
        )
      },

      ensamblado = {
        d <- parsed$data
        if (is.null(d))
          return(div(class="alert alert-warning", "No se pudo parsear el reporte QUAST."))
        df <- data.frame(
          Métrica = c("N° Contigs","N50 (kb)","Tamaño total (Mb)","GC (%)","Contig mayor (bp)"),
          Valor   = c(d$contigs%||%"?", d$n50_kb%||%"?", d$total_mb%||%"?",
                      d$gc%||%"?", d$largest%||%"?"),
          stringsAsFactors=FALSE)
        tagList(
          div(class="alert alert-info", style="padding:6px 10px;",
            icon("puzzle-piece"), strong(" Ensamblado — QUAST")),
          tbl_html(df),
          import_btn("ensamblado", "Importar a tabla de Ensamblado y actualizar gráficas")
        )
      },

      anotacion = {
        d <- parsed$data
        if (is.null(d) || all(is.na(unlist(d[c("cds","rrna","trna")]))))
          return(div(class="alert alert-warning",
            "No se encontraron estadísticas de anotación parseables."))
        df <- data.frame(
          Campo   = c("Total CDS","ARNr","ARNt","% Funcionales"),
          Valor   = c(d$cds%||%"?", d$rrna%||%"?", d$trna%||%"?",
                      if (!is.na(d$funcion_pct%||%NA)) paste0(d$funcion_pct,"%") else "?"),
          stringsAsFactors = FALSE)
        tagList(
          div(class="alert alert-info", style="padding:6px 10px;",
            icon("tags"), strong(" Anotación RAST")),
          tbl_html(df),
          import_btn("anotacion", "Importar campos de Anotación")
        )
      },

      resistoma = {
        d <- parsed$data
        if (is.null(d) || !is.data.frame(d) || nrow(d) == 0)
          return(div(class="alert alert-warning",
            "No se detectaron genes AMR o el archivo de resultados no pudo parsearse."))
        tagList(
          div(class="alert alert-danger", style="padding:6px 10px;",
            icon("shield-alt"),
            strong(paste0(" ", nrow(d), " genes AMR detectados"))),
          tbl_html(d),
          import_btn("resistoma", "Importar genes AMR a tabla de Resistoma", "btn-danger")
        )
      },

      # Herramienta sin parser específico: muestra primeros 500 chars del primer archivo
      {
        div(
          div(class="alert alert-info", style="padding:6px 10px;",
            icon("info-circle"),
            " Job completado. No hay parser específico para esta herramienta. Primer archivo:"),
          if (length(res$contents) > 0)
            tags$pre(style="font-size:10px;max-height:180px;overflow-y:auto;
                            background:#f8f9fa;padding:8px;border-radius:4px;",
              substr(res$contents[[1]], 1, 800))
        )
      }
    )
  })

  observeEvent(input$btn_clear_log, {
    rv$log_entries <- list()
    log_clear()
    showNotification("Log limpiado.", type="message", duration=2)
  })

  observeEvent(input$btn_refresh_log, {
    alog("DEBUG", "Log refrescado manualmente")
  })

  output$dl_log_file <- downloadHandler(
    filename = function() paste0("anbio_log_", format(Sys.Date(), "%Y%m%d_%H%M%S"), ".txt"),
    content  = function(file) {
      txt <- log_tail(2000)
      writeLines(txt, file)
    }
  )

  # ============================================================
  # ANÁLISIS POR LOTES
  # ============================================================

  # Mapa herramienta -> step_id
  BATCH_STEP_MAP <- c(
    FastqUtils                  = "calidad",
    TaxonomicClassification     = "taxonomia",
    Assembly2                   = "ensamblado",
    MetagenomicReadMapping      = "resistoma",
    ComprehensiveGenomeAnalysis = "anotacion"
  )

  # ── Indicador reactivo de volumen ────────────────────────
  output$batch_status_ui <- renderUI({
    n_m <- length(rv$muestras)
    tools_sel <- input$batch_tools %||% character(0)
    use_cga   <- isTRUE(input$batch_cga)
    if (use_cga) {
      tools_sel <- c(tools_sel[!tools_sel %in% c("Assembly2","MetagenomicReadMapping")],
                     "ComprehensiveGenomeAnalysis")
    }
    n_t <- length(tools_sel)
    if (n_m == 0)
      div(class="alert alert-warning",
          style="padding:6px 10px;font-size:12px;margin-top:4px;",
        icon("exclamation-triangle"), " Sin muestras cargadas.")
    else if (n_t == 0)
      div(class="alert alert-warning",
          style="padding:6px 10px;font-size:12px;margin-top:4px;",
        icon("exclamation-triangle"), " Selecciona al menos una etapa.")
    else
      div(class="alert alert-info",
          style="padding:6px 10px;font-size:12px;margin-top:4px;",
        icon("info-circle"),
        paste0(" ", n_jobs <- n_m * n_t, " job(s) a enviar (",
               n_m, " muestra(s) × ", n_t, " etapa(s))"))
  })

  # ── Modal de confirmación ─────────────────────────────
  observeEvent(input$btn_launch_all, {
    if (length(rv$muestras) == 0) {
      showNotification("Carga archivos FASTQ primero.", type="warning", duration=4)
      return()
    }
    if (is.null(rv$bvbrc_token)) {
      showModal(creds_modal())
      showNotification("Inicia sesión BV-BRC primero.", type="warning", duration=5)
      return()
    }

    tools_sel <- input$batch_tools %||% character(0)
    use_cga   <- isTRUE(input$batch_cga)
    muestras  <- rv$muestras

    if (use_cga) {
      tools_sel <- c(tools_sel[!tools_sel %in% c("Assembly2","MetagenomicReadMapping")],
                     "ComprehensiveGenomeAnalysis")
    }
    if (!length(tools_sel)) {
      showNotification("Selecciona al menos una etapa.", type="warning"); return()
    }

    n_jobs    <- length(muestras) * length(tools_sel)
    mins_est  <- ceiling(length(muestras) * 1.5)   # ~90 seg por muestra para upload

    showModal(modalDialog(
      title     = tagList(icon("rocket"), " Confirmar lanzamiento por lotes"),
      size      = "m",
      easyClose = FALSE,

      div(class="alert alert-warning", style="padding:8px 12px;margin-bottom:10px;",
        icon("exclamation-triangle"),
        strong(" Resumen del lanzamiento por lotes:")
      ),
      tags$ul(style="margin:4px 0 10px;font-size:13px;",
        tags$li(paste0(strong(length(muestras)), " muestra(s): ",
                       paste(muestras, collapse=", "))),
        tags$li(paste0(strong(length(tools_sel)), " etapa(s): ",
                       paste(tools_sel, collapse=", "))),
        tags$li(strong(paste0("Total jobs a enviar: ", n_jobs)))
      ),
      div(class="alert alert-info", style="padding:8px 12px;font-size:12px;",
        icon("info-circle"),
        " Los archivos FASTQ se subirán a BV-BRC ",
        strong("una sola vez por muestra"),
        " y luego se enviarán todos los jobs. ",
        "Los resultados se importarán automáticamente al completarse."
      ),
      div(class="alert alert-warning", style="padding:8px 12px;font-size:12px;",
        icon("clock"),
        paste0(" Tiempo estimado de upload: ~", mins_est, " min. "),
        "Puedes navegar en la app mientras se envían los jobs."
      ),

      footer = tagList(
        actionButton("btn_launch_confirm",
          tagList(icon("rocket"), " Confirmar — Lanzar ahora"),
          class = "btn-danger btn-lg"),
        modalButton("Cancelar")
      )
    ))
  }, ignoreInit = TRUE)

  # ── Ejecución batch ────────────────────────────────────
  observeEvent(input$btn_launch_confirm, {
    removeModal()

    tools_sel <- isolate(input$batch_tools) %||% character(0)
    use_cga   <- isolate(isTRUE(input$batch_cga))
    muestras  <- isolate(rv$muestras)
    organism  <- isolate(input$proy_organismo) %||% "Unknown"
    tok       <- isolate(rv$bvbrc_token)
    user      <- isolate(rv$bvbrc_user)

    if (use_cga) {
      tools_sel <- c(tools_sel[!tools_sel %in% c("Assembly2","MetagenomicReadMapping")],
                     "ComprehensiveGenomeAnalysis")
    }

    if (!length(tools_sel) || !length(muestras) || is.null(tok) || is.null(user)) {
      showNotification("Faltan muestras, etapas o sesión BV-BRC.", type="error"); return()
    }

    n_total <- length(muestras) * length(tools_sel)
    n_ok    <- 0L; n_err <- 0L

    alog("INFO", "[BATCH] Iniciando lanzamiento por lotes",
         paste0("Muestras: ", paste(muestras, collapse=", "),
                "\nEtapas: ",  paste(tools_sel, collapse=", "),
                "\nTotal jobs planificados: ", n_total))

    # Crear carpeta de salida una sola vez antes del loop
    tryCatch(bvbrc_mkdir(tok, ws_out(user)), error = function(e) NULL)
    # Sincronizar registro con el workspace (evita re-analizar lo ya hecho)
    tryCatch(bvbrc_sync_registry(tok, user), error = function(e) NULL)

    withProgress(message = "Lanzando análisis por lotes...", value = 0, {

      for (i in seq_along(muestras)) {
        muestra <- muestras[i]

        incProgress(0,
          detail = paste0("Muestra ", i, "/", length(muestras),
                          " — subiendo archivos: ", muestra))

        # ── Upload FASTQ una vez por muestra ─────────────
        sfiles <- get_sample_files(muestra)
        r1 <- sfiles$r1; r2 <- sfiles$r2

        ws_r1 <- paste0(ws_home(user), "/.fastq_uploads/", muestra, "_R1.fastq.gz")
        ws_r2 <- NULL  # NULL = single-end por defecto

        if (!is.null(r1) && file.exists(r1)) {
          alog("API", paste0("[BATCH] Upload R1: ", basename(r1), " — ", muestra))
          up1 <- tryCatch(
            bvbrc_upload_file(tok, user, r1),
            error = function(e) list(success=FALSE, error=e$message))
          if (isTRUE(up1$success)) {
            ws_r1 <- up1$ws_path
            alog("OK", paste0("[BATCH] R1 subido: ", muestra), up1$log_det)
          } else {
            alog("WARN", paste0("[BATCH] R1 upload fallido: ", muestra),
                 up1$error %||% "sin detalle")
          }
        }

        if (!is.null(r2) && file.exists(r2) && !identical(r2, r1)) {
          alog("API", paste0("[BATCH] Upload R2: ", basename(r2), " — ", muestra))
          up2 <- tryCatch(
            bvbrc_upload_file(tok, user, r2),
            error = function(e) list(success=FALSE, error=e$message))
          if (isTRUE(up2$success)) {
            ws_r2 <- up2$ws_path
            alog("OK", paste0("[BATCH] R2 subido: ", muestra), up2$log_det)
          } else {
            alog("WARN", paste0("[BATCH] R2 upload fallido — modo single-end: ", muestra),
                 up2$error %||% "sin detalle")
          }
        } else if (is.null(r2)) {
          alog("DEBUG", paste0("[BATCH] R2 no disponible — modo single-end: ", muestra))
        }

        # ── Enviar un job por etapa ───────────────────────
        for (j in seq_along(tools_sel)) {
          tool_id <- tools_sel[j]
          step_id <- BATCH_STEP_MAP[tool_id] %||% "calidad"
          job_num <- (i - 1L) * length(tools_sel) + j

          # ── Skip si ya existe un job completado ──────────────
          if (!is.null(bvbrc_get_done_job(rv$jobs, tool_id, muestra))) {
            alog("INFO", paste0("[BATCH] Skip (ya completado): ", tool_id, " — ", muestra))
            n_ok <- n_ok + 1L
            next
          }

          incProgress(
            1 / n_total,
            detail = paste0("Enviando ", job_num, "/", n_total,
                            ": ", tool_id, " — ", muestra)
          )

          params <- tryCatch(switch(tool_id,
            FastqUtils                  = bvbrc_params_fastqutil(user, ws_r1, ws_r2, muestra),
            TaxonomicClassification     = bvbrc_params_taxonomy(user, ws_r1, ws_r2, muestra),
            Assembly2                   = bvbrc_params_assembly(user, ws_r1, ws_r2, muestra),
            MetagenomicReadMapping      = bvbrc_params_metareads(user, ws_r1, ws_r2, muestra),
            ComprehensiveGenomeAnalysis = bvbrc_params_cga(user, ws_r1, ws_r2,
                                                            organism, muestra),
            list()
          ), error = function(e) list())

          res <- tryCatch(
            bvbrc_submit_job(tok, tool_id, params, ws_home(user)),
            error = function(e) list(success=FALSE, error=e$message))

          if (isTRUE(res$success)) {
            add_job(res$task_id, "bvbrc", tool_id, step_id,
                    paste0("[Lote] ", tool_id, " — ", muestra),
                    output_path = params$output_path,
                    output_file = params$output_file,
                    muestra     = muestra)
            n_ok <- n_ok + 1L
            alog("OK", paste0("[BATCH] Job enviado: ", tool_id, " — ", muestra),
                 paste0("task_id: ", res$task_id,
                        "\n", res$log_det %||% ""))
          } else {
            n_err <- n_err + 1L
            alog("ERROR", paste0("[BATCH] Submit fallido: ", tool_id, " — ", muestra),
                 res$error %||% "sin detalle")
          }
        }
      }   # fin loop muestras
    })   # fin withProgress

    msg_fin <- paste0(n_ok, "/", n_total, " jobs enviados a BV-BRC")
    if (n_err > 0) msg_fin <- paste0(msg_fin, " (", n_err, " errores — ver Log)")
    alog("OK", "[BATCH] Lanzamiento por lotes completado",
         paste0("Enviados: ", n_ok, " | Errores: ", n_err, " | Total: ", n_total))
    showNotification(
      tagList(icon("rocket"), " ", msg_fin),
      type     = if (n_err == 0) "message" else "warning",
      duration = 12)

  }, ignoreInit = TRUE)

  # ============================================================
  # PIPELINE AUTOMÁTICO (pasos 1→9, orden dependiente de resultados)
  # ============================================================
  #
  # Olas de dependencia:
  #   A. Solo requieren FASTQ (en paralelo, de inmediato):
  #        FastqUtils(1), TaxonomicClassification(2), Assembly2(3),
  #        MetagenomicReadMapping(5)
  #   B. Requiere ensamblado terminado (por muestra):
  #        Annotation(4)  ← usa los contigs del ensamblado
  #   C. Requiere anotación terminada (por muestra):
  #        CoreGenomeMLST(7)  ← usa el genome_id de la anotación
  #   D. Requiere TODAS las anotaciones (global, una vez, ≥2 genomas):
  #        PhylogeneticTree(8) + alineamiento(6, interno) + Gen.Comp.(9)

  # Etapa → step_id de la UI (para marcar rv$completados)
  PIPE_STEP_OF <- c(
    FastqUtils              = "calidad",
    TaxonomicClassification = "taxonomia",
    Assembly2               = "ensamblado",
    Annotation              = "anotacion",
    MetagenomicReadMapping  = "resistoma",
    CoreGenomeMLST          = "mlst",
    PhylogeneticTree        = "filogenia"
  )

  pipe_is_submitted <- function(m, tool) {
    tool %in% (isolate(rv$auto_pipe)$submitted[[m]] %||% character(0))
  }
  pipe_mark_submitted <- function(m, tool) {
    ap <- rv$auto_pipe
    ap$submitted[[m]] <- unique(c(ap$submitted[[m]] %||% character(0), tool))
    rv$auto_pipe <- ap
  }

  # Buscar el archivo de contigs dentro de la salida de un ensamblado
  pipe_find_contigs <- function(tok, user, m) {
    out_dir <- paste0(ws_out(user), "/assembly_", m)
    ls <- tryCatch(bvbrc_ls(tok, out_dir), error = function(e) list(success = FALSE))
    if (!isTRUE(ls$success) || !length(ls$files)) return(NULL)
    is_contig <- function(f) {
      identical(tolower(f$type %||% ""), "contigs") ||
      grepl("contigs.*\\.f(ast)?a$|_contigs\\.f(ast)?a$", f$name, ignore.case = TRUE)
    }
    cand <- Filter(is_contig, ls$files)
    if (!length(cand))
      cand <- Filter(function(f) grepl("\\.f(ast)?a$", f$name, ignore.case = TRUE), ls$files)
    # Buscar en subcarpetas si no aparece en la raíz
    if (!length(cand)) {
      for (f in ls$files) {
        if ((f$type %||% "") %in% c("folder","directory","job_result")) {
          ls2 <- tryCatch(bvbrc_ls(tok, f$path), error = function(e) list(success = FALSE))
          if (isTRUE(ls2$success)) {
            sub <- Filter(is_contig, ls2$files)
            if (length(sub)) { cand <- sub; break }
          }
        }
      }
    }
    if (length(cand)) cand[[1]]$path else NULL
  }

  # Extraer el genome_id del resultado de una anotación
  pipe_find_genome_id <- function(tok, user, m) {
    out_dir <- paste0(ws_out(user), "/annotation_", m)
    ls <- tryCatch(bvbrc_ls(tok, out_dir), error = function(e) list(success = FALSE))
    if (!isTRUE(ls$success) || !length(ls$files)) return(NULL)
    # 1) nombre de archivo/carpeta con patrón genome_id  <taxid>.<n>
    for (f in ls$files) {
      mm <- regmatches(f$name, regexpr("^[0-9]+\\.[0-9]+", f$name))
      if (length(mm) && nchar(mm)) return(mm)
    }
    # 2) archivo pequeño "genome_id" leído directamente
    gidf <- Filter(function(f) grepl("genome_id", f$name, ignore.case = TRUE), ls$files)
    if (length(gidf)) {
      g <- tryCatch(bvbrc_get_file(tok, gidf[[1]]$path), error = function(e) list(success = FALSE))
      if (isTRUE(g$success)) {
        val <- regmatches(g$content, regexpr("[0-9]+\\.[0-9]+", g$content))
        if (length(val) && nchar(val)) return(val)
      }
    }
    # 3) buscar patrón dentro de un .genome / .gto pequeño
    genf <- Filter(function(f) grepl("\\.(genome|gto)$", f$name, ignore.case = TRUE) &&
                     (f$size %||% 0) < 5e6, ls$files)
    if (length(genf)) {
      g <- tryCatch(bvbrc_get_file(tok, genf[[1]]$path), error = function(e) list(success = FALSE))
      if (isTRUE(g$success)) {
        val <- regmatches(g$content, regexpr("\"id\"\\s*:\\s*\"([0-9]+\\.[0-9]+)\"", g$content))
        if (length(val)) {
          num <- regmatches(val, regexpr("[0-9]+\\.[0-9]+", val))
          if (length(num)) return(num)
        }
      }
    }
    NULL
  }

  # Enviar un job del pipeline y registrarlo
  pipe_submit <- function(tok, user, tool_id, params, m, label) {
    res <- tryCatch(
      bvbrc_submit_job(tok, tool_id, params, ws_home(user)),
      error = function(e) list(success = FALSE, error = conditionMessage(e)))
    if (isTRUE(res$success)) {
      add_job(res$task_id, "bvbrc", tool_id,
              PIPE_STEP_OF[tool_id] %||% "calidad", label,
              output_path = params$output_path,
              output_file = params$output_file,
              muestra     = m)
      pipe_mark_submitted(m %||% "_global", tool_id)
      alog("OK", paste0("[PIPE] Enviado: ", tool_id, " — ", m %||% "(global)"),
           paste0("task_id: ", res$task_id))
      TRUE
    } else {
      alog("ERROR", paste0("[PIPE] Falló envío: ", tool_id, " — ", m %||% "(global)"),
           res$error %||% "sin detalle")
      FALSE
    }
  }

  # ── Galaxy dentro del pipeline (taxonomía con Kraken2) ─────
  # ¿Ya hay un job de Galaxy completado para esta herramienta/muestra?
  pipe_galaxy_done <- function(tool, m) {
    any(vapply(isolate(rv$jobs), function(j)
      isTRUE(j$platform == "galaxy") && isTRUE(j$tool == tool) &&
      isTRUE((j$muestra %||% "") == m) &&
      isTRUE(j$status %in% c("completed","ok","done","FINISHED")),
      logical(1)))
  }

  # Sube R1 a Galaxy (una vez por muestra) y lanza Kraken2.
  pipe_submit_kraken2 <- function(m) {
    key <- isolate(rv$galaxy_key)
    if (is.null(key)) {
      alog("WARN", "[PIPE] Sin API key de Galaxy — taxonomía omitida",
           "Configura la clave de Galaxy para clasificar taxonómicamente.")
      return(invisible(FALSE))
    }
    hist <- isolate(rv$galaxy_hist)
    if (is.null(hist)) {
      h <- tryCatch(galaxy_get_or_create_history(key), error = function(e) list(success = FALSE))
      if (!isTRUE(h$success)) {
        alog("ERROR", "[PIPE] No se pudo obtener historial de Galaxy", h$error %||% "")
        return(invisible(FALSE))
      }
      hist <- h$history_id; rv$galaxy_hist <- hist
    }

    ap <- isolate(rv$auto_pipe)
    ds <- ap$galaxy_ds[[m]] %||% NULL
    if (is.null(ds)) {
      sf <- get_sample_files(m); r1 <- sf$r1
      if (is.null(r1) || !file.exists(r1)) {
        alog("WARN", paste0("[PIPE] Sin R1 local para Galaxy: ", m)); return(invisible(FALSE))
      }
      alog("API", paste0("[PIPE] Subiendo a Galaxy para taxonomía: ", m), basename(r1))
      up <- tryCatch(galaxy_upload_file(key, hist, r1),
                     error = function(e) list(success = FALSE, error = conditionMessage(e)))
      if (!isTRUE(up$success)) {
        alog("ERROR", paste0("[PIPE] Upload Galaxy fallido: ", m), up$error %||% "")
        return(invisible(FALSE))
      }
      w <- tryCatch(galaxy_wait_dataset(key, up$dataset_id, hist, max_wait = 900, poll_secs = 10),
                    error = function(e) list(success = FALSE))
      if (!isTRUE(w$success)) {
        alog("WARN", paste0("[PIPE] Dataset Galaxy no listo aún: ", m),
             "Se reintentará en el siguiente ciclo.")
        return(invisible(FALSE))
      }
      ds <- up$dataset_id
      ap$galaxy_ds[[m]] <- ds; rv$auto_pipe <- ap
    }

    res <- tryCatch(
      galaxy_run_tool(key, hist, GALAXY_TOOL_KRAKEN2, galaxy_inputs_kraken2(ds)),
      error = function(e) list(success = FALSE, error = conditionMessage(e)))
    if (isTRUE(res$success) && length(res$job_ids) > 0) {
      add_job(res$job_ids[1], "galaxy", "kraken2", "taxonomia",
              paste0("[Auto] Kraken2 (Galaxy) — ", m),
              output_ids = res$output_ids, muestra = m)
      pipe_mark_submitted(m, "kraken2")
      alog("OK", paste0("[PIPE] Kraken2 enviado a Galaxy — ", m),
           paste0("job_id: ", res$job_ids[1]))
      invisible(TRUE)
    } else {
      alog("ERROR", paste0("[PIPE] Kraken2 falló al enviarse — ", m),
           res$error %||% "sin detalle")
      invisible(FALSE)
    }
  }

  # MLST vía Galaxy: BV-BRC no tiene app CoreGenomeMLST (verificado con
  # enumerate_apps). Se descargan los contigs del ensamblado de BV-BRC,
  # se suben a Galaxy y se corre la herramienta mlst.
  pipe_submit_mlst <- function(m, contigs_ws) {
    key <- isolate(rv$galaxy_key)
    if (is.null(key)) {
      alog("WARN", "[PIPE] Sin API key de Galaxy — MLST omitido"); return(invisible(FALSE))
    }
    hist <- isolate(rv$galaxy_hist)
    if (is.null(hist)) {
      h <- tryCatch(galaxy_get_or_create_history(key), error = function(e) list(success = FALSE))
      if (!isTRUE(h$success)) return(invisible(FALSE))
      hist <- h$history_id; rv$galaxy_hist <- hist
    }
    tok <- isolate(rv$bvbrc_token)
    g <- tryCatch(bvbrc_get_file(tok, contigs_ws), error = function(e) list(success = FALSE))
    if (!isTRUE(g$success) || !nchar(trimws(g$content %||% ""))) {
      alog("WARN", paste0("[PIPE] No se pudieron descargar contigs para MLST: ", m))
      return(invisible(FALSE))
    }
    tmp <- file.path(tempdir(), paste0(m, "_contigs.fasta"))
    writeLines(g$content, tmp)
    up <- tryCatch(galaxy_upload_file(key, hist, tmp, file_type = "fasta"),
                   error = function(e) list(success = FALSE, error = conditionMessage(e)))
    if (!isTRUE(up$success)) {
      alog("ERROR", paste0("[PIPE] Upload de contigs a Galaxy falló: ", m), up$error %||% "")
      return(invisible(FALSE))
    }
    w <- tryCatch(galaxy_wait_dataset(key, up$dataset_id, hist, max_wait = 600, poll_secs = 10),
                  error = function(e) list(success = FALSE))
    if (!isTRUE(w$success)) return(invisible(FALSE))
    res <- tryCatch(galaxy_run_tool(key, hist, GALAXY_TOOL_MLST,
                                    galaxy_inputs_mlst(up$dataset_id)),
                    error = function(e) list(success = FALSE, error = conditionMessage(e)))
    if (isTRUE(res$success) && length(res$job_ids) > 0) {
      add_job(res$job_ids[1], "galaxy", "mlst", "mlst",
              paste0("[Auto] MLST (Galaxy) — ", m),
              output_ids = res$output_ids, muestra = m)
      pipe_mark_submitted(m, "mlst")
      alog("OK", paste0("[PIPE] MLST enviado a Galaxy — ", m),
           paste0("job_id: ", res$job_ids[1]))
      invisible(TRUE)
    } else {
      alog("ERROR", paste0("[PIPE] MLST falló al enviarse — ", m), res$error %||% "")
      invisible(FALSE)
    }
  }

  # ── Motor: avanza el pipeline un paso según el estado actual ──
  pipe_tick <- function() {
    ap <- isolate(rv$auto_pipe)
    if (!isTRUE(ap$active)) return(invisible())
    tok  <- isolate(rv$bvbrc_token)
    user <- isolate(rv$bvbrc_user)
    if (is.null(tok) || is.null(user)) return(invisible())
    jobs <- isolate(rv$jobs)
    organism <- ap$organism %||% "Unknown"

    done <- function(tool, m) !is.null(bvbrc_get_done_job(jobs, tool, m))

    # ── Olas A/B/C por muestra ──
    for (m in ap$samples) {
      wsm   <- ap$ws[[m]] %||% list(r1 = NULL, r2 = NULL)
      ws_r1 <- wsm$r1; ws_r2 <- wsm$r2
      if (is.null(ws_r1)) next   # sin upload no se puede enviar

      # Ola A — solo FASTQ (BV-BRC)
      # NOTA: la taxonomía NO va por BV-BRC: su servicio
      # TaxonomicClassification falla siempre ("wrapper command failed 256").
      # Se resuelve con Kraken2 en Galaxy (ver pipe_submit_kraken2).
      waveA <- list(
        FastqUtils              = bvbrc_params_fastqutil(user, ws_r1, ws_r2, m),
        Assembly2               = bvbrc_params_assembly (user, ws_r1, ws_r2, m),
        MetagenomicReadMapping  = bvbrc_params_metareads(user, ws_r1, ws_r2, m)
      )
      for (tool in names(waveA)) {
        if (!pipe_is_submitted(m, tool) && !done(tool, m))
          pipe_submit(tok, user, tool, waveA[[tool]], m, paste0("[Auto] ", tool, " — ", m))
      }

      # Ola A (Galaxy) — Taxonomía con Kraken2
      if (!pipe_is_submitted(m, "kraken2") && !pipe_galaxy_done("kraken2", m))
        pipe_submit_kraken2(m)

      # Ola B — Anotación tras ensamblado
      if (done("Assembly2", m) && !pipe_is_submitted(m, "Annotation") &&
          !done("Annotation", m)) {
        contigs <- ap$contigs[[m]] %||% pipe_find_contigs(tok, user, m)
        if (!is.null(contigs)) {
          ap$contigs[[m]] <- contigs; rv$auto_pipe <- ap
          pipe_submit(tok, user, "Annotation",
                      bvbrc_params_annotation(user, contigs, organism, m),
                      m, paste0("[Auto] Annotation — ", m))
        } else {
          alog("INFO", paste0("[PIPE] Esperando contigs del ensamblado: ", m))
        }
      }

      # Ola C — MLST (Galaxy) tras el ensamblado: usa los contigs
      if (done("Assembly2", m) && !pipe_is_submitted(m, "mlst") &&
          !pipe_galaxy_done("mlst", m)) {
        ctg <- ap$contigs[[m]] %||% pipe_find_contigs(tok, user, m)
        if (!is.null(ctg)) {
          ap$contigs[[m]] <- ctg; rv$auto_pipe <- ap
          pipe_submit_mlst(m, ctg)
        } else {
          alog("INFO", paste0("[PIPE] Esperando contigs para MLST: ", m))
        }
      }
    }

    # ── Ola D — global: filogenia + comparativa (una sola vez) ──
    ap <- isolate(rv$auto_pipe)   # releer por si genome_ids cambió
    all_annot_done <- length(ap$samples) > 0 &&
      all(vapply(ap$samples, function(m) done("Annotation", m), logical(1)))
    gids <- unique(unlist(ap$genome_ids, use.names = FALSE))

    if (all_annot_done && length(gids) >= 2 &&
        !("PhylogeneticTree" %in% ap$global)) {
      okp <- pipe_submit(tok, user, "PhylogeneticTree",
                  bvbrc_params_phylogeny(user, as.list(gids), "filogenia_auto"),
                  NULL, "[Auto] Filogenia (todas las muestras)")
      if (okp) {
        ap$global <- unique(c(ap$global, "PhylogeneticTree")); rv$auto_pipe <- ap
        alog("INFO", "[PIPE] Alineamiento (paso 6) resuelto dentro de la filogenia",
             "PhylogeneticTree calcula el alineamiento múltiple de genes core.")
      }
    }
    if (all_annot_done && length(gids) >= 2 &&
        !("SeqComparison" %in% ap$global)) {
      ok9 <- pipe_submit(tok, user, "SeqComparison",
                  bvbrc_params_proteome(user, gids[[1]], gids[[min(2,length(gids))]],
                                        "comparativa_auto"),
                  NULL, "[Auto] Gen. Comparativa")
      if (ok9) { ap$global <- unique(c(ap$global, "SeqComparison")); rv$auto_pipe <- ap }
    }

    # ── ¿Pipeline terminado? ──
    ap <- isolate(rv$auto_pipe)
    per_sample_done <- all(vapply(ap$samples, function(m)
      done("FastqUtils", m) && pipe_galaxy_done("kraken2", m) &&
      done("Assembly2", m) && done("Annotation", m) &&
      done("MetagenomicReadMapping", m) && pipe_galaxy_done("mlst", m),
      logical(1)))
    global_done <- all(c("PhylogeneticTree","SeqComparison") %in% ap$global)
    if (length(ap$samples) > 0 && per_sample_done && global_done) {
      ap$active <- FALSE; rv$auto_pipe <- ap
      alog("OK", "[PIPE] Pipeline automático COMPLETADO",
           paste0("Todas las etapas 1→9 finalizadas para ",
                  paste(ap$samples, collapse = ", ")))
      showNotification(
        tagList(icon("check-circle"), " Pipeline automático completado (1→9)."),
        type = "message", duration = 12)
    }
    invisible()
  }

  # ── El motor corre en cada ciclo del poll_timer ──
  observe({
    poll_timer()
    if (isTRUE(isolate(rv$auto_pipe)$active))
      tryCatch(pipe_tick(), error = function(e)
        alog("ERROR", "[PIPE] Error en el motor del pipeline", conditionMessage(e)))
  })

  # ── Arrancar el pipeline (reutilizable: botón y asistente) ──
  pipe_arrancar <- function(organismo = NULL) {
    if (length(rv$muestras) == 0) {
      showNotification("Carga archivos FASTQ primero.", type = "warning", duration = 4); return(invisible(FALSE))
    }
    tok  <- rv$bvbrc_token; user <- rv$bvbrc_user
    if (is.null(tok)) {
      showModal(creds_modal())
      showNotification("Inicia sesión BV-BRC primero.", type = "warning", duration = 5); return(invisible(FALSE))
    }
    muestras <- rv$muestras
    organism <- organismo %||% input$proy_organismo %||% "Unknown"

    showModal(modalDialog(
      title = tagList(icon("robot"), " Iniciando pipeline automático..."),
      tags$p(icon("spinner", class = "fa-spin"),
             " Preparando archivos y lanzando las primeras etapas..."),
      tags$small(style = "color:#888;",
        "Se subirán los FASTQ que aún no estén en BV-BRC. Las etapas ",
        "posteriores se lanzarán solas conforme terminen sus dependencias."),
      footer = NULL, easyClose = FALSE
    ))

    tryCatch(bvbrc_mkdir(tok, ws_out(user)), error = function(e) NULL)
    # Sincronizar registro con el workspace (evita re-analizar lo ya hecho)
    tryCatch(bvbrc_sync_registry(tok, user), error = function(e) NULL)

    ws_map <- list()
    for (m in muestras) {
      # Si ya corrió FastqUtils, el upload ya existe → ruta determinista
      det_r1 <- paste0(ws_home(user), "/.fastq_uploads/", m, "_R1.fastq.gz")
      if (!is.null(bvbrc_get_done_job(rv$jobs, "FastqUtils", m))) {
        ws_map[[m]] <- list(r1 = det_r1, r2 = NULL)
        alog("INFO", paste0("[PIPE] Reutilizando upload existente: ", m), det_r1)
        next
      }
      sfiles <- get_sample_files(m)
      r1 <- sfiles$r1; r2 <- sfiles$r2
      ws_r1 <- det_r1; ws_r2 <- NULL
      if (!is.null(r1) && file.exists(r1)) {
        up1 <- tryCatch(bvbrc_upload_file(tok, user, r1),
                        error = function(e) list(success = FALSE, error = e$message))
        if (isTRUE(up1$success)) ws_r1 <- up1$ws_path
        else alog("WARN", paste0("[PIPE] Upload R1 falló: ", m), up1$error %||% "")
        if (!is.null(r2) && file.exists(r2) && !identical(r2, r1)) {
          up2 <- tryCatch(bvbrc_upload_file(tok, user, r2),
                          error = function(e) list(success = FALSE, error = e$message))
          if (isTRUE(up2$success)) ws_r2 <- up2$ws_path
        }
      }
      ws_map[[m]] <- list(r1 = ws_r1, r2 = ws_r2)
    }

    rv$auto_pipe <- list(
      active = TRUE, samples = muestras, organism = organism,
      ws = ws_map, submitted = list(), contigs = list(),
      genome_ids = list(), global = character(0), started_at = Sys.time()
    )
    tryCatch(state_save(rv), error = function(e) NULL)
    removeModal()
    alog("OK", "[PIPE] Pipeline automático INICIADO",
         paste0("Muestras: ", paste(muestras, collapse = ", "),
                "\nOrganismo: ", organism))
    showNotification(
      tagList(icon("robot"), " Pipeline automático iniciado. Las etapas se ",
              "lanzarán solas según terminen sus dependencias."),
      type = "message", duration = 10)

    # Primer avance inmediato (no esperar 45 s)
    tryCatch(pipe_tick(), error = function(e)
      alog("ERROR", "[PIPE] Error en primer avance", conditionMessage(e)))
    # Guardar la sesión con el pipeline ya arrancado
    tryCatch(guardar_sesion(verbose = FALSE), error = function(e) NULL)
    invisible(TRUE)
  }

  observeEvent(input$btn_auto_pipeline, pipe_arrancar(), ignoreInit = TRUE)

  # ── Botón: detener pipeline automático ──
  observeEvent(input$btn_auto_pipeline_stop, {
    ap <- rv$auto_pipe
    if (!isTRUE(ap$active)) {
      showNotification("El pipeline automático no está activo.", type = "message", duration = 4)
      return()
    }
    ap$active <- FALSE; rv$auto_pipe <- ap
    tryCatch(state_save(rv), error = function(e) NULL)
    alog("WARN", "[PIPE] Pipeline automático DETENIDO por el usuario",
         "Los jobs ya enviados siguen corriendo en BV-BRC.")
    showNotification(
      "Pipeline detenido. Los jobs ya enviados continúan en BV-BRC.",
      type = "warning", duration = 6)
  }, ignoreInit = TRUE)

  # ── Estado del pipeline (UI) ──
  output$auto_pipe_status_ui <- renderUI({
    invalidateLater(8000, session)
    ap <- rv$auto_pipe
    if (!isTRUE(ap$active) && is.null(ap$started_at)) {
      return(div(class = "alert alert-secondary",
                 style = "padding:7px 10px;font-size:12px;margin-bottom:6px;",
        icon("info-circle"), " Pipeline automático inactivo."))
    }
    jobs <- rv$jobs
    done <- function(tool, m) !is.null(bvbrc_get_done_job(jobs, tool, m))
    tools_per <- c("FastqUtils","Assembly2","Annotation","MetagenomicReadMapping")
    # +2 por las etapas en Galaxy (Kraken2 taxonomía y MLST)
    n_cells <- length(ap$samples) * (length(tools_per) + 2L)
    n_ok <- if (length(ap$samples))
      sum(vapply(ap$samples, function(m)
        sum(vapply(tools_per, function(t) done(t, m), logical(1))) +
        as.integer(pipe_galaxy_done("kraken2", m)) +
        as.integer(pipe_galaxy_done("mlst", m)), integer(1)))
      else 0L
    pct <- if (n_cells > 0) round(100 * n_ok / n_cells) else 0
    estado <- if (isTRUE(ap$active)) "EN CURSO" else "detenido/completado"
    cls    <- if (isTRUE(ap$active)) "alert-info" else "alert-success"
    div(class = paste("alert", cls),
        style = "padding:8px 12px;font-size:12px;margin-bottom:6px;",
      icon(if (isTRUE(ap$active)) "spinner" else "check-circle",
           class = if (isTRUE(ap$active)) "fa-spin" else NULL),
      strong(paste0(" Pipeline ", estado, " — ")),
      paste0(n_ok, "/", n_cells, " etapas por muestra completadas (", pct, "%). "),
      if (length(ap$global))
        paste0("Global: ", paste(ap$global, collapse = ", "), ".")
    )
  })

  # ============================================================
  # BACKUP DE RESULTADOS (exportar / importar JSON)
  # ============================================================

  # ── Exportar ─────────────────────────────────────────────
  output$dl_results_backup <- downloadHandler(
    filename = function() {
      paste0("anbio_backup_", format(Sys.Date(), "%Y%m%d_%H%M%S"), ".json")
    },
    content = function(file) {
      backup <- list(
        version   = "1.0",
        app       = "ANbio LESP",
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        proyecto  = list(
          nombre      = isolate(input$proy_nombre)      %||% "",
          organismo   = isolate(input$proy_organismo)   %||% "",
          analista    = isolate(input$proy_analista)    %||% "",
          institucion = isolate(input$proy_institucion) %||% "",
          fecha       = format(isolate(input$proy_fecha) %||% Sys.Date(), "%Y-%m-%d"),
          tipo        = isolate(input$proy_tipo)        %||% ""
        ),
        muestras     = as.list(isolate(rv$muestras)),
        completados  = as.list(isolate(rv$completados)),
        calidad      = if (nrow(isolate(rv$calidad))      > 0) isolate(rv$calidad)      else NULL,
        taxonomia    = if (nrow(isolate(rv$taxonomia))    > 0) isolate(rv$taxonomia)    else NULL,
        ensamblado   = if (nrow(isolate(rv$ensamblado))   > 0) isolate(rv$ensamblado)   else NULL,
        mlst         = if (nrow(isolate(rv$mlst))         > 0) isolate(rv$mlst)         else NULL,
        alineamiento = if (nrow(isolate(rv$alineamiento)) > 0) isolate(rv$alineamiento) else NULL,
        ani          = if (nrow(isolate(rv$ani))          > 0) isolate(rv$ani)          else NULL,
        resfinder_df = isolate(rv$resfinder_df),
        jobs_summary = lapply(isolate(rv$jobs), function(j) {
          campos <- c("job_id","platform","tool","step","label","status","muestra")
          out    <- j[intersect(campos, names(j))]
          out$time <- format(j$time %||% Sys.time(), "%Y-%m-%dT%H:%M:%S")
          out
        })
      )
      jsonlite::write_json(backup, file, auto_unbox = TRUE, pretty = TRUE, na = "null")
      alog("OK", "Backup JSON exportado",
           paste0("Tablas incluidas: calidad/taxonomia/ensamblado/mlst/alineamiento/ani/resfinder",
                  "\nMuestras: ", length(isolate(rv$muestras)),
                  "\nJobs: ", length(isolate(rv$jobs))))
    }
  )

  # ── Importar ──────────────────────────────────────────────
  observeEvent(input$backup_file, {
    req(input$backup_file)
    path <- input$backup_file$datapath

    tryCatch({
      backup <- jsonlite::read_json(path, simplifyVector = TRUE)

      if (is.null(backup$version)) {
        showNotification("Archivo no es un backup ANbio válido.", type="error", duration=6)
        return()
      }

      # Convierte a data.frame si vino como lista/NULL
      df_restore <- function(x) {
        if (is.null(x)) return(NULL)
        if (is.data.frame(x) && nrow(x) > 0) return(x)
        if (is.data.frame(x) && nrow(x) == 0) return(NULL)
        tryCatch(as.data.frame(x, stringsAsFactors = FALSE), error = function(e) NULL)
      }

      cal  <- df_restore(backup$calidad)
      tax  <- df_restore(backup$taxonomia)
      ens  <- df_restore(backup$ensamblado)
      mls  <- df_restore(backup$mlst)
      ali  <- df_restore(backup$alineamiento)
      ani  <- df_restore(backup$ani)
      res  <- df_restore(backup$resfinder_df)
      mus  <- if (!is.null(backup$muestras))
                unique(trimws(unlist(backup$muestras))) else character(0)
      mus  <- mus[nchar(mus) > 0]

      comp_raw <- backup$completados %||% list()
      comp     <- setNames(rep(FALSE, length(PASOS_IDS)), PASOS_IDS)
      for (pid in PASOS_IDS) {
        v <- comp_raw[[pid]]
        if (!is.null(v)) comp[pid] <- isTRUE(as.logical(v))
      }

      # Restaurar reactiveValues
      if (!is.null(cal))   rv$calidad      <- cal
      if (!is.null(tax))   rv$taxonomia    <- tax
      if (!is.null(ens))   rv$ensamblado   <- ens
      if (!is.null(mls))   rv$mlst         <- mls
      if (!is.null(ali))   rv$alineamiento <- ali
      if (!is.null(ani))   rv$ani          <- ani
      if (!is.null(res))   rv$resfinder_df <- res
      if (length(mus) > 0) rv$muestras     <- mus
      for (pid in PASOS_IDS) rv$completados[pid] <- comp[pid]

      n_tablas <- sum(!vapply(list(cal, tax, ens, mls, ali, ani, res),
                               is.null, logical(1)))

      alog("OK", "Backup JSON importado",
           paste0("Archivo: ", input$backup_file$name,
                  "\nProyecto: ", backup$proyecto$nombre %||% "?",
                  " | Timestamp backup: ", backup$timestamp %||% "?",
                  "\nTablas restauradas: ", n_tablas,
                  " | Muestras: ", length(mus)))

      rv$backup_import_status <- paste0(
        "✓ Importado el ", format(Sys.time(), "%d/%m/%Y %H:%M"),
        " — Proyecto: ", backup$proyecto$nombre %||% "(sin nombre)",
        " — ", n_tablas, " tabla(s), ", length(mus), " muestra(s)")

      showNotification(
        tagList(icon("check-circle"),
          paste0(" Backup importado: ", n_tablas, " tabla(s), ",
                 length(mus), " muestra(s).")),
        type = "message", duration = 8)

      # Guardar estado con los datos restaurados
      tryCatch(state_save(rv), error = function(e) NULL)

    }, error = function(e) {
      alog("ERROR", "Error importando backup JSON", e$message)
      showNotification(paste("Error al importar:", e$message), type="error", duration=8)
    })
  }, ignoreInit = TRUE)

  output$backup_import_status <- renderUI({
    st <- rv$backup_import_status
    if (is.null(st)) return(NULL)
    div(class="alert alert-success",
        style="padding:6px 10px;font-size:11px;margin-top:4px;",
      icon("check-circle"), " ", st)
  })

  # ── Recuperación de resultados: escaneo del workspace BV-BRC ───────────────
  # Escanea home/ en 2 niveles de profundidad + ruta manual.
  # NO depende de task IDs.

  TOOL_TO_STEP <- c(
    FastqUtils = "calidad", TaxonomicClassification = "taxonomia",
    Assembly2 = "ensamblado", Annotation = "anotacion",
    ComprehensiveGenomeAnalysis = "ensamblado",
    MetagenomicReadMapping = "resistoma", CoreGenomeMLST = "mlst",
    PhylogeneticTree = "filogenia",
    # "genomica" es el step_id real del paso 9 en PASOS_IDS
    # (antes decía "comparacion", que no existe → el paso nunca se marcaba).
    SeqComparison = "genomica", Homology = "alineamiento"
  )

  # Prefijos de carpeta → nombre de herramienta (según bvbrc_params_*)
  WS_FOLDER_MAP <- list(
    FastqUtils              = "^fastqc_",
    TaxonomicClassification = "^taxonomy_",
    Assembly2               = "^assembly_",
    MetagenomicReadMapping  = "^resistoma_",
    Annotation              = "^annotation_",
    ComprehensiveGenomeAnalysis = "^cga_",
    CoreGenomeMLST          = "^mlst",
    PhylogeneticTree        = "^filogenia",
    SeqComparison           = "^comparativa"
  )

  ws_folder_tool <- function(name) {
    for (tool in names(WS_FOLDER_MAP))
      if (grepl(WS_FOLDER_MAP[[tool]], name, ignore.case = TRUE)) return(tool)
    NULL
  }

  # ── Botón: escanear workspace ──────────────────────────────
  observeEvent(input$btn_recover_bvbrc_ids, {
    tok  <- rv$bvbrc_token
    user <- rv$bvbrc_user
    if (is.null(tok) || is.null(user)) {
      showNotification("Inicia sesión BV-BRC primero.", type = "warning"); return()
    }

    rv$ws_scan_dirs  <- NULL
    rv$ws_scan_error <- NULL

    showModal(modalDialog(
      title     = tagList(icon("cloud-download-alt"),
                          " Recuperar resultados — Workspace BV-BRC"),
      size      = "l", easyClose = TRUE,
      div(
        div(class = "alert alert-info", style = "font-size:12px;padding:6px 10px;",
          icon("folder-open"), " Escaneando home (2 niveles): ",
          code(ws_home(user))),
        uiOutput("bvbrc_ws_scan_ui"),
        hr(),
        # Entrada manual: el usuario pega la ruta que ve en BV-BRC → My Jobs
        div(
          tags$b("Ruta manual"),
          tags$p(style = "font-size:11px;color:#888;margin:2px 0 6px;",
            "En BV-BRC: My Jobs → clic en el job → ver 'Output Folder' o la ",
            "columna 'Output'. Pega esa ruta aquí."),
          div(style = "display:flex;gap:6px;",
            textInput("bvbrc_manual_path", NULL,
              placeholder = paste0(ws_out(user), "   (ej: ", ws_out(user), "/fastqc_MUESTRA)"),
              width = "100%"),
            actionButton("btn_list_manual_path",
              tagList(icon("search"), " Listar"),
              class = "btn-info btn-sm")
          ),
          uiOutput("bvbrc_manual_path_ui")
        )
      ),
      footer = tagList(
        modalButton("Cerrar"),
        actionButton("btn_import_all_ws",
          tagList(icon("download"), " Importar todos"),
          class = "btn-success")
      )
    ))

    # Asegurar que la carpeta PRUEBAAPP exista antes de escanear
    tryCatch(bvbrc_mkdir(tok, ws_out(user)), error = function(e) NULL)

    all_dirs   <- list()
    seen_paths <- character(0)

    add_if_result <- function(f) {
      # BV-BRC usa "folder", "directory" Y "job_result" para carpetas de salida
      if (!f$type %in% c("folder","directory","job_result")) return()
      if (!is.null(ws_folder_tool(f$name)) &&
          !(f$path %||% f$name) %in% seen_paths) {
        all_dirs  <<- c(all_dirs, list(f))
        seen_paths <<- c(seen_paths, f$path %||% f$name)
      }
    }

    scan_errors <- character(0)

    # ── Scan 1: carpeta de salida PRUEBAAPP ──────────────────────────────────
    out_dir_path <- ws_out(user)
    alog("INFO", paste0("BV-BRC scan 1 (PRUEBAAPP): ", out_dir_path))
    scan_out <- tryCatch(bvbrc_ls(tok, out_dir_path),
                         error = function(e) list(success=FALSE, error=conditionMessage(e)))
    if (isTRUE(scan_out$success)) {
      for (f in scan_out$files) add_if_result(f)
      alog("INFO", paste0("BV-BRC scan PRUEBAAPP OK: ", length(scan_out$files), " entradas"))
    } else {
      scan_errors <- c(scan_errors,
        paste0("PRUEBAAPP: ", scan_out$error %||% "error desconocido"))
      alog("WARN", paste0("BV-BRC scan PRUEBAAPP fallido: ", out_dir_path),
           scan_out$error %||% "?")
    }

    # ── Scan 2: /home nivel-1 y nivel-2 ─────────────────────────────────────
    skip_dirs <- c(".fastq_uploads")
    home_path <- ws_home(user)
    alog("INFO", paste0("BV-BRC scan 2 (home): ", home_path))
    scan_lvl1 <- tryCatch(bvbrc_ls(tok, home_path),
                          error = function(e) list(success=FALSE, error=conditionMessage(e)))
    if (isTRUE(scan_lvl1$success)) {
      for (f in scan_lvl1$files) {
        add_if_result(f)
        if (f$type %in% c("folder","directory") && !f$name %in% skip_dirs) {
          sub_scan <- tryCatch(bvbrc_ls(tok, f$path),
                               error = function(e) list(success = FALSE))
          if (!isTRUE(sub_scan$success)) next
          for (sf in sub_scan$files) add_if_result(sf)
        }
      }
      alog("INFO", paste0("BV-BRC scan home OK: ", length(scan_lvl1$files), " entradas"))
    } else {
      scan_errors <- c(scan_errors,
        paste0("home: ", scan_lvl1$error %||% "error desconocido"))
      alog("WARN", "BV-BRC scan home fallido", scan_lvl1$error %||% "?")
    }

    # ── Scan 3: rutas directas desde jobs BV-BRC ya registrados ─────────────
    # bvbrc_ls() ya intenta ambos dominios (@patricbrc.org y @bvbrc) internamente
    for (jj in isolate(rv$jobs)) {
      if (!identical(jj$platform, "bvbrc")) next
      if (is.null(jj$output_path) || is.null(jj$output_file)) next
      dir_path <- paste0(jj$output_path, "/", jj$output_file)
      if (dir_path %in% seen_paths) next
      if (is.null(ws_folder_tool(jj$output_file))) next
      ls_j <- tryCatch(bvbrc_ls(tok, dir_path), error = function(e) list(success=FALSE))
      if (isTRUE(ls_j$success) && length(ls_j$files) > 0) {
        synthetic <- list(name = jj$output_file, type = "folder",
                          size = 0, path = ls_j$path %||% dir_path)
        all_dirs   <- c(all_dirs, list(synthetic))
        seen_paths <- c(seen_paths, ls_j$path %||% dir_path)
        alog("INFO", paste0("BV-BRC scan jobs: ", dir_path),
             paste0(length(ls_j$files), " archivo(s)"))
      }
    }

    if (length(all_dirs) == 0 && length(scan_errors) > 0) {
      rv$ws_scan_error <- paste0(
        "No se encontraron resultados. Errores del workspace: ",
        paste(scan_errors, collapse = " | "),
        " — El token puede estar expirado. Cierra sesión y vuelve a iniciar.")
      alog("ERROR", "BV-BRC workspace scan sin resultados", rv$ws_scan_error)
    }

    rv$ws_scan_dirs <- all_dirs
    alog("INFO", paste0("BV-BRC workspace scan total: ", length(all_dirs),
                        " carpeta(s) de resultados"))
  })

  # ── Observer: listar ruta manual ───────────────────────────
  observeEvent(input$btn_list_manual_path, {
    tok  <- isolate(rv$bvbrc_token)
    path <- trimws(isolate(input$bvbrc_manual_path) %||% "")
    if (!nchar(path)) return()

    # Caso A: la ruta ES una carpeta de resultados
    folder_name <- basename(path)
    if (!is.null(ws_folder_tool(folder_name))) {
      synthetic <- list(name = folder_name, type = "directory", size = 0,
                        path = path)
      exist_paths <- vapply(rv$ws_scan_dirs %||% list(),
        function(f) f$path %||% "", character(1))
      if (!path %in% exist_paths)
        rv$ws_scan_dirs <- c(rv$ws_scan_dirs %||% list(), list(synthetic))
      output$bvbrc_manual_path_ui <- renderUI(
        div(class = "alert alert-success", style = "font-size:11px;padding:5px 8px;",
          icon("check"), paste0(" Carpeta '", folder_name, "' añadida.")))
      return()
    }

    # Caso B: la ruta es una carpeta padre — listar su contenido
    ls_res <- tryCatch(
      bvbrc_ls(tok, path),
      error = function(e) list(success = FALSE, error = conditionMessage(e))
    )
    if (!isTRUE(ls_res$success)) {
      output$bvbrc_manual_path_ui <- renderUI(
        div(class = "alert alert-danger", style = "font-size:11px;padding:5px 8px;",
          icon("times-circle"), " Error: ", ls_res$error %||% "no se pudo listar"))
      return()
    }

    new_dirs <- Filter(function(f)
      f$type %in% c("folder","directory") && !is.null(ws_folder_tool(f$name)),
      ls_res$files)

    if (!length(new_dirs)) {
      output$bvbrc_manual_path_ui <- renderUI(
        div(class = "alert alert-warning", style = "font-size:11px;padding:5px 8px;",
          icon("exclamation-triangle"),
          " No se encontraron carpetas de resultados en esa ruta. ",
          "Las carpetas deben llamarse fastqc_*, taxonomy_*, assembly_*, ",
          "resistoma_*, annotation_*, cga_*, mlst*, filogenia*"))
      return()
    }

    exist_paths <- vapply(rv$ws_scan_dirs %||% list(),
      function(f) f$path %||% "", character(1))
    to_add <- Filter(function(f) !(f$path %||% "") %in% exist_paths, new_dirs)
    rv$ws_scan_dirs <- c(rv$ws_scan_dirs %||% list(), to_add)

    output$bvbrc_manual_path_ui <- renderUI(
      div(class = "alert alert-success", style = "font-size:11px;padding:5px 8px;",
        icon("check"), paste0(" ", length(to_add),
          " carpeta(s) añadida(s). Haz clic en 'Importar todos'.")))
  })

  # ── Render: tabla de carpetas encontradas ──────────────────
  output$bvbrc_ws_scan_ui <- renderUI({
    if (!is.null(rv$ws_scan_error))
      return(div(class = "alert alert-danger",
        icon("times-circle"), " Error: ", rv$ws_scan_error))

    dirs <- rv$ws_scan_dirs
    if (is.null(dirs))
      return(div(style = "padding:10px;color:#888;",
        icon("spinner", class = "fa-spin"), " Escaneando workspace..."))
    if (!length(dirs))
      return(div(class = "alert alert-warning",
        icon("exclamation-triangle"),
        " No se encontraron carpetas de resultados en home ni en sus subcarpetas. ",
        "Usa la 'Ruta manual' abajo: en BV-BRC → My Jobs → clic en el job → ",
        "copia el valor de 'Output Folder' y pégalo aquí."))

    rows <- lapply(dirs, function(f) {
      tool <- ws_folder_tool(f$name) %||% "?"
      step <- TOOL_TO_STEP[tool] %||% "?"
      ya <- any(vapply(rv$jobs, function(j)
        isTRUE(j$platform == "bvbrc") &&
        isTRUE((j$output_file %||% "") == f$name) &&
        !is.null(rv$job_results[[j$job_id]]),
        logical(1)))
      tags$tr(
        tags$td(code(f$name)),
        tags$td(tags$small(tool)),
        tags$td(tags$small(step)),
        tags$td(if (ya)
          span(style = "color:green;font-size:11px;", icon("check"), " importado")
          else
          span(style = "color:#aaa;font-size:11px;",  icon("clock"),  " pendiente")
        )
      )
    })

    tagList(
      div(class = "alert alert-success", style = "font-size:12px;padding:6px 10px;",
        icon("check"), paste0(" ", length(dirs),
          " carpeta(s) encontrada(s). Haz clic en 'Importar todos'.")),
      tags$table(class = "table table-condensed table-hover",
        style = "font-size:12px;margin-bottom:4px;",
        tags$thead(tags$tr(
          tags$th("Carpeta"), tags$th("Herramienta"),
          tags$th("Etapa"), tags$th("Estado"))),
        tags$tbody(rows)),
      uiOutput("bvbrc_ws_import_status")
    )
  })

  # ── Botón: importar todas las carpetas encontradas ─────────
  observeEvent(input$btn_import_all_ws, {
    tok  <- isolate(rv$bvbrc_token)
    user <- isolate(rv$bvbrc_user)
    dirs <- isolate(rv$ws_scan_dirs)
    if (is.null(tok) || !length(dirs)) {
      showNotification("No hay carpetas para importar.", type = "warning"); return()
    }

    n_ok <- 0L; n_err <- 0L; msgs <- character(0)

    withProgress(message = "Importando resultados BV-BRC...", value = 0, {
      for (i in seq_along(dirs)) {
        f <- dirs[[i]]
        incProgress(1 / length(dirs), detail = f$name)

        tool_name <- ws_folder_tool(f$name) %||% "Unknown"
        step_id   <- TOOL_TO_STEP[tool_name] %||% "inicio"
        prefix    <- WS_FOLDER_MAP[[tool_name]] %||% "^[^_]+_"
        muestra   <- sub(prefix, "", f$name)

        job_id <- paste0("ws_", f$name)
        # output_path derivado de f$path real (puede estar en home/ o en home/sub/)
        f_output_path <- if (!is.null(f$path) && nchar(f$path) > nchar(f$name))
          sub(paste0("/", f$name, "$"), "", f$path)
        else
          ws_home(user)
        job <- list(
          job_id      = job_id,
          platform    = "bvbrc",
          tool        = tool_name,
          step        = step_id,
          label       = paste0(tool_name, " — ", muestra),
          status      = "completed",
          time        = Sys.time(),
          output_path = f_output_path,
          output_file = f$name,
          muestra     = muestra,
          output_ids  = NULL
        )

        # Registrar SIEMPRE el job como completado (aunque el import falle):
        # el folder existe → la muestra ya se procesó → no re-analizar.
        ei <- which(vapply(rv$jobs, function(jj) jj$job_id == job_id, logical(1)))
        if (length(ei)) rv$jobs[[ei[1]]] <- job else rv$jobs <- c(rv$jobs, list(job))

        fetch_res <- tryCatch(
          bvbrc_fetch_results(tok, job),
          error = function(e) list(success = FALSE, error = conditionMessage(e))
        )

        if (!isTRUE(fetch_res$success)) {
          msgs  <- c(msgs, paste0("⚠️ ", f$name, ": registrado (sin re-análisis), ",
                                  "pero no se pudieron leer los resultados: ",
                                  fetch_res$error %||% "fetch fallido"))
          n_err <- n_err + 1L
          alog("WARN", paste0("BV-BRC ws registrado sin resultados: ", f$name),
               fetch_res$error %||% "")
          next
        }

        rv$job_results[[job_id]] <- fetch_res
        auto_import_results(fetch_res, job)

        msgs <- c(msgs, paste0("✅ ", f$name, " → ",
                               fetch_res$parsed$type %||% "?", " (",
                               length(fetch_res$files), " archivos)"))
        alog("OK", paste0("BV-BRC ws importado: ", f$name),
             paste0("tipo: ", fetch_res$parsed$type, " | archivos: ",
                    length(fetch_res$files)))
        n_ok <- n_ok + 1L
      }
    })

    tryCatch(state_save(rv), error = function(e) NULL)

    output$bvbrc_ws_import_status <- renderUI(div(
      style = "margin-top:8px;",
      if (n_ok > 0)
        div(class = "alert alert-success", style = "padding:6px 10px;",
          icon("check-circle"), paste0(" ", n_ok, " carpeta(s) importadas.")),
      if (n_err > 0)
        div(class = "alert alert-danger", style = "padding:6px 10px;",
          icon("times-circle"), paste0(" ", n_err, " con errores.")),
      tags$ul(style = "font-size:11px;margin:6px 0 0 0;padding-left:16px;",
        lapply(msgs, tags$li))
    ))

    if (n_ok > 0)
      showNotification(
        paste0(n_ok, " resultado(s) importados del workspace BV-BRC."),
        type = "message", duration = 7)
  })

  output$bvbrc_ws_scan_ui       <- renderUI(NULL)
  output$bvbrc_ws_import_status <- renderUI(NULL)

}
