# ============================================================
# Logger — ANbio  (sin dependencias Shiny)
# Escribe a archivo Y devuelve la entrada para rv$log_entries
# ============================================================

LOG_FILE <- file.path(tempdir(), "anbio_session_log.txt")

# Inicializar archivo al cargar el modulo
.log_init <- function() {
  tryCatch(
    cat(paste0("# ── ANbio log ── ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               " ──────────────────────────────\n"),
        file = LOG_FILE, append = FALSE),
    error = function(e) NULL
  )
}

# ── Funcion principal ──────────────────────────────────────
app_log <- function(nivel, mensaje, detalle = NULL) {
  ts    <- format(Sys.time(), "%H:%M:%S")
  fecha <- format(Sys.time(), "%Y-%m-%d")
  niv   <- toupper(as.character(nivel))

  # Construir linea de archivo
  linea <- paste0("[", fecha, " ", ts, "] [", niv, "] ", mensaje)
  if (!is.null(detalle) && nchar(trimws(as.character(detalle))) > 0) {
    det_str <- gsub("\n", "\n    | ", trimws(as.character(detalle)))
    linea   <- paste0(linea, "\n    | ", det_str)
  }
  tryCatch(
    cat(linea, "\n", file = LOG_FILE, append = TRUE, sep = ""),
    error = function(e) NULL
  )

  # Devuelve lista para rv$log_entries
  list(
    ts      = ts,
    fecha   = fecha,
    nivel   = niv,
    mensaje = as.character(mensaje),
    detalle = if (!is.null(detalle)) substr(as.character(detalle), 1, 600) else ""
  )
}

# ── Shortcuts ──────────────────────────────────────────────
log_info  <- function(msg, det = NULL) app_log("INFO",  msg, det)
log_ok    <- function(msg, det = NULL) app_log("OK",    msg, det)
log_warn  <- function(msg, det = NULL) app_log("WARN",  msg, det)
log_error <- function(msg, det = NULL) app_log("ERROR", msg, det)
log_api   <- function(msg, det = NULL) app_log("API",   msg, det)
log_debug <- function(msg, det = NULL) app_log("DEBUG", msg, det)

# ── Leer archivo (ultimas N lineas) ───────────────────────
log_tail <- function(n = 400) {
  if (!file.exists(LOG_FILE)) return("(sin entradas)")
  ll <- readLines(LOG_FILE, warn = FALSE)
  if (length(ll) == 0) return("(sin entradas)")
  paste(tail(ll, n), collapse = "\n")
}

# ── Borrar log ────────────────────────────────────────────
log_clear <- function() {
  tryCatch(
    cat(paste0("# ── Log limpiado: ", format(Sys.time(), "%H:%M:%S"), "\n"),
        file = LOG_FILE, append = FALSE),
    error = function(e) NULL
  )
}

.log_init()
