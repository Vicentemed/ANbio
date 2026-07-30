# ============================================================
# Helpers de UI reutilizables
# ============================================================

# Tarjeta de herramienta con boton de apertura
tool_card <- function(nombre, descripcion, link, entrada) {
  div(style = "background:#f0f4f8;border-left:4px solid #3c8dbc;padding:10px 14px;margin-bottom:8px;border-radius:4px;",
    div(
      tags$a(href = link, target = "_blank",
             style = "font-weight:bold;font-size:14px;color:#1a5276;",
             icon("external-link-alt"), " ", nombre)
    ),
    p(descripcion, style = "color:#555;font-size:12.5px;margin:5px 0 3px;"),
    tags$span(paste("Entrada:", entrada),
              style = "background:#d6eaf8;color:#1a5276;padding:2px 8px;border-radius:8px;font-size:11px;font-weight:bold;")
  )
}

criterios_panel <- function(criterios) {
  div(style = "background:#e8f4f8;border:1px solid #bee5eb;padding:10px 14px;border-radius:4px;margin-bottom:12px;",
    tags$b(icon("clipboard-check"), " Criterios de calidad:"),
    tags$ul(style = "margin:5px 0 0 0;padding-left:18px;",
      lapply(criterios, function(c) tags$li(style = "font-size:13px;", c))
    )
  )
}

paso_header <- function(numero, titulo, descripcion) {
  div(style = "background:linear-gradient(135deg,#1a5276,#2e86c1);color:white;padding:14px 20px;border-radius:6px;margin-bottom:18px;",
    tags$h3(style = "color:white;margin:0 0 5px;", paste0("Paso ", numero, ": ", titulo)),
    tags$p(style = "color:#d6eaf8;margin:0;font-size:13px;", descripcion)
  )
}

# Panel de estado de un job
job_status_ui <- function(job) {
  if (is.null(job)) return(NULL)

  icon_map <- list(
    queued    = list(i = "clock",         col = "#7f8c8d", cls = "info"),
    running   = list(i = "spinner",       col = "#f39c12", cls = "warning"),
    RUNNING   = list(i = "spinner",       col = "#f39c12", cls = "warning"),
    FINISHED  = list(i = "check-circle",  col = "#27ae60", cls = "success"),
    completed = list(i = "check-circle",  col = "#27ae60", cls = "success"),
    ok        = list(i = "check-circle",  col = "#27ae60", cls = "success"),
    error     = list(i = "times-circle",  col = "#e74c3c", cls = "danger"),
    ERROR     = list(i = "times-circle",  col = "#e74c3c", cls = "danger"),
    failed    = list(i = "times-circle",  col = "#e74c3c", cls = "danger")
  )
  s <- icon_map[[job$status]] %||% list(i = "question-circle", col = "#888", cls = "secondary")

  div(class = paste0("alert alert-", s$cls),
    style = "margin:8px 0 0;padding:8px 12px;",
    if (job$status %in% c("running", "RUNNING", "queued"))
      icon(s$i, class = "fa-spin")
    else
      icon(s$i),
    strong(paste0(" [", job$platform, "] ", job$label, " — ")),
    job$status,
    if (!is.null(job$job_id))
      tags$small(style = "color:#666;margin-left:8px;", paste0("ID: ", substr(job$job_id, 1, 20)))
  )
}

# Selector de muestra a partir del vector rv$muestras
muestra_selector <- function(inputId, label = "Muestra:", choices = character(0)) {
  if (length(choices) == 0) choices <- c("(Sube archivos FASTQ primero)" = "")
  selectInput(inputId, label, choices = choices, width = "100%")
}

# Badge de plataforma
platform_badge <- function(platform) {
  colors <- list(
    bvbrc  = "#c0392b",
    galaxy = "#2980b9",
    ebi    = "#16a085",
    ncbi   = "#8e44ad",
    cge    = "#d35400",
    web    = "#7f8c8d"
  )
  col <- colors[[tolower(platform)]] %||% "#555"
  tags$span(
    style = paste0("background:", col, ";color:white;padding:2px 8px;border-radius:4px;font-size:11px;font-weight:bold;margin-right:4px;"),
    toupper(platform)
  )
}

# Panel de ejecucion por herramienta
exec_panel <- function(step_id, tools_list, muestras = character(0)) {
  div(
    style = "background:#fffbea;border:1px solid #f0d060;border-radius:6px;padding:12px 14px;margin-top:10px;",
    div(style = "font-weight:bold;font-size:13px;color:#856404;margin-bottom:8px;",
        icon("bolt"), " Ejecucion directa via API"),
    fluidRow(
      column(4,
        if (length(muestras) > 0)
          selectInput(paste0("exec_muestra_", step_id), "Muestra:",
                      choices = muestras, width = "100%")
        else
          div(class = "alert alert-warning", style = "padding:5px 8px;font-size:12px;margin:0;",
              icon("exclamation-triangle"), " Sube archivos FASTQ primero")
      ),
      column(5,
        selectInput(paste0("exec_tool_", step_id), "Herramienta:",
                    choices = tools_list, width = "100%")
      ),
      column(3,
        br(),
        actionButton(paste0("exec_btn_", step_id),
                     tagList(icon("play"), " Ejecutar"),
                     class = "btn-warning btn-sm",
                     style = "width:100%;margin-top:5px;")
      )
    ),
    uiOutput(paste0("exec_status_", step_id)),
    uiOutput(paste0("exec_result_", step_id))
  )
}
