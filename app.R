# ============================================================
# ANÁLISIS BIOINFORMÁTICO — LESP Aguascalientes
# App Shiny: interfaz unificada para herramientas NGS en línea
# ============================================================

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# ---- Pandoc ----
.setup_pandoc <- function() {
  if (rmarkdown::pandoc_available()) return(invisible(TRUE))
  candidatos <- unique(c(
    Sys.getenv("RSTUDIO_PANDOC"),
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "C:/Program Files/RStudio/bin/pandoc",
    file.path(Sys.getenv("LOCALAPPDATA"), "Pandoc"),
    dirname(Sys.which("pandoc"))
  ))
  for (dir in candidatos[nchar(candidatos) > 0]) {
    exe <- file.path(dir, if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc")
    if (file.exists(exe)) {
      Sys.setenv(RSTUDIO_PANDOC = dir)
      rmarkdown::find_pandoc(cache = FALSE, dir = dir)
      return(invisible(TRUE))
    }
  }
  invisible(FALSE)
}
PANDOC_OK <- .setup_pandoc()

# ---- Librerias ----
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyFiles)
library(DT)
library(ggplot2)
library(plotly)
library(dplyr)
library(rmarkdown)
library(httr2)
library(jsonlite)

# ---- Modulos de API ----
source("R/logger.R",          local = TRUE)
source("R/api_bvbrc.R",       local = TRUE)
source("R/api_ebi.R",         local = TRUE)
source("R/api_ncbi.R",        local = TRUE)
source("R/api_galaxy.R",      local = TRUE)
source("R/api_cge.R",         local = TRUE)
source("R/parsers_results.R",  local = TRUE)
source("R/session_state.R",   local = TRUE)
source("R/tree_plot.R",        local = TRUE)
source("R/ui_helpers.R",       local = TRUE)

# Directorio temporal de sesion para archivos subidos
SESSION_TMPDIR <- file.path(tempdir(), "anbio_session")
dir.create(SESSION_TMPDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# DATOS GLOBALES: herramientas por etapa (fuente: PDF LESP 2026)
# ============================================================
herramientas_info <- list(
  calidad = list(
    numero = 1, titulo = "Control de Calidad", icono = "check-circle", color = "green",
    descripcion = "Evaluación de calidad de lecturas FASTQ: scores Phred, contenido GC, duplicados y adaptadores.",
    criterios = c("% bases ≥ Q30 ≥ 75-85%","Clusters Passing Filter ≥ 80%","Error Rate < 1%","% PhiX 1-5%"),
    herramientas = data.frame(
      Herramienta = "Fastq Utilities (BV-BRC)",
      Descripcion = "Evaluación de calidad tipo FastQC y trimming de lecturas",
      Link        = "https://www.bv-brc.org/app/FastqUtil",
      Entrada     = "FASTQ (.fastq.gz)",
      Plataforma  = "bvbrc",
      AppID       = "FastqUtils",
      stringsAsFactors = FALSE
    )
  ),
  taxonomia = list(
    numero = 2, titulo = "Análisis Taxonómico", icono = "leaf", color = "olive",
    descripcion = "Clasificación taxonómica de secuencias para confirmar identidad del organismo y detectar contaminación.",
    criterios = c("Contaminación < 5%","Identificación a nivel de especie con confianza alta"),
    herramientas = data.frame(
      Herramienta = c("Genome Distance (BV-BRC)","Taxonomic Classification (BV-BRC)"),
      Descripcion = c("Comparación rápida entre genomas usando distancias MinHash","Clasificación taxonómica por k-mers (Kraken2)"),
      Link        = c("https://www.bv-brc.org/app/GenomeDistance","https://www.bv-brc.org/app/TaxonomicClassification"),
      Entrada     = c("FASTQ/FASTA","FASTQ"),
      Plataforma  = c("bvbrc","bvbrc"),
      AppID       = c("GenomeDistance","TaxonomicClassification"),
      stringsAsFactors = FALSE
    )
  ),
  ensamblado = list(
    numero = 3, titulo = "Ensamblado", icono = "puzzle-piece", color = "blue",
    descripcion = "Ensamblado de novo del genoma a partir de lecturas FASTQ usando grafos de Bruijn.",
    criterios = c("Cobertura ≥ 95%","Profundidad ≥ 30X","Contigs < 200","N50 > 50 kb","Bases ambiguas mínimas"),
    herramientas = data.frame(
      Herramienta = c("SPAdes (Galaxy)","Unicycler (Galaxy)","Assembly Service (BV-BRC)"),
      Descripcion = c("Ensamblador basado en grafos de Bruijn","Ensamblado híbrido optimizado para bacterias","Ensamblado de novo con trimming y polishing integrados"),
      Link        = c(
        "https://galaxy-main.usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fnml%2Fspades%2Fspades%2F4.2.0%2Bgalaxy0&version=latest",
        "https://galaxy-main.usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Funicycler%2Funicycler%2F0.5.1%2Bgalaxy0&version=latest",
        "https://www.bv-brc.org/app/Assembly2"
      ),
      Entrada    = c("FASTQ","FASTQ","FASTQ"),
      Plataforma = c("galaxy","galaxy","bvbrc"),
      AppID      = c(GALAXY_TOOL_SPADES, GALAXY_TOOL_UNICYCLER, "Assembly2"),
      stringsAsFactors = FALSE
    )
  ),
  anotacion = list(
    numero = 4, titulo = "Anotación", icono = "tags", color = "yellow",
    descripcion = "Identificación y anotación funcional de genes usando RASTtk (bacterias) y VIGOR4 (virus).",
    criterios = c("Cobertura de anotación > 80%","Proporción aceptable genes hipotéticos vs funcionales"),
    herramientas = data.frame(
      Herramienta = c("Comprehensive Genome Analysis (BV-BRC)","Annotation Service (BV-BRC)"),
      Descripcion = c("Anotación integral: genes estructurales, funcionales y elementos móviles","Anotación automática con RASTtk/VIGOR4"),
      Link        = c("https://www.bv-brc.org/app/ComprehensiveGenomeAnalysis","https://www.bv-brc.org/app/Annotation"),
      Entrada     = c("FASTQ","FASTA (contigs)"),
      Plataforma  = c("bvbrc","bvbrc"),
      AppID       = c("ComprehensiveGenomeAnalysis","Annotation"),
      stringsAsFactors = FALSE
    )
  ),
  resistoma = list(
    numero = 5, titulo = "Resistoma", icono = "shield-alt", color = "red",
    descripcion = "Identificación de genes AMR y clasificación MDR/XDR/PDR del perfil de resistencia.",
    criterios = c("MDR: > 3 grupos","XDR: > 10 grupos","PDR: todos los grupos"),
    herramientas = data.frame(
      Herramienta = c("Comprehensive Genome Analysis (BV-BRC)","Metagenomic Read Mapping (BV-BRC)","ResFinder (DTU/CGE)","WHO GLASS","Mykrobe"),
      Descripcion = c("Identificación integral de genes AMR en genoma completo","Detección de genes AMR desde lecturas sin ensamblado","Identificación de genes AMR (API pública DTU)","Tendencias globales de resistencia","Predicción AMR desde lecturas o genomas"),
      Link        = c(
        "https://www.bv-brc.org/app/ComprehensiveGenomeAnalysis",
        "https://www.bv-brc.org/app/MetagenomicReadMapping",
        "https://genepi.food.dtu.dk/resfinder",
        "https://worldhealthorg.shinyapps.io/glassdashboard/#!/amr",
        "https://www.mykrobe.com/"
      ),
      Entrada    = c("FASTQ","FASTQ","FASTA/FASTQ","Web","FASTQ/FASTA"),
      Plataforma = c("bvbrc","bvbrc","cge","web","web"),
      AppID      = c("ComprehensiveGenomeAnalysis","MetagenomicReadMapping","resfinder","",""),
      stringsAsFactors = FALSE
    )
  ),
  alineamiento = list(
    numero = 6, titulo = "Alineamiento", icono = "align-left", color = "purple",
    descripcion = "Alineamiento de secuencias para detectar variaciones, mutaciones y similitudes.",
    criterios = c("Identidad > 95% misma especie","E-value < 1e-5","Cobertura ≥ 90%"),
    herramientas = data.frame(
      Herramienta = c("BLAST (BV-BRC)","MUSCLE (EBI JDispatcher)","Clustal Omega (EBI)","MAFFT (EBI)","COBALT (NCBI)","T-Coffee","Expresso (T-Coffee)"),
      Descripcion = c("Búsqueda similitud entre secuencias","Alineamiento múltiple de secuencias (API pública)","Alineamiento múltiple rápido (API pública)","Alineamiento múltiple preciso (API pública)","Alineamiento por dominios conservados","Alineamiento múltiple de alta precisión","Alineamiento basado en estructura proteica"),
      Link        = c(
        "https://www.bv-brc.org/app/Homology",
        "https://www.ebi.ac.uk/jdispatcher/msa/muscle",
        "https://www.ebi.ac.uk/jdispatcher/msa/clustalo",
        "https://mafft.cbrc.jp/alignment/server/",
        "https://www.ncbi.nlm.nih.gov/tools/cobalt/",
        "https://tcoffee.crg.eu/",
        "https://tcoffee.crg.eu/apps/tcoffee/do:expresso"
      ),
      Entrada    = c("FASTA","FASTA","FASTA","FASTA","FASTA","FASTA","FASTA"),
      Plataforma = c("bvbrc","ebi","ebi","ebi","web","web","web"),
      AppID      = c("Homology","muscle","clustalo","mafft","","",""),
      stringsAsFactors = FALSE
    )
  ),
  mlst = list(
    numero = 7, titulo = "MLST", icono = "barcode", color = "teal",
    descripcion = "Tipificación multilocus (MLST) para asignación de ST y rastreo epidemiológico.",
    criterios = c("Asignación ST con todos los loci cubiertos","Cobertura de cada locus ≥ 95%"),
    herramientas = data.frame(
      Herramienta = c("Core Genome MLST (BV-BRC)","MLST (Galaxy)"),
      Descripcion = c("Tipificación cgMLST basada en genes core","Asignación de ST usando esquemas PubMLST/Pasteur"),
      Link        = c(
        "https://www.bv-brc.org/app/CoreGenomeMLST",
        "https://galaxy-main.usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fmlst%2Fmlst%2F2.22.0&version=latest"
      ),
      Entrada    = c("FASTA","FASTA"),
      Plataforma = c("bvbrc","galaxy"),
      AppID      = c("CoreGenomeMLST", GALAXY_TOOL_MLST),
      stringsAsFactors = FALSE
    )
  ),
  filogenia = list(
    numero = 8, titulo = "Filogenia", icono = "code-branch", color = "maroon",
    descripcion = "Construcción de árboles filogenéticos para análisis evolutivo y rastreo de brotes.",
    criterios = c("Bootstrap ≥ 70%","Modelo evolutivo apropiado (GTR+G, HKY+G)"),
    herramientas = data.frame(
      Herramienta = c("Phylogenetic Tree (BV-BRC)","IQ-TREE (Galaxy)","iTOL (EMBL)","Phylogeny.fr"),
      Descripcion = c("Árbol filogenético integrado en BV-BRC","Máxima verosimilitud con selección automática de modelo","Visualización y anotación de árboles (subir árbol Newick)","Plataforma en línea para construir árboles"),
      Link        = c(
        "https://www.bv-brc.org/app/PhylogeneticTree",
        "https://galaxy-main.usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fiqtree%2Fiqtree%2F2.4.0%2Bgalaxy1&version=latest",
        "https://itol.embl.de/",
        "http://www.phylogeny.fr/"
      ),
      Entrada    = c("FASTA/Genomas","FASTA/MSA","Newick","FASTA/Alineamiento"),
      Plataforma = c("bvbrc","galaxy","web","web"),
      AppID      = c("PhylogeneticTree", GALAXY_TOOL_IQTREE, "", ""),
      stringsAsFactors = FALSE
    )
  ),
  genomica = list(
    numero = 9, titulo = "Genómica Comparativa", icono = "dna", color = "navy",
    descripcion = "Comparación de genomas: ANI, pangenoma y redes funcionales entre aislamientos.",
    criterios = c("ANI ≥ 95% = misma especie","ANI ≥ 99% = mismo clon/serotipo"),
    herramientas = data.frame(
      Herramienta = c("ANI (EZBioCloud)","Proteome Comparison (BV-BRC)","Proteomica (KBase)","Cytoscape (Galaxy)"),
      Descripcion = c("Average Nucleotide Identity entre aislamientos","Comparación de proteomas entre múltiples genomas","Comparación de proteomas (genes compartidos/diferencias)","Visualización de redes biológicas"),
      Link        = c(
        "https://www.ezbiocloud.net/tools/ani",
        "https://www.bv-brc.org/app/SeqComparison",
        "https://narrative.kbase.us/legacy/catalog/apps/GenomeProteomeComparison/compare_two_proteomes/release",
        "https://galaxy-main.usegalaxy.org/visualizations/create/cytoscape"
      ),
      Entrada    = c("FASTA","FASTA","FASTA/Genomas KBase","Redes JSON/SIF"),
      Plataforma = c("web","bvbrc","web","web"),
      AppID      = c("","SeqComparison","",""),
      stringsAsFactors = FALSE
    )
  )
)

PASOS_IDS    <- c("calidad","taxonomia","ensamblado","anotacion","resistoma","alineamiento","mlst","filogenia","genomica")
PASOS_LABELS <- c("Calidad","Taxonómica","Ensamblado","Anotación","Resistoma","Alineamiento","MLST","Filogenia","Gen. Comp.")
PASOS_ICONS  <- c("check-circle","leaf","puzzle-piece","tags","shield-alt","align-left","barcode","code-branch","dna")

# Leyenda al pie del reporte final. Es solo el valor POR DEFECTO:
# el analista puede cambiarla desde la pestaña "Generar Reporte" y
# el texto queda guardado en la sesión .anbio.
LEYENDA_REPORTE_DEFAULT <- paste(
  "Reporte Generado Automáticamente por ANbio-Vicente Esparza-Villalpando",
  "a partir de ejecución de datos en BV-BRC y Galaxy."
)

# Opciones de ejecucion por paso (plataforma -> label)
EXEC_TOOLS <- list(
  calidad     = c("FastqUtils via BV-BRC" = "bvbrc:FastqUtils"),
  taxonomia   = c("Kraken2 via Galaxy (recomendado)" = "galaxy:kraken2",
                  "TaxonomicClassification via BV-BRC" = "bvbrc:TaxonomicClassification",
                  "GenomeDistance via BV-BRC" = "bvbrc:GenomeDistance"),
  ensamblado  = c("SPAdes via Galaxy" = "galaxy:spades",
                  "Unicycler via Galaxy" = "galaxy:unicycler",
                  "Assembly2 via BV-BRC" = "bvbrc:Assembly2"),
  anotacion   = c("ComprehensiveGenomeAnalysis via BV-BRC" = "bvbrc:ComprehensiveGenomeAnalysis",
                  "Annotation via BV-BRC" = "bvbrc:Annotation"),
  resistoma   = c("ResFinder via CGE/DTU (publica)" = "cge:resfinder",
                  "MetagenomicReadMapping via BV-BRC" = "bvbrc:MetagenomicReadMapping",
                  "ComprehensiveGenomeAnalysis via BV-BRC" = "bvbrc:ComprehensiveGenomeAnalysis"),
  alineamiento = c("MUSCLE via EBI (publica)" = "ebi:muscle",
                   "Clustal Omega via EBI (publica)" = "ebi:clustalo",
                   "MAFFT via EBI (publica)" = "ebi:mafft",
                   "BLAST via NCBI (publica)" = "ncbi:blastn",
                   "BLAST via BV-BRC" = "bvbrc:Homology"),
  mlst        = c("CoreGenomeMLST via BV-BRC" = "bvbrc:CoreGenomeMLST",
                  "MLST via Galaxy" = "galaxy:mlst"),
  filogenia   = c("PhylogeneticTree via BV-BRC" = "bvbrc:PhylogeneticTree",
                  "IQ-TREE via Galaxy" = "galaxy:iqtree"),
  genomica    = c("SeqComparison via BV-BRC" = "bvbrc:SeqComparison")
)

# ---- Boton de navegacion con check ----
nav_buttons <- function(tab_id) {
  idx      <- which(PASOS_IDS == tab_id)
  prev_id  <- if (idx == 1) "inicio" else PASOS_IDS[idx - 1]
  next_id  <- if (idx == length(PASOS_IDS)) "reporte" else PASOS_IDS[idx + 1]
  prev_lbl <- if (idx == 1) "Inicio" else PASOS_LABELS[idx - 1]
  next_lbl <- if (idx == length(PASOS_IDS)) "Ver Reporte" else PASOS_LABELS[idx + 1]

  div(style = "margin-top:14px;",
    tags$hr(style = "margin:8px 0 12px;"),
    div(style = "display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:8px;padding:0 4px 6px;",
      actionButton(paste0("btn_ant_", tab_id),
        tags$span(icon("arrow-left"), " ", prev_lbl), class = "btn-default btn-sm"),
      div(style = "transform:scale(1.1);",
        checkboxInput(paste0("chk_", tab_id),
          tags$span(style = "font-size:13px;font-weight:bold;color:#27ae60;",
                    icon("check-circle"), " Etapa completada"), value = FALSE)
      ),
      actionButton(paste0("btn_sig_", tab_id),
        tags$span(next_lbl, " ", icon("arrow-right")), class = "btn-primary btn-sm")
    )
  )
}

# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = tags$span(tags$b("LESP — Análisis Bioinformático NGS"), style = "font-size:16px;"),
    titleWidth = 320,
    tags$li(class = "dropdown",
      actionLink("btn_credentials", tagList(icon("key"), " Credenciales"),
                 style = "padding:15px 12px;color:white;cursor:pointer;"))
  ),

  dashboardSidebar(
    width = 270,
    sidebarMenu(
      id = "sidebar",
      menuItem("Inicio / Muestras", tabName = "inicio", icon = icon("home")),
      menuItemOutput("sidebar_step_calidad"),
      menuItemOutput("sidebar_step_taxonomia"),
      menuItemOutput("sidebar_step_ensamblado"),
      menuItemOutput("sidebar_step_anotacion"),
      menuItemOutput("sidebar_step_resistoma"),
      menuItemOutput("sidebar_step_alineamiento"),
      menuItemOutput("sidebar_step_mlst"),
      menuItemOutput("sidebar_step_filogenia"),
      menuItemOutput("sidebar_step_genomica"),
      tags$hr(),
      uiOutput("sidebar_active_jobs"),
      tags$div(style = "padding:4px 14px 10px;", uiOutput("sidebar_progreso")),
      menuItem("Generar Reporte", tabName = "reporte", icon = icon("file-alt"),
               badgeLabel = "MD", badgeColor = "green"),
      menuItem("Log / Consola",  tabName = "log_panel", icon = icon("terminal"),
               badgeLabel = "DEV", badgeColor = "purple")
    )
  ),

  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML("
      .content-wrapper,.right-side{background:#f5f7fa}
      .metric-value{font-size:28px;font-weight:bold;color:#2c3e50;text-align:center}
      .metric-label{font-size:11px;color:#7f8c8d;text-align:center}
      .muestra-chip{display:inline-block;background:#d5e8d4;border:1px solid #82b366;padding:3px 10px;border-radius:10px;margin:3px;font-size:13px}
      #barra_progreso_global{position:sticky;top:0;z-index:999;background:white;border-bottom:2px solid #e0e7ef;padding:10px 20px 8px;margin-bottom:14px;box-shadow:0 2px 6px rgba(0,0,0,.06)}
      .step-dot{display:flex;flex-direction:column;align-items:center;min-width:58px;cursor:pointer}
      .step-circle{width:32px;height:32px;border-radius:50%;display:flex;align-items:center;justify-content:center;font-size:12px;font-weight:bold;transition:all .25s}
      .step-circle.pending{background:#e9ecef;color:#868e96}
      .step-circle.active{background:#2e86c1;color:white;box-shadow:0 0 0 4px #aed6f1}
      .step-circle.done{background:#27ae60;color:white}
      .step-label{font-size:10px;margin-top:3px;text-align:center;line-height:1.2}
      .step-label.active{font-weight:bold;color:#2e86c1}
      .step-label.done{color:#27ae60}
      .step-label.pending{color:#999}
      .step-connector{height:3px;flex:1;margin-bottom:18px;min-width:6px;border-radius:2px;transition:background .3s}
      .sidebar-prog-bar{height:5px;border-radius:3px;background:#2a4870;margin-top:4px;overflow:hidden}
      .sidebar-prog-fill{height:100%;background:#27ae60;border-radius:3px;transition:width .4s}
      .job-card{border-left:4px solid #f39c12;background:#fffbea;padding:8px 12px;border-radius:4px;margin-top:6px;font-size:13px}
      .job-done{border-color:#27ae60;background:#eafaf1}
      .job-error{border-color:#e74c3c;background:#fdf0ef}
      .sidebar-jobs-panel{padding:2px 14px 8px;border-top:1px solid rgba(255,255,255,.12);margin:0 0 2px}
      .sj-header{font-size:10px;color:#95a5a6;text-transform:uppercase;letter-spacing:.4px;padding:4px 0 3px;display:flex;align-items:center;gap:4px}
      .sj-row{display:flex;align-items:center;gap:5px;padding:3px 0;font-size:11px;color:#d5dbdb;line-height:1.3;border-bottom:1px solid rgba(255,255,255,.06)}
      .sj-row:last-child{border-bottom:none}
      .sj-icon{font-size:10px;flex-shrink:0;width:12px;text-align:center}
      .sj-label{flex:1;overflow:hidden;white-space:nowrap;text-overflow:ellipsis;max-width:140px}
      .sj-step{color:#85c1e9;font-weight:bold}
      .sj-time{color:#7f8c8d;font-size:10px;white-space:nowrap;margin-left:2px}
      .sj-done .sj-icon{color:#a9dfbf}
      .sj-done .sj-label{color:#a9dfbf}
      .sj-err  .sj-icon{color:#f1948a}
      .sj-err  .sj-label{color:#f1948a}
      .skin-blue .sidebar-menu>li>a>.badge{min-width:36px;text-align:center;font-size:10px;padding:2px 5px}
    "))),

    # Modal de credenciales
    tags$div(id = "credentials_modal_placeholder"),

    uiOutput("barra_progreso_global"),

    tabItems(

      # ============================================================
      # INICIO
      # ============================================================
      tabItem(tabName = "inicio",
        fluidRow(
          column(12,
            div(style = "background:linear-gradient(135deg,#1a5276,#2e86c1);color:white;padding:18px 24px;border-radius:8px;margin-bottom:20px;",
              tags$h2(style = "color:white;margin:0 0 6px;", icon("dna"), " Análisis Bioinformático de Patógenos"),
              tags$p(style = "color:#d6eaf8;margin:0;",
                "LESP Aguascalientes | Interfaz unificada para herramientas NGS: BV-BRC, Galaxy, EBI, NCBI, ResFinder")
            )
          )
        ),

        # ══════════ ASISTENTE DE 3 PASOS ══════════
        fluidRow(
          box(title = tagList(icon("play-circle"), " Iniciar análisis en 3 pasos"),
              status = "success", solidHeader = TRUE, width = 12,
            fluidRow(
              # Paso 1 — carpeta
              column(4,
                div(style = "border:2px solid #27ae60;border-radius:8px;padding:12px;min-height:190px;",
                  tags$h4(style = "margin:0 0 8px;color:#1e8449;",
                          tags$span(style="background:#27ae60;color:white;border-radius:50%;padding:2px 9px;margin-right:6px;", "1"),
                          "Carpeta con datos crudos"),
                  shinyDirButton("wiz_carpeta", "Seleccionar carpeta...",
                                 "Carpeta con los FASTQ", icon = icon("folder-open"),
                                 class = "btn-success", style = "width:100%;"),
                  br(), br(),
                  uiOutput("wiz_carpeta_info")
                )
              ),
              # Paso 2 — credenciales
              column(4,
                div(style = "border:2px solid #2e86c1;border-radius:8px;padding:12px;min-height:190px;",
                  tags$h4(style = "margin:0 0 8px;color:#1a5276;",
                          tags$span(style="background:#2e86c1;color:white;border-radius:50%;padding:2px 9px;margin-right:6px;", "2"),
                          "Credenciales"),
                  actionButton("wiz_creds", tagList(icon("key"), " Configurar BV-BRC y Galaxy"),
                               class = "btn-primary", style = "width:100%;"),
                  br(), br(),
                  uiOutput("wiz_creds_info")
                )
              ),
              # Paso 3 — iniciar
              column(4,
                div(style = "border:2px solid #e67e22;border-radius:8px;padding:12px;min-height:190px;",
                  tags$h4(style = "margin:0 0 8px;color:#ba4a00;",
                          tags$span(style="background:#e67e22;color:white;border-radius:50%;padding:2px 9px;margin-right:6px;", "3"),
                          "Ejecutar las 9 etapas"),
                  textInput("wiz_organismo", NULL, value = "Acinetobacter baumannii",
                            placeholder = "Organismo (nombre científico)"),
                  actionButton("wiz_iniciar",
                               tagList(icon("rocket"), " Iniciar análisis completo"),
                               class = "btn-warning btn-lg", style = "width:100%;font-weight:bold;"),
                  br(), br(),
                  uiOutput("wiz_estado")
                )
              )
            )
          )
        ),

        # ══════════ SESIÓN DE TRABAJO ══════════
        fluidRow(
          box(title = tagList(icon("save"), " Sesión de trabajo (reanudar análisis)"),
              status = "info", solidHeader = TRUE, width = 12, collapsible = TRUE,
            div(class = "alert alert-info", style = "padding:8px 12px;font-size:13px;",
              icon("info-circle"),
              " El análisis continúa en los servidores aunque cierres la app. La sesión ",
              tags$b("se guarda automáticamente"), " al cerrar; al volver a abrirla, carga el ",
              "archivo de sesión y las etapas se actualizarán solas con el avance real."
            ),
            fluidRow(
              column(4,
                textInput("sesion_nombre", "Nombre de la sesión:", value = "sesion_activa"),
                actionButton("btn_sesion_guardar", tagList(icon("save"), " Guardar ahora"),
                             class = "btn-info btn-sm", style = "width:100%;")
              ),
              column(4,
                tags$label("Sesiones guardadas:"),
                uiOutput("sesion_lista_ui"),
                actionButton("btn_sesion_cargar", tagList(icon("folder-open"), " Cargar y reanudar"),
                             class = "btn-success btn-sm", style = "width:100%;margin-top:6px;")
              ),
              column(4,
                tags$label("Archivo portable (.anbio):"),
                br(),
                downloadButton("dl_sesion", tagList(icon("download"), " Descargar sesión"),
                               class = "btn-default btn-sm", style = "width:100%;"),
                br(), br(),
                fileInput("up_sesion", NULL, accept = ".anbio",
                          buttonLabel = "Abrir archivo...", placeholder = "sesion.anbio")
              )
            ),
            uiOutput("sesion_estado_ui")
          )
        ),
        fluidRow(
          box(title = "Archivos FASTQ", status = "primary", solidHeader = TRUE, width = 5,
            div(class = "alert alert-info", style = "margin-bottom:10px;padding:8px 12px;font-size:13px;",
              icon("info-circle"),
              " Sube los archivos FASTQ aquí. La app los usará para enviarlos a las plataformas de análisis."
            ),
            fileInput("fastq_upload", "Seleccionar archivos FASTQ:",
                      multiple = TRUE,
                      accept = c(".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fasta", ".fa", ".fna"),
                      placeholder = "Seleccionar .fastq.gz o .fasta..."),
            uiOutput("muestras_detectadas"),
            hr(),
            div(style = "font-size:12px;color:#888;",
              icon("folder-open"),
              " También puedes seleccionar carpeta local si los archivos son muy grandes:",
              br(),
              shinyDirButton("carpeta_fastq", "Explorar carpeta...", "Seleccionar carpeta FASTQ",
                             icon = icon("folder"), class = "btn-default btn-sm")
            ),
            uiOutput("carpeta_seleccionada")
          ),
          box(title = "Información del proyecto", status = "info", solidHeader = TRUE, width = 7,
            fluidRow(
              column(6,
                textInput("proy_nombre",    "Proyecto",     placeholder = "Análisis ACINE 2025"),
                textInput("proy_organismo", "Organismo",    placeholder = "Acinetobacter baumannii"),
                textInput("proy_analista",  "Analista",     placeholder = "Nombre del analista")
              ),
              column(6,
                textInput("proy_institucion", "Institución", value = "LESP Aguascalientes"),
                dateInput("proy_fecha", "Fecha", value = Sys.Date()),
                selectInput("proy_tipo", "Tipo de muestra",
                  choices = c("Bacterias no tuberculosas","Mycobacterium tuberculosis","Virus respiratorios"))
              )
            )
          )
        ),
        fluidRow(
          box(title = "Estado de credenciales API", status = "warning", solidHeader = TRUE, width = 12,
            fluidRow(
              column(4,
                div(style = "text-align:center;padding:10px;",
                  uiOutput("cred_status_bvbrc"),
                  actionButton("btn_open_creds", tagList(icon("key"), " Configurar credenciales"),
                               class = "btn-warning btn-sm", style = "margin-top:8px;")
                )
              ),
              column(4,
                div(style = "text-align:center;padding:10px;",
                  uiOutput("cred_status_galaxy")
                )
              ),
              column(4,
                div(style = "background:#eaf4fb;border-radius:6px;padding:10px;font-size:12px;",
                  tags$b(icon("check-circle", style = "color:#27ae60;"), " APIs públicas (sin cuenta):"),
                  tags$ul(style = "margin:5px 0 0;padding-left:16px;",
                    tags$li("EBI JDispatcher — MUSCLE, Clustal Omega, MAFFT"),
                    tags$li("NCBI BLAST — búsqueda de similitud"),
                    tags$li("ResFinder (CGE/DTU) — genes AMR")
                  )
                )
              )
            )
          )
        ),
        fluidRow(
          box(title = "Parámetros de la corrida NGS (BaseSpace)", status = "warning", solidHeader = TRUE, width = 12,
            div(class = "alert alert-info", style = "margin-bottom:14px;",
              tags$strong(icon("info-circle"), " ¿De dónde vienen estos valores?"),
              " Estos parámetros son métricas reportadas por ",
              tags$strong("Illumina BaseSpace Sequence Hub"), " una vez que la corrida ha finalizado. Consultarlos en la plataforma y registrarlos aquí.",
              tags$br(),
              tags$a(href = "https://basespace.illumina.com", target = "_blank",
                icon("external-link-alt"), " Abrir Illumina BaseSpace Sequence Hub")
            ),
            fluidRow(
              column(2, numericInput("bs_cpf",  "Clusters PF (%)", value = NA, min=0, max=100, step=0.1)),
              column(2, numericInput("bs_q30",  "% ≥ Q30",         value = NA, min=0, max=100, step=0.1)),
              column(2, numericInput("bs_err",  "Error Rate (%)",  value = NA, min=0, max=10,  step=0.01)),
              column(2, numericInput("bs_phix", "% PhiX",          value = NA, min=0, max=10,  step=0.1)),
              column(4, br(), uiOutput("bs_status"))
            )
          )
        ),
        fluidRow(
          box(title = tagList(icon("rocket"), " Análisis por lotes (BV-BRC)"),
              status = "danger", solidHeader = TRUE, width = 8,
            div(class = "alert alert-warning",
                style = "margin-bottom:10px;padding:8px 12px;font-size:13px;",
              icon("exclamation-triangle"),
              " Lanza todas las etapas seleccionadas para ",
              strong("todas las muestras"), " simultáneamente. ",
              "Los jobs se encolan en BV-BRC y se monitorizan automáticamente."
            ),
            fluidRow(
              column(6,
                checkboxGroupInput("batch_tools",
                  "Etapas a lanzar (requieren FASTQ):",
                  choices = c(
                    "Calidad — FastqUtils"                  = "FastqUtils",
                    "Taxonómica — TaxonomicClassification"  = "TaxonomicClassification",
                    "Ensamblado — Assembly2"                = "Assembly2",
                    "Resistoma — MetagenomicReadMapping"    = "MetagenomicReadMapping"
                  ),
                  selected = c("FastqUtils", "TaxonomicClassification", "Assembly2")
                )
              ),
              column(6,
                checkboxInput("batch_cga",
                  tagList(icon("layer-group"),
                    " CGA — ComprehensiveGenomeAnalysis"),
                  value = FALSE),
                div(class = "alert alert-info",
                    style = "font-size:11px;padding:7px 10px;margin-bottom:8px;",
                  icon("info-circle"),
                  " CGA combina Ensamblado + Anotación + Resistoma en un job integral."
                ),
                uiOutput("batch_status_ui")
              )
            ),
            actionButton("btn_launch_all",
              tagList(icon("rocket"), " Lanzar análisis completo"),
              class = "btn-danger btn-lg",
              style = "width:100%;font-size:15px;padding:11px;margin-top:4px;"),

            hr(style = "margin:14px 0;"),
            div(class = "alert alert-info",
                style = "padding:8px 12px;font-size:13px;margin-bottom:8px;",
              icon("robot"), " ",
              strong("Pipeline automático (pasos 1→9): "),
              "ejecuta todas las etapas en el orden óptimo según van surgiendo ",
              "los resultados. El ensamblado alimenta la anotación; la anotación ",
              "alimenta MLST y filogenia. No necesitas intervenir entre pasos."
            ),
            uiOutput("auto_pipe_status_ui"),
            fluidRow(
              column(8,
                actionButton("btn_auto_pipeline",
                  tagList(icon("robot"), " Ejecutar pipeline automático (1→9)"),
                  class = "btn-primary btn-lg",
                  style = "width:100%;font-size:15px;padding:11px;margin-top:4px;")
              ),
              column(4,
                actionButton("btn_auto_pipeline_stop",
                  tagList(icon("stop"), " Detener"),
                  class = "btn-outline-secondary btn-lg",
                  style = "width:100%;font-size:15px;padding:11px;margin-top:4px;")
              )
            )
          ),
          box(title = tagList(icon("archive"), " Backup de resultados"),
              status = "success", solidHeader = TRUE, width = 4,
            div(style = "margin-bottom:10px;",
              tags$b("Exportar resultados actuales:"), br(), br(),
              downloadButton("dl_results_backup",
                tagList(icon("file-download"), " Descargar backup JSON"),
                class = "btn-success btn-sm", style = "width:100%;")
            ),
            hr(),
            div(
              tags$b("Importar backup guardado:"), br(), br(),
              fileInput("backup_file", NULL,
                accept   = ".json",
                placeholder = "Seleccionar .json...",
                buttonLabel = tagList(icon("folder-open"), " Cargar")
              ),
              uiOutput("backup_import_status")
            ),
            hr(),
            div(
              tags$b("Recuperar resultados de BV-BRC:"), br(),
              tags$small(style = "color:#888;",
                "Para jobs completados que no se importaron automáticamente."),
              br(), br(),
              actionButton("btn_recover_bvbrc_ids",
                tagList(icon("folder-open"), " Escanear workspace BV-BRC"),
                class = "btn-warning btn-sm", style = "width:100%;")
            )
          )
        ),
        fluidRow(
          column(12, style = "text-align:center;padding:20px 0 10px;",
            uiOutput("inicio_validacion"),
            br(),
            actionButton("btn_iniciar_analisis",
              label = tagList(icon("play-circle"), " Iniciar Análisis Paso a Paso"),
              class = "btn-success btn-lg",
              style = "padding:12px 40px;font-size:16px;border-radius:6px;")
          )
        )
      ),

      # ============================================================
      # PASO 1: CONTROL DE CALIDAD
      # ============================================================
      tabItem(tabName = "calidad",
        paso_header(1, "Control de Calidad", herramientas_info$calidad$descripcion),
        fluidRow(
          column(5,
            box(title = "Herramientas", status = "primary", solidHeader = TRUE, width = NULL,
              criterios_panel(herramientas_info$calidad$criterios),
              lapply(seq_len(nrow(herramientas_info$calidad$herramientas)), function(i) {
                h <- herramientas_info$calidad$herramientas[i,]
                tool_card(h$Herramienta, h$Descripcion, h$Link, h$Entrada)
              }),
              uiOutput("exec_panel_calidad")
            )
          ),
          column(7,
            box(title = "Registro de resultados", status = "success", solidHeader = TRUE, width = NULL,
              fluidRow(
                column(3, textInput("cq_muestra",     "Muestra",    placeholder = "ACINE1")),
                column(3, numericInput("cq_q30",      "% ≥ Q30",    value = NA, min=0, max=100, step=0.1)),
                column(3, numericInput("cq_gc",       "% GC",       value = NA, min=0, max=100, step=0.1)),
                column(3, numericInput("cq_reads",    "Reads (M)",  value = NA, min=0, step=0.1))
              ),
              fluidRow(
                column(3, numericInput("cq_dup",      "% Dup.",     value = NA, min=0, max=100, step=0.1)),
                column(3, numericInput("cq_longitud", "Long. media (bp)", value = NA, min=0, step=1)),
                column(6, br(), actionButton("add_calidad", "Agregar muestra", icon=icon("plus"),
                               class="btn-success btn-sm", style="margin-top:5px;"))
              ),
              br(),
              DTOutput("tabla_calidad"),
              textAreaInput("calidad_notas", "Observaciones / Acciones de trimming:", rows = 2,
                placeholder = "Ej: Trimming aplicado con Trimmomatic Q>30...")
            )
          )
        ),
        fluidRow(
          box(title = "Calidad por muestra — % Bases ≥ Q30", status = "info", solidHeader = TRUE, width = 8,
              plotlyOutput("plot_calidad", height = "320px")),
          box(title = "Contenido GC", status = "info", solidHeader = TRUE, width = 4,
              plotlyOutput("plot_gc", height = "320px"))
        ),
        nav_buttons("calidad")
      ),

      # ============================================================
      # PASO 2: ANÁLISIS TAXONÓMICO
      # ============================================================
      tabItem(tabName = "taxonomia",
        paso_header(2, "Análisis Taxonómico", herramientas_info$taxonomia$descripcion),
        fluidRow(
          column(5,
            box(title = "Herramientas", status = "primary", solidHeader = TRUE, width = NULL,
              criterios_panel(herramientas_info$taxonomia$criterios),
              lapply(seq_len(nrow(herramientas_info$taxonomia$herramientas)), function(i) {
                h <- herramientas_info$taxonomia$herramientas[i,]
                tool_card(h$Herramienta, h$Descripcion, h$Link, h$Entrada)
              }),
              uiOutput("exec_panel_taxonomia")
            )
          ),
          column(7,
            box(title = "Resultados taxonómicos", status = "success", solidHeader = TRUE, width = NULL,
              fluidRow(
                column(4, textInput("tax_muestra",   "Muestra",           placeholder = "ACINE1")),
                column(4, textInput("tax_especie",   "Especie identificada", placeholder = "A. baumannii")),
                column(4, numericInput("tax_conf",   "Confianza (%)",      value = NA, min=0, max=100, step=0.1))
              ),
              fluidRow(
                column(4, numericInput("tax_contam", "% Contaminación",  value = NA, min=0, max=100, step=0.1)),
                column(4, textInput("tax_org_contam","Organismo contam.", placeholder = "Homo sapiens")),
                column(4, br(), actionButton("add_taxonomia", "Agregar", icon=icon("plus"),
                               class="btn-success btn-sm", style="margin-top:5px;"))
              ),
              br(),
              DTOutput("tabla_taxonomia"),
              textAreaInput("tax_notas", "Interpretación taxonómica:", rows = 2)
            )
          )
        ),
        fluidRow(
          box(title = "Distribución de especies", status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_taxonomia_pie", height = "280px")),
          box(title = "% Contaminación por muestra", status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_contaminacion", height = "280px"))
        ),
        nav_buttons("taxonomia")
      ),

      # ============================================================
      # PASO 3: ENSAMBLADO
      # ============================================================
      tabItem(tabName = "ensamblado",
        paso_header(3, "Ensamblado Genómico", herramientas_info$ensamblado$descripcion),
        fluidRow(
          column(4,
            box(title = "Herramientas", status = "primary", solidHeader = TRUE, width = NULL,
              criterios_panel(herramientas_info$ensamblado$criterios),
              lapply(seq_len(nrow(herramientas_info$ensamblado$herramientas)), function(i) {
                h <- herramientas_info$ensamblado$herramientas[i,]
                tool_card(h$Herramienta, h$Descripcion, h$Link, h$Entrada)
              }),
              uiOutput("exec_panel_ensamblado")
            )
          ),
          column(8,
            box(title = "Métricas del ensamblado", status = "success", solidHeader = TRUE, width = NULL,
              fluidRow(
                column(3, textInput("ens_muestra",   "Muestra",        placeholder = "ACINE1")),
                column(3, numericInput("ens_cobert", "Cobertura (%)",  value = NA, min=0, max=100, step=0.1)),
                column(3, numericInput("ens_prof",   "Profundidad (X)",value = NA, min=0, step=0.1)),
                column(3, numericInput("ens_contigs","N° Contigs",     value = NA, min=0))
              ),
              fluidRow(
                column(3, numericInput("ens_n50",    "N50 (kb)",       value = NA, min=0, step=1)),
                column(3, numericInput("ens_tam",    "Tamaño (Mb)",    value = NA, min=0, step=0.01)),
                column(3, numericInput("ens_gc",     "GC (%)",         value = NA, min=0, max=100, step=0.1)),
                column(3, br(), actionButton("add_ensamblado", "Agregar", icon=icon("plus"),
                               class="btn-success btn-sm", style="margin-top:5px;"))
              ),
              br(),
              DTOutput("tabla_ensamblado"),
              textAreaInput("ens_notas", "Observaciones:", rows = 2,
                placeholder = "Ej: SPAdes v4.2 con k-mers 21,33,55,77...")
            )
          )
        ),
        fluidRow(
          box(title = "Cobertura y N50", status = "info", solidHeader = TRUE, width = 8,
              plotlyOutput("plot_ensamblado", height = "320px")),
          box(title = "Evaluación", status = "warning", solidHeader = TRUE, width = 4,
              uiOutput("eval_ensamblado"))
        ),
        nav_buttons("ensamblado")
      ),

      # ============================================================
      # PASO 4: ANOTACIÓN
      # ============================================================
      tabItem(tabName = "anotacion",
        paso_header(4, "Anotación Genómica", herramientas_info$anotacion$descripcion),
        fluidRow(
          column(4,
            box(title = "Herramientas", status = "primary", solidHeader = TRUE, width = NULL,
              criterios_panel(herramientas_info$anotacion$criterios),
              lapply(seq_len(nrow(herramientas_info$anotacion$herramientas)), function(i) {
                h <- herramientas_info$anotacion$herramientas[i,]
                tool_card(h$Herramienta, h$Descripcion, h$Link, h$Entrada)
              }),
              uiOutput("exec_panel_anotacion")
            )
          ),
          column(8,
            box(title = "Resultados de anotación", status = "success", solidHeader = TRUE, width = NULL,
              fluidRow(
                column(3, numericInput("anot_cds",       "Total CDS",         value = NA, min=0)),
                column(3, numericInput("anot_funcion",   "Con función (%)",    value = NA, min=0, max=100, step=0.1)),
                column(3, numericInput("anot_hipot",     "Hipotéticos (%)",    value = NA, min=0, max=100, step=0.1)),
                column(3, numericInput("anot_virulencia","Genes virulencia",   value = NA, min=0))
              ),
              fluidRow(
                column(3, numericInput("anot_rrna",  "ARNr",   value = NA, min=0)),
                column(3, numericInput("anot_trna",  "ARNt",   value = NA, min=0)),
                column(6, numericInput("anot_movil", "Elem. móviles", value = NA, min=0))
              ),
              h5("Categorías funcionales COG (Categoría: N genes, una por línea):"),
              textAreaInput("anot_cog", NULL, rows = 5,
                placeholder = "Metabolismo energético: 312\nTransporte: 245\nReplicación/Reparación: 189"),
              textAreaInput("anot_notas", "Observaciones:", rows = 2)
            )
          )
        ),
        fluidRow(
          box(title = "Distribución funcional COG", status = "info", solidHeader = TRUE, width = 8,
              plotlyOutput("plot_anotacion", height = "320px")),
          box(title = "Resumen", status = "success", solidHeader = TRUE, width = 4,
              uiOutput("resumen_anotacion"))
        ),
        nav_buttons("anotacion")
      ),

      # ============================================================
      # PASO 5: RESISTOMA
      # ============================================================
      tabItem(tabName = "resistoma",
        paso_header(5, "Resistoma", herramientas_info$resistoma$descripcion),
        fluidRow(
          column(4,
            box(title = "Herramientas", status = "primary", solidHeader = TRUE, width = NULL,
              criterios_panel(herramientas_info$resistoma$criterios),
              lapply(seq_len(nrow(herramientas_info$resistoma$herramientas)), function(i) {
                h <- herramientas_info$resistoma$herramientas[i,]
                tool_card(h$Herramienta, h$Descripcion, h$Link, h$Entrada)
              }),
              uiOutput("exec_panel_resistoma")
            )
          ),
          column(8,
            box(title = "Perfil de resistencia", status = "danger", solidHeader = TRUE, width = NULL,
              h5("Mapa S/I/R por grupo de antibióticos:"),
              DTOutput("tabla_resistencia"),
              br(),
              fluidRow(
                column(6,
                  textAreaInput("res_genes", "Genes AMR detectados (separados por coma):", rows = 3,
                    placeholder = "blaNDM-1, aac(6')-Ib, tetA, sul1, qnrS..."),
                  textAreaInput("res_mutaciones", "Mutaciones cromosómicas:", rows = 2,
                    placeholder = "gyrA S83L, parC E84V, rpoB D516V...")
                ),
                column(6,
                  h5("Clasificación MDR:"),
                  checkboxGroupInput("res_clasif", NULL,
                    choices = c("MDR — >3 grupos" = "MDR",
                                "XDR — >10 grupos" = "XDR",
                                "PDR — todos" = "PDR")),
                  uiOutput("clasif_badge"),
                  br(),
                  uiOutput("resfinder_result_display")
                )
              )
            )
          )
        ),
        fluidRow(
          box(title = "Mapa de calor — Resistencia antimicrobiana", status = "danger", solidHeader = TRUE, width = 12,
              plotlyOutput("plot_resistoma", height = "380px"))
        ),
        nav_buttons("resistoma")
      ),

      # ============================================================
      # PASO 6: ALINEAMIENTO
      # ============================================================
      tabItem(tabName = "alineamiento",
        paso_header(6, "Alineamiento de Secuencias", herramientas_info$alineamiento$descripcion),
        fluidRow(
          column(5,
            box(title = "Herramientas", status = "primary", solidHeader = TRUE, width = NULL,
              criterios_panel(herramientas_info$alineamiento$criterios),
              lapply(seq_len(nrow(herramientas_info$alineamiento$herramientas)), function(i) {
                h <- herramientas_info$alineamiento$herramientas[i,]
                tool_card(h$Herramienta, h$Descripcion, h$Link, h$Entrada)
              }),
              uiOutput("exec_panel_alineamiento")
            )
          ),
          column(7,
            box(title = "Resultados del alineamiento", status = "success", solidHeader = TRUE, width = NULL,
              fluidRow(
                column(5, textInput("ali_par",   "Comparación", placeholder = "ACINE1 vs ATCC17978")),
                column(4, selectInput("ali_tool", "Herramienta:",
                  choices = c("BLAST","MUSCLE","Clustal Omega","MAFFT","COBALT","T-Coffee"))),
                column(3, br(), actionButton("add_ali", "Agregar", icon=icon("plus"),
                               class="btn-success btn-sm", style="margin-top:5px;"))
              ),
              fluidRow(
                column(3, numericInput("ali_ident",  "% Identidad",  value = NA, min=0, max=100, step=0.01)),
                column(3, numericInput("ali_cobert", "% Cobertura",  value = NA, min=0, max=100, step=0.01)),
                column(3, numericInput("ali_gaps",   "% Gaps",       value = NA, min=0, max=100, step=0.01)),
                column(3, numericInput("ali_snps",   "SNPs",         value = NA, min=0))
              ),
              br(),
              DTOutput("tabla_alineamiento"),
              uiOutput("ebi_result_display"),
              textAreaInput("ali_notas", "Interpretación:", rows = 2)
            )
          )
        ),
        fluidRow(
          box(title = "Matriz de identidad", status = "info", solidHeader = TRUE, width = 12,
              plotlyOutput("plot_alineamiento", height = "340px"))
        ),
        nav_buttons("alineamiento")
      ),

      # ============================================================
      # PASO 7: MLST
      # ============================================================
      tabItem(tabName = "mlst",
        paso_header(7, "MLST — Tipificación Multilocus", herramientas_info$mlst$descripcion),
        fluidRow(
          column(4,
            box(title = "Herramientas", status = "primary", solidHeader = TRUE, width = NULL,
              criterios_panel(herramientas_info$mlst$criterios),
              lapply(seq_len(nrow(herramientas_info$mlst$herramientas)), function(i) {
                h <- herramientas_info$mlst$herramientas[i,]
                tool_card(h$Herramienta, h$Descripcion, h$Link, h$Entrada)
              }),
              uiOutput("exec_panel_mlst")
            )
          ),
          column(8,
            box(title = "Resultados MLST", status = "success", solidHeader = TRUE, width = NULL,
              fluidRow(
                column(5, textInput("mlst_esquema", "Esquema MLST:", placeholder = "A. baumannii Oxford")),
                column(4, selectInput("mlst_db", "Base de datos:", choices = c("PubMLST","Pasteur","Oxford","Otro"))),
                column(3)
              ),
              fluidRow(
                column(3, textInput("mlst_muestra", "Muestra",    placeholder = "ACINE1")),
                column(2, textInput("mlst_st",      "ST",         placeholder = "79")),
                column(5, textInput("mlst_alelos",  "Perfil alélico", placeholder = "1-1-1-1-1-1-1")),
                column(2, br(), actionButton("add_mlst", "Agregar", icon=icon("plus"),
                               class="btn-success btn-sm", style="margin-top:5px;"))
              ),
              br(),
              DTOutput("tabla_mlst"),
              textAreaInput("mlst_notas", "Interpretación epidemiológica:", rows = 3,
                placeholder = "Ej: ST79 — clon endémico en hospitales mexicanos...")
            )
          )
        ),
        fluidRow(
          box(title = "Distribución de Sequence Types", status = "info", solidHeader = TRUE, width = 12,
              plotlyOutput("plot_mlst", height = "280px"))
        ),
        nav_buttons("mlst")
      ),

      # ============================================================
      # PASO 8: FILOGENIA
      # ============================================================
      tabItem(tabName = "filogenia",
        paso_header(8, "Filogenia", herramientas_info$filogenia$descripcion),
        fluidRow(
          column(4,
            box(title = "Herramientas", status = "primary", solidHeader = TRUE, width = NULL,
              criterios_panel(herramientas_info$filogenia$criterios),
              lapply(seq_len(nrow(herramientas_info$filogenia$herramientas)), function(i) {
                h <- herramientas_info$filogenia$herramientas[i,]
                tool_card(h$Herramienta, h$Descripcion, h$Link, h$Entrada)
              }),
              uiOutput("exec_panel_filogenia")
            )
          ),
          column(8,
            box(title = "Resultados filogenéticos", status = "success", solidHeader = TRUE, width = NULL,
              fluidRow(
                column(4, selectInput("fil_metodo", "Método:",
                  choices = c("Máxima verosimilitud (IQ-TREE)","Neighbor-Joining","Bayesiano (MrBayes)"))),
                column(4, selectInput("fil_modelo", "Modelo:",
                  choices = c("GTR+G","GTR+G+I","HKY+G","TN93+G","Auto (ModelTest)"))),
                column(4, numericInput("fil_bootstrap", "Bootstrap mín. (%)", value = NA, min=0, max=100))
              ),
              textAreaInput("fil_newick", "Árbol Newick (opcional):", rows = 3,
                placeholder = "((ACINE1:0.0001,ACINE2:0.0001)99:0.0002,...);"),
              textAreaInput("fil_notas", "Interpretación filogenética:", rows = 3,
                placeholder = "Ej: ACINE1, 2 y 4 forman clado monofilético...")
            )
          )
        ),
        fluidRow(
          box(title = "Árbol filogenético (esquemático)", status = "info", solidHeader = TRUE, width = 8,
              plotlyOutput("plot_filogenia", height = "380px")),
          box(title = "Recursos visualización", status = "warning", solidHeader = TRUE, width = 4,
            p("Para árbol interactivo con anotaciones:"),
            lapply(seq_len(nrow(herramientas_info$filogenia$herramientas)), function(i) {
              h <- herramientas_info$filogenia$herramientas[i,]
              div(style = "margin-bottom:8px;",
                tags$a(href = h$Link, target = "_blank",
                       icon("external-link-alt"), " ", h$Herramienta,
                       style = "color:#1a5276;font-weight:bold;"), br(),
                tags$small(h$Descripcion, style = "color:#555;")
              )
            })
          )
        ),
        nav_buttons("filogenia")
      ),

      # ============================================================
      # PASO 9: GENÓMICA COMPARATIVA
      # ============================================================
      tabItem(tabName = "genomica",
        paso_header(9, "Genómica Comparativa", herramientas_info$genomica$descripcion),
        fluidRow(
          column(4,
            box(title = "Herramientas", status = "primary", solidHeader = TRUE, width = NULL,
              criterios_panel(herramientas_info$genomica$criterios),
              lapply(seq_len(nrow(herramientas_info$genomica$herramientas)), function(i) {
                h <- herramientas_info$genomica$herramientas[i,]
                tool_card(h$Herramienta, h$Descripcion, h$Link, h$Entrada)
              }),
              uiOutput("exec_panel_genomica")
            )
          ),
          column(8,
            box(title = "Resultados — ANI y pangenoma", status = "success", solidHeader = TRUE, width = NULL,
              h5("Matriz ANI:"),
              fluidRow(
                column(4, textInput("ani_m1",   "Muestra 1", placeholder = "ACINE1")),
                column(4, textInput("ani_m2",   "Muestra 2", placeholder = "ACINE2")),
                column(3, numericInput("ani_val","ANI (%)",   value = NA, min=90, max=100, step=0.01)),
                column(1, br(), actionButton("add_ani", "OK", icon=icon("plus"),
                               class="btn-success btn-sm", style="margin-top:5px;"))
              ),
              DTOutput("tabla_ani"),
              br(),
              h5("Pangenoma:"),
              fluidRow(
                column(3, numericInput("gc_core",  "Core",        value = NA, min=0)),
                column(3, numericInput("gc_acc",   "Accesorios",  value = NA, min=0)),
                column(3, numericInput("gc_uniq",  "Únicos",      value = NA, min=0)),
                column(3, numericInput("gc_total", "Total pan.",  value = NA, min=0))
              ),
              textAreaInput("gc_notas", "Interpretación:", rows = 2)
            )
          )
        ),
        fluidRow(
          box(title = "Heatmap ANI", status = "info", solidHeader = TRUE, width = 7,
              plotlyOutput("plot_ani", height = "380px")),
          box(title = "Pangenoma", status = "success", solidHeader = TRUE, width = 5,
              plotlyOutput("plot_pangenoma", height = "380px"))
        ),
        nav_buttons("genomica")
      ),

      # ============================================================
      # REPORTE FINAL
      # ============================================================
      tabItem(tabName = "reporte",
        div(style = "background:linear-gradient(135deg,#1a5276,#2e86c1);color:white;padding:14px 20px;border-radius:6px;margin-bottom:18px;",
          tags$h3(style = "color:white;margin:0 0 5px;", icon("file-alt"), " Reporte Final del Análisis Bioinformático"),
          tags$p(style = "color:#d6eaf8;margin:0;font-size:13px;",
            "Genera un documento con todos los resultados, gráficas y observaciones de cada etapa.")
        ),
        fluidRow(
          box(title = "Configuración", status = "primary", solidHeader = TRUE, width = 4,
            textInput("rep_titulo",  "Título",    value = "Reporte de Análisis Bioinformático — VEV"),
            textInput("rep_autor",   "Autor(es)", placeholder = "Nombre del analista"),
            dateInput("rep_fecha",   "Fecha",     value = Sys.Date()),
            textAreaInput("rep_resumen", "Resumen ejecutivo:", rows = 4,
              placeholder = "Breve descripción de resultados principales..."),
            selectInput("rep_formato", "Formato:",
              choices = c("HTML (recomendado)"="html_document",
                          "PDF (requiere LaTeX)"="pdf_document",
                          "Word (.docx)"="word_document")),
            hr(),
            h5("Secciones a incluir:"),
            checkboxGroupInput("rep_secciones", NULL,
              choices = c("Control de Calidad"="calidad","Análisis Taxonómico"="taxonomia",
                          "Ensamblado"="ensamblado","Anotación"="anotacion",
                          "Resistoma"="resistoma","Alineamiento"="alineamiento",
                          "MLST"="mlst","Filogenia"="filogenia",
                          "Genómica Comparativa"="genomica"),
              selected = PASOS_IDS),
            br(),
            actionButton("generar_reporte", "Generar Reporte", icon = icon("download"),
                         class = "btn-success btn-lg", style = "width:100%;"),
            br(), br(),
            uiOutput("reporte_descarga")
          ),
          box(title = "Vista previa", status = "info", solidHeader = TRUE, width = 8,
              uiOutput("preview_reporte"))
        ),
        # ══════════ REPORTE FINAL COMPLETO ══════════
        fluidRow(
          box(title = tagList(icon("file-medical-alt"), " Reporte final completo (HTML)"),
              status = "success", solidHeader = TRUE, width = 12,
            div(class = "alert alert-success", style = "padding:9px 13px;font-size:13px;",
              icon("check-circle"),
              " Documento con ", tags$b("las 9 etapas ordenadas"),
              ", tablas, gráficas e interpretación automática, generado a partir de los ",
              tags$b("resultados reales"), " importados de BV-BRC y Galaxy."),
            uiOutput("etapas_check_ui"),
            br(),
            fluidRow(
              column(8,
                textAreaInput("rep_notas_final", "Notas del analista (opcional):",
                              rows = 2, width = "100%",
                              placeholder = "Observaciones, contexto epidemiológico, decisiones..."),
                textAreaInput("rep_leyenda", "Leyenda al pie del reporte:",
                              rows = 2, width = "100%",
                              value = LEYENDA_REPORTE_DEFAULT),
                div(style = "font-size:11px;color:#888;margin-top:-8px;",
                    icon("info-circle"),
                    " Aparece al final del documento, bajo el bloque de validación. ",
                    "Se guarda con la sesión.")),
              column(4, br(),
                actionButton("btn_reporte_final",
                  tagList(icon("file-export"), " Generar reporte final"),
                  class = "btn-success btn-lg", style = "width:100%;"),
                br(), br(),
                uiOutput("reporte_final_dl"))
            )
          )
        )
      ),

      # ============================================================
      # LOG / CONSOLA
      # ============================================================
      tabItem(tabName = "log_panel",
        div(style = "background:linear-gradient(135deg,#2c3e50,#4a6274);color:white;padding:14px 20px;border-radius:6px;margin-bottom:18px;",
          tags$h3(style = "color:white;margin:0 0 4px;", icon("terminal"), " Consola de Procesos y Errores"),
          tags$p(style = "color:#b0c4d8;margin:0;font-size:13px;",
            "Registro en tiempo real de llamadas API, estados de jobs y errores del sistema.")
        ),
        fluidRow(
          column(12,
            div(style = "display:flex;gap:8px;margin-bottom:12px;flex-wrap:wrap;",
              actionButton("btn_sync_now",
                           tagList(icon("cloud-download-alt"), " Sincronizar resultados"),
                           class = "btn-success btn-sm"),
              actionButton("btn_diag_bvbrc", tagList(icon("stethoscope"), " Diagnóstico BV-BRC"),
                           class = "btn-primary btn-sm"),
              actionButton("btn_clear_log", tagList(icon("trash-alt"), " Limpiar"),
                           class = "btn-danger btn-sm"),
              downloadButton("dl_log_file",  tagList(icon("download"),  " Descargar .log"),
                             class = "btn-default btn-sm"),
              actionButton("btn_refresh_log", tagList(icon("sync-alt"), " Refrescar"),
                           class = "btn-default btn-sm"),
              div(style = "margin-left:auto;font-size:12px;color:#888;padding-top:6px;",
                  uiOutput("log_stats_ui", inline = TRUE))
            ),
            # Resultado del diagnóstico end-to-end
            uiOutput("diag_bvbrc_result"),
            # Leyenda de niveles
            div(style = "display:flex;gap:6px;margin-bottom:10px;flex-wrap:wrap;",
              lapply(list(
                list(n="API",   c="#8e44ad"), list(n="OK",    c="#27ae60"),
                list(n="INFO",  c="#2980b9"), list(n="WARN",  c="#e67e22"),
                list(n="ERROR", c="#e74c3c"), list(n="DEBUG", c="#95a5a6")
              ), function(x)
                tags$span(style=paste0("background:",x$c,";color:white;padding:1px 7px;border-radius:3px;font-size:11px;font-weight:bold;"), x$n)
              )
            ),
            # Consola de log
            div(style = paste0(
                  "background:#1e272e;color:#dfe6e9;font-family:'Courier New',monospace;",
                  "font-size:12px;border-radius:6px;padding:12px;",
                  "min-height:420px;max-height:620px;overflow-y:auto;"),
              uiOutput("log_console")
            )
          )
        )
      )

    )
  )
)

# ---- Credenciales (sobreescritas por credentials.R si existe) ----
CRED_BVBRC_USER <- ""
CRED_BVBRC_PASS <- ""
CRED_GALAXY_KEY <- ""
if (file.exists("credentials.R")) source("credentials.R", local = environment())

# ============================================================
# SERVER — cargado desde R/server_main.R
# ============================================================
source("R/server_main.R", local = TRUE)

shinyApp(ui = ui, server = server)
