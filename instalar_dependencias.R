# Instalar todos los paquetes necesarios para la app bioinformática LESP
paquetes <- c(
  "shiny",
  "shinydashboard",
  "shinyjs",
  "DT",
  "ggplot2",
  "plotly",
  "dplyr",
  "tidyr",
  "rmarkdown",
  "knitr",
  "kableExtra"
)

cat("Verificando e instalando paquetes...\n\n")

for (pkg in paquetes) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("Instalando: ", pkg, "\n"))
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(paste0("[OK] ", pkg, "\n"))
  }
}

cat("\n¡Todos los paquetes están listos!\n")
cat("Para iniciar la app, ejecuta:\n")
cat("  shiny::runApp('app.R')\n")
cat("O bien: Ctrl+Enter en RStudio con app.R abierto.\n")
