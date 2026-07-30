# ============================================================
# Árbol filogenético — parser de Newick y dibujo con ggplot2
# Sin dependencias externas (no requiere ape/ggtree).
# ============================================================

# ---- Parser de Newick ----
# Devuelve data.frame: id | parent | label | brlen | soporte
# Soporta etiquetas entre comillas, longitudes de rama y valores de
# soporte (bootstrap) en los nodos internos: (A:0.1,B:0.2)100:0.3;
newick_parse <- function(txt) {
  if (is.null(txt) || !nchar(trimws(txt))) return(NULL)
  s <- gsub("[\r\n\t]", "", txt)
  # Proteger los espacios dentro de etiquetas entre comillas antes de
  # eliminar el resto del espaciado (si no, "A. baumannii X" se pega).
  ph <- "§ESP§"
  m <- gregexpr('"[^"]*"|\'[^\']*\'', s)
  regmatches(s, m) <- lapply(regmatches(s, m), function(q) gsub(" ", ph, q, fixed = TRUE))
  s <- gsub(" ", "", s, fixed = TRUE)
  s <- gsub(ph, " ", s, fixed = TRUE)
  s <- sub(";.*$", "", s)
  if (!nchar(s)) return(NULL)

  toks <- regmatches(s, gregexpr("\\(|\\)|,|[^(),]+", s))[[1]]
  parent <- integer(0); lab <- character(0); bl <- numeric(0)
  add <- function(p) {
    parent[[length(parent) + 1L]] <<- p
    lab[[length(lab) + 1L]]       <<- ""
    bl[[length(bl) + 1L]]         <<- NA_real_
    length(parent)
  }
  stack <- integer(0); target <- NA_integer_; after_close <- FALSE

  for (t in toks) {
    if (t == "(") {
      nd <- add(if (length(stack)) stack[length(stack)] else 0L)
      stack <- c(stack, nd); after_close <- FALSE
    } else if (t == ",") {
      after_close <- FALSE
    } else if (t == ")") {
      if (!length(stack)) next
      target <- stack[length(stack)]
      stack  <- stack[-length(stack)]
      after_close <- TRUE
    } else {
      nd <- if (after_close && !is.na(target)) target
            else add(if (length(stack)) stack[length(stack)] else 0L)
      partes <- strsplit(t, ":", fixed = TRUE)[[1]]
      nm <- if (length(partes)) partes[1] else ""
      b  <- if (length(partes) > 1) suppressWarnings(as.numeric(partes[2])) else NA_real_
      nm <- gsub('^"|"$', "", nm)
      nm <- gsub("^'|'$", "", nm)
      if (nchar(nm)) lab[nd] <- nm
      if (!is.na(b)) bl[nd] <- b
      after_close <- FALSE
    }
  }
  if (!length(parent)) return(NULL)
  data.frame(id = seq_along(parent), parent = parent,
             label = lab, brlen = bl, stringsAsFactors = FALSE)
}

# ---- Calcular coordenadas de un filograma rectangular ----
# x = distancia acumulada desde la raíz; y = orden de las puntas
newick_layout <- function(df) {
  if (is.null(df) || !nrow(df)) return(NULL)
  hijos <- split(df$id, df$parent)
  es_punta <- !(df$id %in% df$parent)
  raiz <- df$id[df$parent == 0L]
  if (!length(raiz)) raiz <- df$id[1]
  raiz <- raiz[1]

  # x acumulada (recorrido en anchura desde la raíz)
  x <- rep(NA_real_, nrow(df)); x[raiz] <- 0
  cola <- raiz
  while (length(cola)) {
    nd <- cola[1]; cola <- cola[-1]
    hs <- hijos[[as.character(nd)]]
    if (is.null(hs)) next
    for (h in hs) {
      b <- df$brlen[h]; if (is.na(b)) b <- 0
      x[h] <- x[nd] + b
      cola <- c(cola, h)
    }
  }

  # y: puntas en orden de recorrido en profundidad; internos = media de hijos
  orden <- integer(0)
  dfs <- function(nd) {
    hs <- hijos[[as.character(nd)]]
    if (is.null(hs)) { orden <<- c(orden, nd); return(invisible()) }
    for (h in hs) dfs(h)
  }
  dfs(raiz)
  y <- rep(NA_real_, nrow(df))
  y[orden] <- seq_along(orden)
  # asignar internos de las hojas hacia la raíz
  pendientes <- setdiff(df$id, orden)
  for (k in seq_len(length(pendientes) + 2L)) {
    for (nd in pendientes) {
      hs <- hijos[[as.character(nd)]]
      if (!is.null(hs) && all(!is.na(y[hs]))) y[nd] <- mean(range(y[hs]))
    }
    if (all(!is.na(y))) break
  }
  df$x <- x; df$y <- y; df$es_punta <- es_punta
  df
}

# ---- Dibujar el árbol (ggplot2) ----
# resaltar: vector de patrones (regex) que se colorean como "muestras propias"
plot_newick <- function(txt, titulo = "Árbol filogenético",
                        resaltar = NULL, escala = TRUE, renombrar = NULL) {
  df <- newick_layout(newick_parse(txt))
  if (is.null(df) || !nrow(df) || all(is.na(df$y))) return(NULL)

  # renombrar: vector con nombres = etiqueta original (o genome_id),
  # valores = nombre a mostrar. Útil porque BV-BRC etiqueta todos los
  # genomas propios con el mismo nombre científico.
  if (!is.null(renombrar) && length(renombrar)) {
    for (k in names(renombrar)) {
      hit <- which(df$label == k)
      if (length(hit)) df$label[hit] <- renombrar[[k]]
    }
  }

  raiz <- df$id[df$parent == 0L][1]
  segs <- df[df$id != raiz & !is.na(df$parent) & df$parent > 0, ]
  if (!nrow(segs)) return(NULL)
  segs$x0 <- df$x[segs$parent]; segs$y0 <- df$y[segs$parent]

  puntas <- df[df$es_punta, ]
  puntas$etiqueta <- ifelse(nchar(puntas$label), puntas$label, paste0("nodo_", puntas$id))
  puntas$propia <- FALSE
  if (!is.null(resaltar) && length(resaltar))
    puntas$propia <- vapply(puntas$etiqueta, function(e)
      any(vapply(resaltar, function(p) grepl(p, e, ignore.case = TRUE), logical(1))),
      logical(1))

  # Soporte (bootstrap) en nodos internos: etiquetas numéricas
  internos <- df[!df$es_punta & nchar(df$label) > 0, ]
  if (nrow(internos)) {
    sup <- suppressWarnings(as.numeric(internos$label))
    internos <- internos[!is.na(sup), ]
    internos$sup <- sup[!is.na(sup)]
  }

  rango_x <- max(df$x, na.rm = TRUE)
  margen  <- if (rango_x > 0) rango_x * 0.55 else 1

  p <- ggplot2::ggplot() +
    # ramas horizontales
    ggplot2::geom_segment(data = segs,
      ggplot2::aes(x = x0, xend = x, y = y, yend = y),
      linewidth = 0.7, color = "#2C3E50") +
    # conectores verticales
    ggplot2::geom_segment(data = segs,
      ggplot2::aes(x = x0, xend = x0, y = y0, yend = y),
      linewidth = 0.7, color = "#2C3E50") +
    ggplot2::geom_point(data = puntas,
      ggplot2::aes(x = x, y = y, color = propia), size = 2.6) +
    ggplot2::geom_text(data = puntas,
      ggplot2::aes(x = x, y = y, label = paste0("  ", etiqueta), color = propia),
      hjust = 0, size = 3.4, fontface = "bold") +
    ggplot2::scale_color_manual(values = c(`FALSE` = "#5D6D7E", `TRUE` = "#C0392B"),
                                guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-rango_x * 0.02, rango_x + margen)) +
    ggplot2::labs(title = titulo, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13))

  if (nrow(internos))
    p <- p + ggplot2::geom_text(data = internos,
      ggplot2::aes(x = x, y = y, label = round(sup)),
      hjust = 1.15, vjust = -0.45, size = 2.7, color = "#7F8C8D")

  # Barra de escala de distancia evolutiva
  if (escala && rango_x > 0) {
    paso <- signif(rango_x / 4, 1)
    y0 <- -0.6
    p <- p +
      ggplot2::annotate("segment", x = 0, xend = paso, y = y0, yend = y0,
                        linewidth = 0.8, color = "#34495E") +
      ggplot2::annotate("text", x = paso / 2, y = y0 - 0.45,
                        label = format(paso, scientific = FALSE),
                        size = 3, color = "#34495E") +
      ggplot2::expand_limits(y = y0 - 0.9)
  }
  p
}
