# Publicar ANbio en GitHub

Guía para subir el proyecto **sin exponer credenciales ni datos de pacientes**.

---

## ⚠️ Antes de empezar — lee esto

Estos archivos de tu carpeta **contienen contraseñas o datos sensibles** y **no deben subirse**:

| Archivo | Qué contiene |
|---|---|
| `credentials.R` | Usuario y contraseña de BV-BRC, API key de Galaxy |
| `runapp.txt` | Contraseñas de BV-BRC y Galaxy en texto plano |
| `.RData` / `.Rhistory` | Pueden conservar credenciales de sesiones de R |
| `Datos/` | Secuencias FASTQ (datos de pacientes) |
| `*.anbio` | Sesiones de trabajo con resultados |
| `*.docx`, `*.pdf` | Documentos internos del proyecto |

El archivo **`.gitignore`** ya excluye todo lo anterior. Aun así, **verifica antes de publicar** (paso 3).

> Si alguna de esas contraseñas llegó a subirse alguna vez, **cámbiala**: borrarla del repositorio no la elimina del historial ni de las copias que otros hayan clonado.

---

## Paso 1 — Inicializar el repositorio

Desde la carpeta del proyecto:

```bash
git init
git add .gitignore
git commit -m "Añadir .gitignore"
```

> Añadir `.gitignore` **primero** evita que Git registre archivos sensibles en el primer commit.

---

## Paso 2 — Añadir los archivos del proyecto

```bash
git add .
git status
```

---

## Paso 3 — Verificar que no se cuele nada sensible

**Este paso es el importante.** Revisa la lista que muestra `git status`:

```bash
git status --short
```

Confirma que **NO** aparecen: `credentials.R`, `runapp.txt`, `Datos/`, `.RData`, `.Rhistory`, `*.anbio`.

Comprobación adicional — busca contraseñas en lo que vas a subir:

```bash
git grep -I -n -i -e "password" -e "contrase" -e "api_key" -e "apikey" -- :^.gitignore :^*.md
```

Si aparece algo real (no un ejemplo ni un comentario), quítalo antes de continuar:

```bash
git rm --cached ARCHIVO
```

---

## Paso 4 — Primer commit y publicación

```bash
git commit -m "ANbio: app Shiny para análisis bioinformático de patógenos NGS"
git branch -M main
git remote add origin https://github.com/USUARIO/ANbio.git
git push -u origin main
```

---

## Qué queda en el repositorio

```
ANbio/
├── app.R                     ← Interfaz y carga de módulos
├── credentials.example.R     ← Plantilla (sin claves reales)
├── reporte_final.Rmd         ← Reporte HTML de las 9 etapas
├── reporte_template.Rmd      ← Reporte clásico HTML/PDF/Word
├── instalar_dependencias.R
├── README.md
├── PUBLICAR_EN_GITHUB.md     ← Este archivo
├── .gitignore
└── R/
    ├── server_main.R         ├── api_bvbrc.R      ├── api_galaxy.R
    ├── api_ebi.R             ├── api_ncbi.R       ├── api_cge.R
    ├── parsers_results.R     ├── tree_plot.R      ├── session_state.R
    ├── ui_helpers.R          └── logger.R
```

---

## Instrucciones para quien clone el repositorio

Conviene que el README lo indique (ya lo hace en la sección 3):

```bash
git clone https://github.com/USUARIO/ANbio.git
cd ANbio
```

```r
source("instalar_dependencias.R")     # instala los paquetes
file.copy("credentials.example.R", "credentials.R")   # y rellenar
shiny::runApp("app.R")
```

---

## Sugerencias para el repositorio

| Elemento | Recomendación |
|---|---|
| **Visibilidad** | Privado mientras contenga referencias a datos de pacientes o al laboratorio |
| **Licencia** | MIT si quieres uso libre; sin licencia el código queda "todos los derechos reservados" |
| **Descripción** | "Aplicación Shiny para análisis bioinformático de patógenos NGS" |
| **Temas (topics)** | `bioinformatics`, `shiny`, `ngs`, `bvbrc`, `galaxy`, `amr`, `r` |

---

## Limpieza opcional antes de publicar

Archivos de trabajo que quizá no quieras en el repositorio:

| Archivo | Comentario |
|---|---|
| `config.txt`, `filogenetica R.txt`, `Snakefile secuenciación.txt` | Notas sueltas; revisa si aportan |
| `Analisis bioinformatico.docx` | Resultados del proyecto — excluido por `.gitignore` |
| `GUIA PRACTICA CURSO BIOINFORMATICA 2026.pdf` | Documento del curso — excluido por `.gitignore` |

Si quieres publicar alguno de los `.docx`/`.pdf`, quita esa regla del `.gitignore` y añádelo de forma explícita.
