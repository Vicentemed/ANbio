# LESP — Análisis Bioinformático de Patógenos NGS

**Laboratorio Estatal de Salud Pública — Aguascalientes**  
Aplicación Shiny para el análisis bioinformático de secuencias NGS (archivos FASTQ) de patógenos bacterianos.

---

## Tabla de contenido

1. [Descripción general](#1-descripción-general)
2. [Requisitos del sistema](#2-requisitos-del-sistema)
3. [Instalación y arranque](#3-instalación-y-arranque)
   · [**Cómo ejecutar un análisis (paso a paso)**](#3-bis-cómo-ejecutar-un-análisis-paso-a-paso)
4. [Estructura de archivos](#4-estructura-de-archivos)
5. [Inicio rápido — asistente de 3 pasos](#5-inicio-rápido--asistente-de-3-pasos)
   · [Resto de la pantalla de Inicio](#5-bis-resto-de-la-pantalla-de-inicio)
6. [Credenciales de plataformas API](#6-credenciales-de-plataformas-api)
7. [Pipeline automático (etapas 1→9)](#7-pipeline-automático-etapas-19)
8. [Sesión de trabajo reanudable (.anbio)](#8-sesión-de-trabajo-reanudable-anbio)
9. [Indicadores de progreso en la barra lateral](#9-indicadores-de-progreso-en-la-barra-lateral)
10. [Las 9 etapas de análisis](#10-las-9-etapas-de-análisis)
11. [Reporte final HTML](#11-reporte-final-html)
12. [Flujo de trabajo recomendado](#12-flujo-de-trabajo-recomendado)
13. [Diagnóstico y solución de problemas](#13-diagnóstico-y-solución-de-problemas)
14. [Notas técnicas de las APIs](#14-notas-técnicas-de-las-apis)
15. [Preguntas frecuentes](#15-preguntas-frecuentes)

---

## 1. Descripción general

La app guía al analista a través de las **9 etapas del flujo bioinformático** definido en la Guía Práctica del Curso de Bioinformática LESP 2026. Está diseñada para trabajar con archivos FASTQ provenientes de secuenciadores Illumina y organismos como *Acinetobacter baumannii*, *Mycobacterium tuberculosis* u otros patógenos prioritarios.

**Lo que hace la app:**

- **Asistente de 3 pasos**: seleccionar carpeta con los FASTQ → capturar credenciales → iniciar. Nada más.
- Detecta automáticamente las muestras y sus pares R1/R2 a partir de los nombres de archivo, y avisa de descargas incompletas (`.part`).
- **Pipeline automático dependency-aware**: lanza las 9 etapas en el orden óptimo, esperando a que cada dependencia termine (el ensamblado alimenta la anotación; la anotación alimenta MLST y filogenia).
- **Envía jobs directamente** a BV-BRC, Galaxy, EBI JDispatcher, NCBI BLAST y CGE/ResFinder vía API REST — sin abrir los sitios web.
- **Monitorea y descarga** los resultados en segundo plano cada 45 s, incluidos los archivos grandes almacenados en Shock.
- **No re-analiza lo ya hecho**: sincroniza con el workspace de BV-BRC antes de enviar y omite las muestras con resultados válidos (los jobs fallidos sí se pueden reintentar).
- **Sesión reanudable (`.anbio`)**: se autoguarda al cerrar la app; al reabrirla y cargar el archivo, consulta el avance real en los servidores y actualiza cada etapa.
- **Solo muestra resultados reales** — sin datos de ejemplo; las etapas sin resultados aparecen vacías con un aviso.
- Genera **gráficas interactivas** (plotly), **árbol filogenético** dibujado desde el Newick, y un **reporte final HTML** con las 9 etapas.
- **Diagnóstico integrado** que prueba la cadena completa de la API (token → carpeta → subida → listado → descarga) y una consola de log con niveles.

> Los jobs se ejecutan en los servidores de BV-BRC / Galaxy / EBI. **El análisis continúa aunque cierres la app**: los identificadores de tarea quedan guardados en la sesión y al reabrir se recupera el avance.

---

## 2. Requisitos del sistema

| Componente | Versión mínima |
|---|---|
| R | 4.4.0 o superior |
| RStudio | 2023.09 o superior (recomendado) |
| Pandoc | 1.12.3 o superior (incluido en RStudio) |
| Sistema operativo | Windows 10/11, macOS, Linux |
| Conexión a internet | Requerida para enviar jobs a plataformas externas |

### Paquetes R requeridos

```r
# Interfaz
shiny, shinydashboard, shinyjs, shinyFiles

# Tablas y gráficas
DT, ggplot2, plotly, dplyr

# Reporte
rmarkdown, knitr, kableExtra

# API / red (instalación adicional necesaria)
httr2, jsonlite, curl
```

Instala todos con:

```r
source("instalar_dependencias.R")
```

---

## 3. Instalación y arranque

### Paso 1 — Instalar dependencias

```r
source("instalar_dependencias.R")
```

Verifica e instala automáticamente todos los paquetes necesarios.

### Paso 2 — Configurar credenciales

Copia la plantilla y rellénala:

```bash
cp credentials.example.R credentials.R
```

```r
CRED_BVBRC_USER <- "tu_usuario"           # Usuario o correo BV-BRC
CRED_BVBRC_PASS <- "tu_contraseña"        # Contraseña BV-BRC
CRED_GALAXY_KEY <- "tu_api_key_galaxy"    # API key de usegalaxy.org
```

> ⚠️ **Seguridad:** `credentials.R` guarda contraseñas en **texto plano**. Está incluido en `.gitignore` — **nunca lo subas a un repositorio**.
>
> Si prefieres no guardarlas en disco, deja las cadenas vacías y captúralas desde la app (*Configurar credenciales*): no se escriben en ningún archivo.

**Dónde obtener las claves:**
- Cuenta BV-BRC gratuita: [www.bv-brc.org](https://www.bv-brc.org)
- API key Galaxy: usegalaxy.org → User → Preferences → Manage API Key

### Paso 3 — Iniciar la app

**Desde RStudio:** abre `app.R` y haz clic en **Run App**.

**Desde la consola de R:**
```r
setwd("C:/ruta/a/ANbio")
shiny::runApp("app.R")
```

La app se abrirá en `http://127.0.0.1:PUERTO`. Al iniciarse, realiza **auto-login** a BV-BRC y Galaxy con las credenciales guardadas.

---

## 3-bis. Cómo ejecutar un análisis (paso a paso)

### A. Análisis nuevo

1. **Prepara los datos.** Coloca los FASTQ de todas las cepas en una carpeta. Nombres típicos de Illumina:
   ```
   MiCarpeta/
   ├── CEPA1_S1_L001_R1_001.fastq.gz
   ├── CEPA1_S1_L001_R2_001.fastq.gz
   ├── CEPA2_S2_L001_R1_001.fastq.gz
   └── CEPA2_S2_L001_R2_001.fastq.gz
   ```
   > Verifica que **no haya archivos `.part`** (descargas incompletas). La app avisa si los encuentra.

2. **Abre la app** (`Run App` en RStudio).

3. **Paso 1 del asistente** — clic en *"Seleccionar carpeta..."* y elige la carpeta. Confirma que el recuadro muestre el número de muestras esperado.

4. **Paso 2** — clic en *"Configurar BV-BRC y Galaxy"*. Si ya están en `credentials.R`, aparecerán conectados (✓ verde) sin hacer nada.

5. **Paso 3** — escribe el organismo (p. ej. `Acinetobacter baumannii`) y pulsa **"Iniciar análisis completo"**.

6. **Espera.** La app sube los FASTQ y lanza las primeras etapas. A partir de ahí trabaja sola: cada 45 s revisa qué terminó y lanza lo siguiente.

7. **Sigue el avance** en la barra lateral (jobs activos) o en la pestaña *Log / Consola*.

8. **Al terminar** — menú *Generar Reporte* → **"Generar reporte final"** → descarga el HTML.

**Tiempos orientativos** (3 muestras, servidores con carga normal):

| Etapa | Duración aproximada |
|---|---|
| Subida de FASTQ | 1–3 min por muestra |
| Calidad | 10–20 s por muestra |
| Ensamblado | 5–20 min por muestra |
| Anotación | 5–15 min por muestra |
| Resistoma | < 1 min por muestra |
| Taxonomía (Kraken2) | 30–60 min |
| MLST | 2–5 min por muestra |
| Filogenia / Comparativa | 5–20 min |

> El total depende de la cola de BV-BRC y Galaxy. **No es necesario dejar la app abierta.**

### B. Cerrar y retomar después

1. **Cierra la app** cuando quieras — se guarda sola. Los jobs siguen corriendo en los servidores.
2. **Vuelve a abrirla** más tarde (horas o días).
3. Pestaña Inicio → caja **"Sesión de trabajo"** → elige la sesión en la lista → **"Cargar y reanudar"**.
4. La app consulta el avance real y actualiza las etapas:
   *"Sesión reanudada — 4 resultado(s) nuevos importados, 2 aún en proceso."*
5. Si el pipeline seguía activo, continúa lanzando las etapas pendientes.

> Para pasar el análisis a otro equipo: **"Descargar sesión"** → copia el `.anbio` → en el otro equipo **"Abrir archivo…"**.

### C. Repetir una sola etapa

Si una etapa falló o quieres re-ejecutarla:

1. Ve a la pestaña de esa etapa en el menú lateral.
2. Panel *"Ejecución directa vía API"* → selecciona muestra y herramienta → **Enviar**
   (o **"Enviar para todas las muestras"**).

> Si esa muestra ya tiene un resultado **válido**, la app avisa y reimporta en vez de re-analizar. Los resultados **fallidos o vacíos sí se reintentan**.

### D. Antes de empezar: comprobar la conexión

Si es la primera vez o algo no funciona, ejecuta el diagnóstico:

**Menú lateral → Log / Consola → botón "Diagnóstico BV-BRC"**

Prueba la cadena completa en 6 pasos y muestra dónde falla exactamente. Ver [Sección 13](#13-diagnóstico-y-solución-de-problemas).

---

## 4. Estructura de archivos

```
ANbio/
│
├── app.R                        ← UI + carga de módulos
├── credentials.R                ← Credenciales API (NO compartir)
├── reporte_final.Rmd            ← Reporte final HTML — 9 etapas, gráficas, árbol
├── reporte_template.Rmd         ← Plantilla del reporte clásico (HTML/PDF/Word)
├── instalar_dependencias.R      ← Script de instalación de paquetes
├── README.md                    ← Este documento
│
├── R/
│   ├── server_main.R            ← Lógica principal del servidor Shiny
│   ├── api_bvbrc.R              ← Cliente REST BV-BRC (auth, upload, submit, poll, Shock)
│   ├── api_galaxy.R             ← Cliente REST Galaxy (upload, run, wait, download)
│   ├── api_ebi.R                ← Cliente EBI JDispatcher (MUSCLE, Clustal, MAFFT)
│   ├── api_ncbi.R               ← Cliente NCBI BLAST
│   ├── api_cge.R                ← Cliente CGE/ResFinder (DTU)
│   ├── parsers_results.R        ← Parsers de resultados por herramienta
│   ├── tree_plot.R              ← Parser de Newick + dibujo del árbol (sin ape/ggtree)
│   ├── session_state.R          ← Sesión .anbio: guardar, cargar, autoguardado
│   ├── ui_helpers.R             ← Componentes de interfaz reutilizables
│   └── logger.R                 ← Sistema de log de la app
│
├── Datos/                       ← Carpeta con archivos FASTQ de las muestras
│   ├── ACINE1_S2_L001_R1_001.fastq.gz
│   ├── ACINE1_S2_L001_R2_001.fastq.gz
│   └── ...
│
└── ~/.anbio_session/            ← (fuera del proyecto) estado y sesiones
    ├── anbio_state.rds          ← Autoguardado interno
    └── sesiones/
        └── sesion_activa.anbio  ← Sesiones reanudables
```

> Los archivos FASTQ pueden tener extensiones `.fastq`, `.fastq.gz` o `.fq.gz`.  
> La app extrae los nombres de muestra eliminando automáticamente sufijos de Illumina (`_S#_L00#_R#_001`, `_R1`, `_R2`, etc.).

---

## 5. Inicio rápido — asistente de 3 pasos

Al abrir la app, lo primero que aparece es el asistente. **Con estos tres pasos se lanza el análisis completo.**

```
┌─ 1 ─────────────────┐  ┌─ 2 ─────────────────┐  ┌─ 3 ─────────────────┐
│ Carpeta con datos   │  │ Credenciales        │  │ Ejecutar 9 etapas   │
│ crudos              │  │                     │  │                     │
│ [Seleccionar        │→ │ [Configurar BV-BRC  │→ │ Organismo: [______] │
│  carpeta...]        │  │  y Galaxy]          │  │ [Iniciar análisis]  │
│ ✓ 3 muestra(s)      │  │ ✓ BV-BRC ✓ Galaxy   │  │                     │
└─────────────────────┘  └─────────────────────┘  └─────────────────────┘
```

### Paso 1 — Carpeta con datos crudos

Botón **"Seleccionar carpeta..."**. La app:

- Detecta todos los `.fastq` / `.fastq.gz` / `.fq.gz` de la carpeta.
- Empareja automáticamente R1 y R2 de cada muestra.
- Deriva el nombre de muestra quitando sufijos de Illumina (`_S#`, `_L00#`, `_R1/_R2`, `_001`).
- **Avisa si hay descargas incompletas** (`.part`), que romperían el análisis.

Puede seleccionarse **una o varias cepas**: todas las muestras encontradas entran al pipeline.

### Paso 2 — Credenciales

Un solo botón abre el modal con las dos plataformas. El indicador muestra ✓ o ✗ para cada una:

| Plataforma | Necesaria para |
|---|---|
| **BV-BRC** | Calidad, ensamblado, anotación, resistoma, filogenia, comparativa |
| **Galaxy** | Taxonomía (Kraken2) y MLST |

> Si Galaxy no está configurado, el pipeline avisa y continúa con el resto de las etapas.

### Paso 3 — Ejecutar

Se escribe el organismo (p. ej. *Acinetobacter baumannii*) y se pulsa **"Iniciar análisis completo"**. A partir de ahí el pipeline trabaja solo — ver [Sección 7](#7-pipeline-automático-etapas-19).

---

## 5-bis. Resto de la pantalla de Inicio

### Información del proyecto

| Campo | Ejemplo |
|---|---|
| Proyecto | Análisis ACINE 2025 |
| Organismo | *Acinetobacter baumannii* |
| Analista | María García López |
| Institución | LESP Aguascalientes |
| Fecha | 27/05/2026 |
| Tipo de muestra | Bacterias no tuberculosas |

Esta información se incluye automáticamente en el reporte final.

### Sesión de trabajo

Caja para guardar, cargar y descargar la sesión reanudable. Ver [Sección 8](#8-sesión-de-trabajo-reanudable-anbio).

### Estado de credenciales API

Muestra si BV-BRC y Galaxy están conectados (verde) o no (rojo). El botón **Configurar credenciales** abre el modal para introducir o actualizar las claves sin reiniciar la app.

### Parámetros de la corrida NGS (BaseSpace)

> Estos valores **NO se calculan a partir de los FASTQ**. Son métricas reportadas por **Illumina BaseSpace** tras la corrida. El operador los consulta allí y los registra aquí.

| Parámetro | Criterio aceptable |
|---|---|
| Clusters PF (%) | ≥ 80% |
| % ≥ Q30 | ≥ 75–85% |
| Error Rate (%) | < 1% |
| % PhiX | 1–5% |

### Análisis por lotes

Caja roja que permite lanzar etapas concretas para todas las muestras. Útil cuando solo se quiere repetir una etapa; para el flujo normal se usa el **pipeline automático** del asistente.

### Backup de resultados

Caja verde con dos funciones:
- **Descargar backup JSON** — exporta todas las tablas de resultados actuales.
- **Importar backup** — restaura resultados desde un JSON guardado previamente.

> El backup JSON archiva **resultados**; la sesión `.anbio` archiva además los **jobs en curso** para poder reanudar.

---

## 6. Credenciales de plataformas API

La app se conecta a cinco plataformas externas:

| Plataforma | Tipo de acceso | Cómo obtener |
|---|---|---|
| **BV-BRC** | Usuario + contraseña | Cuenta gratuita en [bv-brc.org](https://www.bv-brc.org) |
| **Galaxy** | API key | usegalaxy.org → User → Preferences → Manage API Key |
| **EBI JDispatcher** | Pública (sin cuenta) | Sin configuración |
| **NCBI BLAST** | Pública (sin cuenta) | Sin configuración |
| **CGE/ResFinder** | Pública (sin cuenta) | Sin configuración |

### Configuración en `credentials.R`

```r
CRED_BVBRC_USER  <- "usuario@correo.com"
CRED_BVBRC_PASS  <- "contraseña"
CRED_GALAXY_KEY  <- "abc123...api_key"
```

Las credenciales se cargan al arrancar la app y se usan automáticamente. Si hay un error de autenticación, el modal de credenciales (icono de llave en la cabecera) permite actualizar las claves sin reiniciar.

### Auto-login al iniciar

Al arrancar, la app:
1. Intenta autenticarse en BV-BRC con `CRED_BVBRC_USER` / `CRED_BVBRC_PASS`.
2. Verifica la API key de Galaxy y crea (o recupera) el historial `ANbio_LESP`.
3. Muestra el estado en el panel de credenciales de la pestaña Inicio.

---

## 7. Pipeline automático (etapas 1→9)

Es el modo de trabajo principal. Se lanza desde el **paso 3 del asistente** (o con el botón *"Ejecutar pipeline automático"* del panel de lotes).

### 7.1 Cómo decide el orden

El motor revisa cada 45 s qué terminó y lanza lo que ya puede correr. No es una secuencia rígida: **cada etapa arranca en cuanto su dependencia está lista.**

| Ola | Etapas | Se lanza cuando… |
|---|---|---|
| **A** | 1 Calidad · 2 Taxonomía · 3 Ensamblado · 5 Resistoma | de inmediato (solo necesitan los FASTQ) |
| **B** | 4 Anotación | terminó el ensamblado de esa muestra (usa sus *contigs*) |
| **C** | 7 MLST | hay *contigs* disponibles (se envían a Galaxy) |
| **D** | 8 Filogenia · 6 Alineamiento · 9 Comparativa | terminaron **todas** las anotaciones (≥ 2 genomas) |

> El **paso 6 (Alineamiento)** queda cubierto por los *hits* bidireccionales de la genómica comparativa y por el alineamiento múltiple de genes core de la filogenia.

### 7.2 Reparto entre plataformas

| Etapa | Plataforma | Herramienta |
|---|---|---|
| 1 Calidad | BV-BRC | `FastqUtils` (FastQC) |
| 2 Taxonomía | **Galaxy** | `Kraken2` |
| 3 Ensamblado | BV-BRC | `GenomeAssembly2` (Unicycler) |
| 4 Anotación | BV-BRC | `GenomeAnnotation` (RASTtk) |
| 5 Resistoma | BV-BRC | `MetagenomicReadMapping` (CARD) |
| 7 MLST | **Galaxy** | `mlst` (esquema Pasteur) |
| 8 Filogenia | BV-BRC | `CodonTree` (RAxML) |
| 9 Comparativa | BV-BRC | `GenomeComparison` |

> La taxonomía y el MLST **no** usan BV-BRC: su `TaxonomicClassification` falla de forma sistemática y `CoreGenomeMLST` no existe como app. Ver [Sección 14](#14-notas-técnicas-de-las-apis).

### 7.3 Sin trabajo duplicado

Antes de enviar cualquier job, la app **sincroniza con el workspace de BV-BRC**:

- Registra como completadas las carpetas de resultados que ya existen.
- **Omite** esas muestras/etapas — no vuelve a subir los FASTQ ni a re-analizar.
- Los jobs **fallidos o vacíos se ignoran** en esa comprobación, de modo que sí pueden reintentarse.

Esto funciona incluso en una sesión nueva, sin historial local.

### 7.4 Seguimiento

- Panel lateral con los jobs en curso y su tiempo transcurrido.
- Barra de progreso del pipeline (`n/N etapas por muestra completadas`).
- Al completar todo, notificación de **pipeline finalizado**.
- Botón **Detener**: corta el envío de nuevas etapas; los jobs ya enviados siguen en los servidores.

### 7.5 Enviar para todas las muestras (por etapa)

En el panel de ejecución de **cada etapa** aparece un segundo botón azul:

```
[ ▶ Enviar ]                                  ← una muestra
[ ⊞ Enviar para todas las muestras (N) ]      ← todas las muestras
```

El número `(N)` se actualiza automáticamente. Al hacer clic, la app envía el job para cada muestra detectada mostrando una barra de progreso.

### 7.3 Comportamiento por plataforma

| Plataforma | Comportamiento |
|---|---|
| **BV-BRC** | Upload R1+R2 por muestra → submit job → sin espera de bloqueo |
| **Galaxy** | Upload → espera `state=ok` → submit. Más lento por la espera por muestra |
| **EBI / NCBI / CGE** | Submit directo (sin upload de archivos) — muy rápido |

### 7.4 Monitoreo automático

Los jobs enviados quedan registrados y el sistema de polling (cada 45 s) los monitorea automáticamente. Al completarse:
- Los resultados se descargan del workspace.
- Se parsean y se importan directamente a las tablas y gráficas de cada etapa.
- El badge en la barra lateral cambia de amarillo (activo) a verde (completado).

---

## 8. Sesión de trabajo reanudable (.anbio)

La app tiene **tres mecanismos** de persistencia, complementarios:

| Mecanismo | Para qué sirve |
|---|---|
| **Autoguardado interno (RDS)** | Red de seguridad ante cierres inesperados |
| **Sesión `.anbio`** | **Cerrar y retomar** el análisis; mover entre equipos |
| **Backup JSON** | Archivar resultados finales o compartirlos |

### 8.1 Guardado automático (RDS)

**La app guarda el estado automáticamente.** No es necesario hacer nada.

**Archivo:** `C:/Users/<usuario>/.anbio_session/anbio_state.rds`

Se guarda en tres momentos:
- Cada **45 segundos** durante el ciclo de monitoreo
- Cada vez que se **envía un job** (inmediatamente)
- Cada vez que se **importan resultados** de un job completado

**Qué incluye el guardado:**

| Dato | Guardado |
|---|---|
| Lista de jobs (task_id, job_id, estado) | ✅ |
| Tablas de resultados (calidad, taxonomía, ensamblado, MLST, ANI…) | ✅ |
| Muestras detectadas y carpeta FASTQ | ✅ |
| Etapas marcadas como completadas | ✅ |
| Usuario BV-BRC | ✅ |
| Contraseñas o tokens | ❌ (se re-autentican al abrir) |

**Expiración:** el archivo se descarta automáticamente si tiene más de **7 días** (coincide con la caducidad de los jobs en EBI y NCBI).

### 8.2 Restaurar sesión al reabrir

Al iniciar la app, si existe un estado guardado reciente, aparece automáticamente un **modal de restauración**:

```
┌─────────────────────────────────────────────────┐
│  Sesión guardada encontrada                      │
│  Guardada el 27/05/2026 a las 14:32 (3 min atrás)│
│                                                   │
│  Jobs registrados:                                │
│  • 12 jobs en total                              │
│  • 4 activos — retomarán monitoreo automático    │
│  • 8 completados                                  │
│                                                   │
│  [ Restaurar sesión ]   [ Nueva sesión ]          │
└─────────────────────────────────────────────────┘
```

Al hacer clic en **Restaurar sesión**:
- Se recuperan todas las tablas de resultados.
- Los jobs que seguían corriendo en BV-BRC/Galaxy **reanudan el monitoreo** automáticamente — siguen vivos en el servidor externo aunque la app haya estado cerrada.

### 8.2-bis Archivo de sesión `.anbio` (reanudable)

Además del guardado interno, la app maneja **archivos de sesión portables** que permiten **cerrar la app y retomar el análisis después** — incluso días más tarde o en otro equipo.

**Dónde está:** pestaña Inicio → caja **"Sesión de trabajo (reanudar análisis)"**.

| Acción | Qué hace |
|---|---|
| **Guardar ahora** | Escribe `~/.anbio_session/sesiones/<nombre>.anbio` |
| **Cargar y reanudar** | Lista las sesiones guardadas; al cargar, **consulta el avance real** |
| **Descargar sesión** | Exporta el `.anbio` para archivarlo o moverlo de equipo |
| **Abrir archivo…** | Importa un `.anbio` externo |

**Autoguardado:** al cerrar la app se guarda sola, sin intervención. También se guarda al arrancar el pipeline.

#### Qué ocurre al reanudar

Este es el punto clave: **al cargar la sesión, la app no muestra una foto vieja** — va a los servidores y actualiza:

```
Cargar sesión
      │
      ├─ 1. Restaura muestras, jobs, resultados y metadatos del proyecto
      ├─ 2. Sincroniza con el workspace de BV-BRC (descubre resultados nuevos)
      ├─ 3. Consulta el estado de cada job pendiente en BV-BRC y Galaxy
      └─ 4. Descarga e importa lo que se completó mientras la app estaba cerrada
             ↓
      "Sesión reanudada — 4 resultado(s) nuevos importados, 2 aún en proceso."
```

Cada etapa queda actualizada al avance real, **independientemente de cuánto tiempo estuvo cerrada la app**.

#### Contenido del `.anbio`

| Dato | Incluido |
|---|---|
| Identificadores de tarea (task_id) de BV-BRC y Galaxy | ✅ |
| Todas las tablas de resultados (9 etapas) | ✅ |
| Antibiograma predicho, pangenoma, árbol Newick | ✅ |
| Estado del pipeline automático (qué se envió, contigs, genome_ids) | ✅ |
| Muestras, carpeta FASTQ y metadatos del proyecto | ✅ |
| Contraseñas o tokens | ❌ (se re-autentica al abrir) |

> Los FASTQ originales no se incluyen. Si cambias de equipo, vuelve a seleccionar la carpeta para poder lanzar etapas nuevas; los resultados ya obtenidos se restauran completos.

### 8.3 Backup JSON manual

Para guardar resultados de forma permanente (más allá de los 7 días del RDS, o para archivar un análisis terminado):

**Exportar:**  
Pestaña Inicio → caja "Backup de resultados" → **Descargar backup JSON**

El archivo generado tiene nombre `anbio_backup_YYYYMMDD_HHMMSS.json` e incluye:

```json
{
  "version": "1.0",
  "timestamp": "2026-05-27T14:32:00",
  "proyecto": { "nombre": "...", "organismo": "...", "analista": "..." },
  "muestras": ["ACINE1", "ACINE2", "ACINE3", "ACINE4"],
  "calidad":      [ ... ],
  "taxonomia":    [ ... ],
  "ensamblado":   [ ... ],
  "mlst":         [ ... ],
  "alineamiento": [ ... ],
  "ani":          [ ... ],
  "resfinder_df": [ ... ],
  "jobs_summary": [ ... ]
}
```

**Importar:**  
Pestaña Inicio → caja "Backup de resultados" → **Cargar** → seleccionar `.json`

Al importar se restauran todas las tablas. Las gráficas se actualizan automáticamente.

### 8.4 Resumen de escenarios

| Situación | Solución |
|---|---|
| Cerré la app accidentalmente con jobs corriendo | Reabrir → modal automático → "Restaurar sesión" |
| Corte de luz / apagado con jobs en cola | Reabrir → restaurar → los jobs siguen corriendo en BV-BRC |
| Quiero archivar resultados de un análisis terminado | "Descargar backup JSON" → guardar junto con los FASTQ |
| Pasaron más de 7 días y el RDS expiró | Cargar el JSON que se descargó antes |
| Comparto el análisis con otro analista | Enviar el JSON; el destinatario lo importa en su app |

> **Nota:** los archivos FASTQ originales no se guardan en el RDS ni en el JSON. Si se cambia de PC o de carpeta, hay que volver a indicar la ruta FASTQ en la pestaña Inicio para poder lanzar nuevos jobs. Los resultados ya obtenidos se restauran completamente.

---

## 9. Indicadores de progreso en la barra lateral

### 9.1 Badge por etapa

Cada ítem del menú lateral refleja el estado de los jobs de esa etapa:

| Badge | Significado |
|---|---|
| 🔄 Spinner amarillo `activo` | Hay jobs corriendo en esa etapa |
| 🟢 Verde `N ok` | Todos completados sin errores |
| 🟠 Naranja `N ok` | Algunos completados, algunos con error |
| 🔴 Rojo `error` | Jobs con error, ninguno completado |
| Sin badge | Sin jobs registrados aún |

### 9.2 Panel de jobs en curso

Debajo de las etapas en la barra lateral aparece un mini-panel con los jobs activos y los completados en los últimos 5 minutos:

```
EN PROCESO (2/3)
🔄 Calidad   FastqUtils     3m
🔄 Taxon.    TaxonomicCl.   2m
✅ Ensam.    Assembly2      1m
```

Muestra etapa, herramienta y tiempo transcurrido desde el envío.

### 9.3 Barra de progreso global

En la parte superior de cada etapa: stepper de 9 círculos (gris → azul activo → verde completado) y barra de porcentaje. Los círculos son clicables para saltar directamente a cualquier etapa.

### 9.4 Log / Consola

La pestaña **"Log / Consola"** en el menú lateral registra cada acción de la app con marca de tiempo y nivel:

| Nivel | Color | Uso |
|---|---|---|
| `API` | Morado | Llamadas a servicios externos |
| `OK` | Verde | Operación exitosa |
| `INFO` | Azul | Eventos informativos |
| `WARN` | Amarillo | Advertencias no críticas |
| `ERROR` | Rojo | Errores que requieren atención |
| `DEBUG` | Gris | Detalle técnico |

El log se puede descargar como `.txt` para diagnóstico.

---

## 10. Las 9 etapas de análisis

Cada etapa tiene la misma estructura:

1. **Encabezado** — número, nombre y descripción (fondo azul degradado).
2. **Panel de herramientas** — criterios de calidad + tarjetas con links a herramientas externas.
3. **Panel de ejecución directa vía API** — selecciona muestra + herramienta y envía el job con un clic. Botón adicional **"Enviar para todas las muestras (N)"**.
4. **Panel de resultados** — tabla editable para ingresar datos manualmente o importados automáticamente.
5. **Gráfica interactiva** — plotly, se actualiza al llegar resultados.
6. **Panel de evaluación** — cuántas muestras cumplen los criterios.
7. **Notas** — campo libre para observaciones.
8. **Barra de navegación** — botones Anterior / Siguiente y checkbox "Etapa completada".

---

### Paso 1 — Control de Calidad

**Herramienta API:** `FastqUtils` vía BV-BRC  
**Plataforma:** [bvbrc.org/app/FastqUtil](https://www.bvbrc.org/app/FastqUtil)

**Criterios de aprobación:**

| Parámetro | Mínimo |
|---|---|
| % Q30 | ≥ 75–85% |
| Clusters PF | ≥ 80% |
| Error Rate | < 1% |
| % PhiX | 1–5% |

**Resultados importados automáticamente:** % Q30, % GC, Reads (M), % Duplicados, longitud media por muestra.

---

### Paso 2 — Análisis Taxonómico

**Herramienta principal:** `Kraken2` vía **Galaxy** (base RefSeq Standard-Full)

> **Por qué Galaxy:** el servicio `TaxonomicClassification` de BV-BRC falla de forma sistemática (`wrapper command failed 256`) con todas las combinaciones de parámetros probadas. Kraken2 en Galaxy es la ruta que sí funciona.

**Criterio:** Contaminación < 5%.

**Resultados importados:** género dominante y su %, especie principal, % no clasificado, top 5 de especies.

> **Lectura del porcentaje:** la app reporta el **% a nivel de género**, no de especie. Kraken2 deja la mayoría de las lecturas en el género cuando la especie no es discriminable (p. ej. complejo *A. calcoaceticus/baumannii*), por lo que el % de especie subestimaría la identificación. Un 95 % de género con la especie dominante muy por encima de la siguiente es una identificación sólida.

**Tiempo:** puede tardar **30–60 min** en el Galaxy público por el tamaño de la base de datos.

---

### Paso 3 — Ensamblado

**Herramientas API:**
- `Assembly2` vía BV-BRC (principal para lanzamiento por lotes)
- `SPAdes` vía Galaxy
- `Unicycler` vía Galaxy

**Criterios:**

| Parámetro | Criterio |
|---|---|
| Cobertura | ≥ 95% |
| Profundidad | ≥ 30X |
| Contigs | < 200 |
| N50 | > 50 kb |

**Resultados importados:** N° contigs, N50 (kb), tamaño total (Mb), GC (%), **profundidad** y **nº de bases ambiguas (Ns)**.

> BV-BRC no genera un reporte QUAST; la app **calcula las métricas directamente del `contigs.fasta`**, aprovechando que Unicycler anota longitud y cobertura en cada encabezado.

---

### Paso 4 — Anotación

**Herramienta:** `GenomeAnnotation` (RASTtk) vía BV-BRC
**Alternativa:** `ComprehensiveGenomeAnalysis` (CGA) — cubre Ensamblado + Anotación + Resistoma en un job

**Criterio:** Cobertura de anotación > 80%.

**Resultados importados:** total CDS, ARNr, ARNt, genes hipotéticos, % con función asignada, y **completitud / contaminación (CheckM)**.

> La anotación también produce el **antibiograma predicho** (`amr-sir.txt`), que la app importa para el paso 5.

---

### Paso 5 — Resistoma

**Herramientas:**
- `MetagenomicReadMapping` (CARD) vía BV-BRC — cribado desde lecturas (formato KMA)
- **RGI** dentro de la anotación — asignación de alelos sobre el genoma ensamblado
- `ResFinder` vía CGE/DTU (API pública)

**Clasificación:**

| Categoría | Criterio |
|---|---|
| MDR | Resistencia a > 3 grupos |
| XDR | Resistencia a > 10 grupos |
| PDR | Todos los grupos |

**Resultados importados:**

1. **Genes AMR** con identidad, cobertura, profundidad y **clase de antibiótico inferida** del nombre del gen.
2. **Antibiograma predicho (S/I/R)** por antibiótico, con su confianza (F1), generado por la anotación de BV-BRC.

> **Qué método usar en el reporte clínico:** RGI (basado en la proteína ensamblada) da la asignación de alelo más precisa — p. ej. distingue `tet(B)` de `tet(A)` o `ADC-73` de `ADC-25`. El mapeo de lecturas (KMA) sirve como cribado sensible. Ambos se importan.

---

### Paso 6 — Alineamiento

**Herramientas API:**
- `MUSCLE`, `Clustal Omega`, `MAFFT` vía EBI JDispatcher (públicas, sin cuenta)
- `BLAST` vía NCBI (público)
- `Homology` vía BV-BRC

**Criterios:** Identidad > 95%, E-value < 1e-5, cobertura ≥ 90%.

---

### Paso 7 — MLST

**Herramienta:** `mlst` vía **Galaxy** (esquemas PubMLST/Pasteur)

> **Por qué Galaxy:** `CoreGenomeMLST` **no existe** como app en BV-BRC (verificado con `enumerate_apps`). El pipeline descarga los *contigs* del ensamblado de BV-BRC, los sube a Galaxy y ejecuta `mlst`.

**Criterio:** ST único con todos los loci cubiertos al ≥ 95%.

**Resultados importados:** esquema detectado, ST y perfil alélico por muestra.

> El esquema también confirma la especie (p. ej. `abaumannii_2` → *A. baumannii*), sirviendo como verificación independiente de la taxonomía.

---

### Paso 8 — Filogenia

**Herramienta:** `CodonTree` vía BV-BRC (RAxML sobre proteínas concatenadas de genes core)

> El app `PhylogeneticTree` exige un grupo externo (`in_genome_ids` + `out_genome_ids`); la app usa `CodonTree`, que trabaja con una sola lista de genomas.

**Resultados importados:** árbol Newick + estadísticas del alineamiento (nº de genes, aminoácidos alineados).

**Criterio:** Bootstrap ≥ 70%.

**Visualización:** la app **dibuja el árbol** (pestaña Filogenia y reporte final):

- Filograma con longitudes de rama proporcionales.
- Valores de bootstrap sobre los nodos.
- Barra de escala de distancia evolutiva.
- **Muestras propias en rojo**, referencias en gris.
- Renombrado automático de `genome_id` → nombre de muestra (BV-BRC etiqueta todos los genomas propios con el mismo nombre científico).

> **Recomendación:** con pocos aislados clonales el árbol aporta poca resolución. Conviene **añadir genomas de referencia** de BV-BRC (de la misma ST y de otras) para dar contexto y poder discriminante.

---

### Paso 9 — Genómica Comparativa

**Herramientas API:**
- `SeqComparison` (Proteome Comparison) vía BV-BRC

**Herramientas externas** (sin API integrada):
- [ANI — EZBioCloud](https://www.ezbiocloud.net/tools/ani)
- [KBase](https://narrative.kbase.us)

**Interpretación ANI:**

| Valor | Interpretación |
|---|---|
| ≥ 95% | Misma especie |
| ≥ 99% | Mismo clon / serotipo |
| < 95% | Especies distintas |

---

## 11. Reporte final HTML

Menú lateral → **"Generar Reporte"** → caja verde **"Reporte final completo (HTML)"**.

### 11.1 Indicador de etapas

Antes de generar, la app muestra el estado real de las 9 etapas según los datos importados:

```
✓ 1 · Control de calidad   ✓ 2 · Clasificación taxonómica   ✓ 3 · Ensamblado
✓ 4 · Anotación            ✓ 5 · Resistoma                  ✓ 6 · Alineamiento
✓ 7 · MLST                 ✓ 8 · Filogenia                  ✓ 9 · Comparativa

9 de 9 etapas con resultados
```

Se puede generar el reporte aunque falten etapas: las pendientes salen marcadas como tales.

### 11.2 Contenido del documento

| Sección | Incluye |
|---|---|
| Resumen | Proyecto, organismo, muestras, tabla de estado de las 9 etapas |
| Paso 1 · Calidad | Tabla + gráfica de lecturas por muestra |
| Paso 2 · Taxonomía | Tabla + gráfica de contaminación con umbral 5 % |
| Paso 3 · Ensamblado | Tabla + gráfica de contigs / N50 / profundidad frente a umbrales |
| Paso 4 · Anotación | Tabla + gráfica de CDS con función vs hipotéticos |
| Paso 5 · Resistoma | **Antibiograma S/I/R coloreado** + genes por clase + gráfica + clasificación MDR/XDR |
| Paso 6 · Alineamiento | Tabla o resumen de los genes alineados por BBH |
| Paso 7 · MLST | Tabla de ST y perfil alélico + gráfica de distribución |
| Paso 8 · Filogenia | **Árbol dibujado** con bootstrap y escala + Newick plegable |
| Paso 9 · Comparativa | Tabla y gráfica del pangenoma (core / accesorios / únicos) + ANI |
| Interpretación | Observaciones automáticas y notas del analista |
| Herramientas | Qué herramienta y plataforma se usó en cada etapa |
| Validación | Bloque de firma del responsable |

### 11.3 Interpretación automática

El reporte deriva observaciones de los propios datos:

- **Clonalidad** — si todas las muestras comparten ST.
- **Proporción del genoma core** — indicador de cercanía entre aislados.
- **Nº de clases de antibióticos** con genes de resistencia → MDR / XDR.
- ⚠️ **Avisos** de ensamblados fragmentados (> 200 contigs) y contaminación > 5 %, nombrando las muestras afectadas.

### 11.4 Reporte clásico (HTML / PDF / Word)

La caja superior de la misma pestaña mantiene el reporte configurable por secciones, con salida en HTML, PDF (requiere LaTeX) o Word.

> **Requisito Pandoc:** la app lo detecta automáticamente desde RStudio. Si falla la generación, abre la app **desde RStudio**.

---

## 12. Flujo de trabajo recomendado

### Flujo recomendado (pipeline automático)

```
INICIO — asistente de 3 pasos
  │
  ├─ 1. Seleccionar carpeta con los FASTQ   → muestras detectadas
  ├─ 2. Configurar credenciales             → BV-BRC ✓  Galaxy ✓
  └─ 3. Organismo + "Iniciar análisis completo"
         │
         ▼ Sube los FASTQ (una vez por muestra) y lanza la ola A
         │
         ▼ El motor revisa cada 45 s y encadena las etapas:
         │
         │   Calidad ─┐
         │   Taxonomía├─ inmediatas
         │   Ensamblado ──→ Anotación ──→ MLST
         │   Resistoma ─┘        └──→ Filogenia + Comparativa
         │                             (cuando TODAS las anotaciones terminan)
         │
         ▼ Los resultados se importan solos a tablas y gráficas
         │
         ▼ Puedes CERRAR LA APP — los jobs siguen en los servidores
         │      └─ al reabrir: "Cargar y reanudar" → se actualiza el avance
         │
         ▼ REPORTE FINAL
              → "Generar reporte final" → HTML con las 9 etapas
```

### Flujo paso a paso (manual por etapa)

```
Etapa N
  │  → Seleccionar herramienta y muestra (o "Enviar para todas las muestras")
  │  → Job enviado → spinner en barra lateral
  │  → Al completar: resultados importados automáticamente
  │  → Revisar tabla y gráfica; añadir notas
  └─ [✓ Etapa completada] → [Siguiente →]
```

### Si se cierra la app entre etapas

```
Reabrir la app
  │
  ├─ Modal automático: "Sesión guardada encontrada"
  │    ├─ "Restaurar sesión" → jobs activos reanudan monitoreo
  │    └─ "Nueva sesión"     → sincroniza igualmente con el workspace,
  │                             detecta lo ya analizado y NO lo repite
  │
  └─ O bien: Inicio → "Sesión de trabajo" → elegir sesión → "Cargar y reanudar"
       └─ consulta el avance real en BV-BRC y Galaxy e importa lo nuevo
```

---

## 13. Diagnóstico y solución de problemas

### 13.1 Diagnóstico BV-BRC (self-test)

**Menú lateral → Log / Consola → botón "Diagnóstico BV-BRC"**

Ejecuta la cadena completa contra tu sesión real y reporta cada paso con ✓/✗:

| # | Paso | Qué verifica |
|---|---|---|
| 1 | Sesión / token | Que hay token válido y qué dominio de workspace usa |
| 2 | Crear carpeta | `Workspace.create` de la carpeta de salida |
| 3 | Subir archivo | Subida real de un archivo de prueba por Shock |
| 4 | Listar workspace | `Workspace.ls` y que el archivo de prueba aparezca |
| 5 | Descargar archivo | `Workspace.get` y verificación del contenido |
| 6 | Listar resultados | Listado de la carpeta de salida de los jobs |

Cada paso queda también en el log con su detalle, de modo que se ve **en qué eslabón exacto** se rompe.

### 13.2 Consola de log

Pestaña **Log / Consola**. Registra cada llamada a las APIs con niveles de color:

| Nivel | Uso |
|---|---|
| `API` | Llamada saliente a una plataforma |
| `OK` | Operación exitosa |
| `INFO` | Información de progreso |
| `WARN` | Aviso — no bloquea |
| `ERROR` | Fallo — requiere atención |
| `DEBUG` | Detalle técnico |

Botones: **Limpiar**, **Descargar .log** (para compartir al pedir soporte) y **Refrescar**.

### 13.3 Problemas frecuentes

| Síntoma | Causa probable | Solución |
|---|---|---|
| Una etapa nunca muestra resultados | El job falló en el servidor | Revisa el log; busca `JobFailed` en el workspace de BV-BRC |
| `503 Service Unavailable` | Caída temporal de la API de BV-BRC | Esperar y reintentar; puede durar horas |
| Un job "exitoso" sin archivos | Faltó un parámetro obligatorio | Ver [Sección 14](#14-notas-técnicas-de-las-apis) |
| La taxonomía tarda muchísimo | Normal: Kraken2 en el Galaxy público | 30–60 min; no es un error |
| Error al generar el reporte | Pandoc no encontrado | Abrir la app **desde RStudio** |
| El pipeline no lanza taxonomía ni MLST | Sin API key de Galaxy | Configurar Galaxy en credenciales |
| Se re-analiza algo ya hecho | El resultado anterior falló o quedó vacío | Es el comportamiento esperado: solo se omiten resultados válidos |

> Tras 3 intentos fallidos de descarga, la app marca el job como `sin_resultados` y deja de consultarlo, para no saturar el log ni bloquear la interfaz.

---

## 14. Notas técnicas de las APIs

Detalles verificados contra los servicios reales. Útiles si se amplía o depura la app.

### 14.1 Identificadores de app en BV-BRC

Los nombres "comerciales" **no** son los IDs de la API (verificado con `AppService.enumerate_apps`). Enviar el ID equivocado hace que el job se acepte pero nunca produzca resultados:

| Nombre habitual | ID real |
|---|---|
| Assembly2 | `GenomeAssembly2` |
| Annotation | `GenomeAnnotation` |
| SeqComparison | `GenomeComparison` |
| PhylogeneticTree | `CodonTree` |
| CoreGenomeMLST | **no existe** → se usa Galaxy |

La traducción está centralizada en `BVBRC_APP_IDS` / `bvbrc_app_id()` (`R/api_bvbrc.R`).

### 14.2 Parámetros obligatorios fáciles de omitir

| App | Parámetro | Si falta |
|---|---|---|
| `FastqUtils` | `recipe` (lista, p. ej. `["fastqc"]`) | Corre unos segundos y **no genera nada** |
| `GenomeAnnotation` | `domain` (`Bacteria`/`Archaea`), `code` | Falla |
| `ComprehensiveGenomeAnalysis` | `input_type`, `domain` | Falla |
| `GenomeComparison` | `reference_genome_index` | `Not an ARRAY reference` |
| `GenomeComparison` | `genome_ids` debe ser **array**, no cadena | `Can't use string as ARRAY ref` |

### 14.3 Estructura de salida y descarga

Cada job crea **dos objetos** en la carpeta de salida:

- `<nombre>` — tipo `job_result`: JSON con `success` y `output_files`. **Se lee con `Workspace.get`, no se lista.**
- `.<nombre>` — carpeta oculta con los archivos reales (o `JobFailed.txt` si falló).

**Archivos grandes (Shock):** para archivos de cierto tamaño, `Workspace.get` devuelve una **URL de nodo Shock**, no el contenido. Hay que descargar aparte con `GET <url>?download` y cabecera `Authorization: OAuth <token>`. La app lo hace de forma transparente en `bvbrc_get_file()`.

### 14.4 Formatos de resultado por herramienta

| Herramienta | Archivo | Formato |
|---|---|---|
| FastqUtils | `*_fastqc.html` | HTML de FastQC (no MultiQC) |
| GenomeAssembly2 | `*_contigs.fasta` | FASTA con `length` y `coverage` en el encabezado |
| GenomeAnnotation | `*.features.txt` | TSV; el tipo de RNA va en la columna *function* |
| GenomeAnnotation | `amr-sir.txt` | Antibiograma predicho S/I/R con F1 |
| MetagenomicReadMapping | `kma.res` | KMA; plantilla `CARD\|<acc> <gen> [organismo]` |
| GenomeComparison | `genome_comparison.txt` | TSV; fila 1 = genomas, fila 2 = cabecera real |
| CodonTree | `*_tree.nwk` | Newick con `genome_id` (preferible para renombrar) |

### 14.5 Reintentos y estados obsoletos

Al reenviar un job con el **mismo** `output_file`, el objeto `job_result` conserva el resultado **anterior** hasta que la nueva tarea escribe. Para no leer un estado viejo hay que comparar el campo `id` del `job_result` con el `task_id` esperado.

---

## 15. Preguntas frecuentes

**¿Cómo ejecuto un análisis completo?**  
Tres pasos en la pantalla de Inicio: seleccionar la carpeta con los FASTQ → configurar credenciales → *Iniciar análisis completo*. Ver [Sección 3-bis](#3-bis-cómo-ejecutar-un-análisis-paso-a-paso) para el detalle.

**¿Puedo analizar varias cepas a la vez?**  
Sí. Todas las muestras que haya en la carpeta entran al pipeline; no hay que seleccionarlas una por una.

**¿Se guardan los datos si cierro la app?**  
**Sí**, y de tres formas: autoguardado interno cada 45 s, **sesión `.anbio`** al cerrar (para reanudar) y backup JSON manual (para archivar).

**¿Los jobs siguen corriendo si cierro la app?**  
Sí. Viven en los servidores de BV-BRC y Galaxy. La app solo guarda los identificadores de tarea; al reabrir los vuelve a consultar.

**¿Cómo retomo un análisis al día siguiente?**  
Inicio → *Sesión de trabajo* → elegir la sesión → **"Cargar y reanudar"**. La app consulta el avance real e importa lo que se completó mientras estuvo cerrada.

**¿Qué pasa si se apaga la PC a media corrida?**  
Nada se pierde. Los jobs siguen en los servidores; al reabrir se recupera el avance.

**¿Puedo mover el análisis a otro equipo?**  
Sí: *Descargar sesión* genera un `.anbio` que se abre en el otro equipo con *Abrir archivo…*. Los FASTQ hay que volver a indicarlos solo si se van a lanzar etapas nuevas.

**¿La app repite trabajo ya hecho?**  
No. Antes de enviar sincroniza con el workspace y omite lo que ya tiene resultados **válidos**. Los jobs fallidos o vacíos sí se reintentan.

**¿Los valores de Q30, GC, etc. se calculan automáticamente?**  
GC, longitud y número de lecturas se importan del reporte de FastQC. **El % Q30 no lo reporta FastQC** — según la guía LESP proviene de BaseSpace, por lo que se captura a mano en la pantalla de Inicio.

**¿Por qué la taxonomía usa Galaxy y no BV-BRC?**  
Porque el servicio `TaxonomicClassification` de BV-BRC falla de forma sistemática. Kraken2 en Galaxy es la ruta que funciona. Lo mismo ocurre con MLST: `CoreGenomeMLST` no existe como app en BV-BRC.

**¿Por qué la taxonomía tarda tanto?**  
Kraken2 en el Galaxy público carga una base de datos muy grande: 30–60 min es normal.

**¿Puedo ver el árbol filogenético dentro de la app?**  
Sí, en la pestaña Filogenia y en el reporte final: se dibuja desde el Newick con bootstrap, escala y las muestras propias en rojo.

**¿Por qué no aparece el botón "Descargar Reporte" después de generar?**  
Verifica que no haya un error rojo en pantalla. El error más común es Pandoc no encontrado — abre la app desde RStudio (no desde la consola de Windows).

**¿Puedo usar la app con organismos distintos a *Acinetobacter*?**  
Sí. El flujo es genérico para cualquier bacteria. Asegúrate de indicar el organismo correcto en la información del proyecto (se usa para los parámetros de anotación en BV-BRC) y selecciona el esquema MLST adecuado en el Paso 7.

**¿Por qué algunos círculos de la barra de progreso siguen grises aunque hay resultados?**  
El estado "completado" se activa al marcar el checkbox **"Etapa completada"** en la barra inferior de cada paso, o al presionar **"Siguiente →"**. Los badges de la barra lateral reflejan el estado de los jobs API; el stepper refleja la decisión del analista de dar por cerrada la etapa.

**¿Qué diferencia hay entre usar CGA y lanzar Assembly2 + MetagenomicReadMapping por separado?**  
`ComprehensiveGenomeAnalysis` (CGA) corre ensamblado, anotación y resistoma en un único job integral de BV-BRC. Es más completo pero puede tardar 2–4 horas. Lanzar las herramientas por separado es más rápido para obtener resultados parciales y es lo que hace el pipeline automático.

**Una etapa no muestra nada aunque el job terminó. ¿Qué reviso?**  
Ejecuta **Log / Consola → "Diagnóstico BV-BRC"**: prueba la cadena completa e indica en qué paso falla. Luego revisa el log; si el job falló en el servidor habrá un `JobFailed` en el workspace.

**¿Por qué algunos círculos de la barra de progreso siguen grises aunque hay resultados?**  
Los badges laterales reflejan el estado de los jobs; el stepper refleja que el analista dio por cerrada la etapa (checkbox *"Etapa completada"* o botón *Siguiente*). El indicador del **reporte final** sí se basa solo en si hay datos reales importados.

**¿La app inventa datos si una etapa no corrió?**  
No. Todas las tablas y gráficas quedan vacías con un aviso hasta que llegan resultados reales.

---

*Desarrollado por Vicente Esparza Villalpando, 2026.*  
*App construida con R Shiny, shinydashboard, httr2, plotly, DT y R Markdown.*
