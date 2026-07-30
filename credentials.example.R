# ============================================================
# ANbio — plantilla de credenciales
# ============================================================
#
# CÓMO USARLA
#   1. Copia este archivo y renómbralo como  credentials.R
#   2. Rellena tus datos reales
#   3. NO subas credentials.R al repositorio (ya está en .gitignore)
#
# DÓNDE OBTENER LAS CLAVES
#   BV-BRC : cuenta gratuita en https://www.bv-brc.org
#   Galaxy : https://usegalaxy.org → User → Preferences → Manage API Key
#
# La app hace auto-login con estos valores al arrancar. Si prefieres no
# guardarlos en disco, deja las cadenas vacías: podrás capturarlos desde
# la propia app (botón "Configurar credenciales"), y no se escribirán
# en ningún archivo.
# ============================================================

# ---- BV-BRC ----
CRED_BVBRC_USER <- ""      # usuario o correo de BV-BRC
CRED_BVBRC_PASS <- ""      # contraseña

# ---- Galaxy (usegalaxy.org) ----
CRED_GALAXY_KEY <- ""      # API key

# ------------------------------------------------------------
# Las siguientes plataformas son públicas y NO requieren clave:
#   · EBI JDispatcher (MUSCLE, Clustal Omega, MAFFT)
#   · NCBI BLAST
#   · CGE / ResFinder (DTU)
# ------------------------------------------------------------
