#!/usr/bin/env bash
set -euo pipefail

PORT="${JUPYTER_PORT:-8888}"
DISABLE_TOKEN="${JUPYTER_DISABLE_TOKEN:-0}"
CUSTOM_TOKEN="${JUPYTER_TOKEN:-}"
EXTRA_ARGS="${JUPYTER_LAB_ARGS:-}"
PLAIN_URL="${JUPYTER_PLAIN_URL:-0}"    # If 1, print a single plain URL line for scripting
LINK_ONLY="${JUPYTER_LINK_ONLY:-0}"    # If 1, print only the OPEN THIS URL line (colored) then silence server output
DISABLE_LSP="${JUPYTER_DISABLE_LSP:-1}" # Disable jupyter_lsp extension (removes noisy skipped server list)
JUPYTER_UI="${JUPYTER_UI:-lab}"       # 'lab' (default) or 'notebook'
ROOT=/workspace
IP=0.0.0.0

# Persistence root (allows per-image version segregation)
PERSIST_TAG="${RSAFD_IMAGE_SHA:-unknown}"
PERSIST_BASE="${ROOT}/notebooks/.rsafd-docker-${PERSIST_TAG}"
PERSIST_PY="${PERSIST_BASE}/py"
PERSIST_R="${PERSIST_BASE}/R"

mkdir -p "${PERSIST_PY}" "${PERSIST_R}"

# Prepend custom Python user site (simple .pth via PYTHONPATH)
export PYTHONPATH="${PERSIST_PY}:${PYTHONPATH:-}"

# Configure pip to default to persistent dir if user runs '!pip install pkg'
export PIP_TARGET="${PERSIST_PY}"
export PIP_DISABLE_PIP_VERSION_CHECK=1

# R user library path (so install.packages() without lib= writes here and survives rebuild)
export R_LIBS_USER="${PERSIST_R}"
mkdir -p "${R_LIBS_USER}"

# Helpful MOTD line (only when banner mode)
PERSIST_NOTE="Custom Python -> ${PERSIST_PY} | R -> ${PERSIST_R}"

# Colors
B="\e[1m"; G="\e[32m"; Y="\e[33m"; C="\e[36m"; R="\e[31m"; NC="\e[0m"

if [[ "$DISABLE_TOKEN" == "1" ]]; then
  TOKEN_ARG=(--IdentityProvider.token='' --ServerApp.token='') # include old flag for backward compat
  DISPLAY_TOKEN_MSG="${Y}NO TOKEN (DISABLED)${NC}"
  TOKEN_QUERY=""
else
  # If user did not supply a token, generate a cryptographically strong one
  if [[ -z "$CUSTOM_TOKEN" ]]; then
    if command -v openssl >/dev/null 2>&1; then
      CUSTOM_TOKEN=$(openssl rand -hex 24)
    else
      CUSTOM_TOKEN=$(python3 - <<'PY'
import secrets; print(secrets.token_hex(24))
PY
)
    fi
  fi
  TOKEN_ARG=(--IdentityProvider.token="${CUSTOM_TOKEN}")
  DISPLAY_TOKEN_MSG="${G}${CUSTOM_TOKEN}${NC}"
  TOKEN_QUERY="?token=${CUSTOM_TOKEN}"
fi

ARGS=(
  --ServerApp.ip=${IP}
  --ServerApp.port=${PORT}
  --ServerApp.root_dir="${ROOT}"
  --ServerApp.open_browser=False
  --ServerApp.allow_origin='*'
  --ServerApp.allow_remote_access=True
  --ServerApp.password=''
  "${TOKEN_ARG[@]}"
)

if [[ "$JUPYTER_UI" == "notebook" ]]; then
  # Add classic NotebookApp token flags for compatibility
  if [[ "$DISABLE_TOKEN" == "1" ]]; then
    ARGS+=(--NotebookApp.token='')
  else
    ARGS+=(--NotebookApp.token="${CUSTOM_TOKEN}")
  fi
fi

if [[ "$DISABLE_LSP" == "1" ]]; then
  # Create a minimal config overriding extensions
  export JUPYTER_CONFIG_DIR="${HOME}/.jupyter"
  mkdir -p "${JUPYTER_CONFIG_DIR}"
  python3 - <<'PY'
import json,os,io
cfg_dir=os.path.join(os.environ.get('JUPYTER_CONFIG_DIR',''), 'labconfig')
os.makedirs(cfg_dir, exist_ok=True)
# jupyter_lsp registers as a server extension; disable via overrides.json
override_path=os.path.join(cfg_dir,'overrides.json')
data={'disabledExtensions':['@jupyter-server/jupyter-lsp-extension']}
with open(override_path,'w') as f: json.dump(data,f)
PY
fi

# Add any user supplied additional args
if [[ -n "$EXTRA_ARGS" ]]; then
  # shellcheck disable=SC2206
  ARGS+=( $EXTRA_ARGS )
fi

HOSTNAME_DISPLAY=$(hostname)

if [[ "$JUPYTER_UI" == "notebook" ]]; then
  URL_PATH="/tree"
else
  URL_PATH="/lab"
fi

LOCAL_URL="http://127.0.0.1:${PORT}${URL_PATH}${TOKEN_QUERY}"
CONTAINER_URL="http://${HOSTNAME_DISPLAY}:${PORT}${URL_PATH}${TOKEN_QUERY}"

if [[ "$LINK_ONLY" == "1" ]]; then
  # Only show the colored OPEN THIS URL line; silence subsequent server logs
  echo -e "${B}OPEN THIS URL:${NC} ${LOCAL_URL}"
elif [[ "$PLAIN_URL" == "1" ]]; then
  # Plain single line (no colors) for copy/paste or automation
  if [[ -n "$TOKEN_QUERY" ]]; then
    echo "${LOCAL_URL}"
  else
    echo "${LOCAL_URL}" # no token case
  fi
else
  cat <<EOB
${C}┌──────────────────────────────────────────────────────────────┐${NC}
${C}│${NC} ${B}JupyterLab starting${NC}
${C}│${NC} Host:     ${HOSTNAME_DISPLAY}
${C}│${NC} Port:     ${PORT}
${C}│${NC} Mode:     ${JUPYTER_UI}
${C}│${NC} Token:    ${DISPLAY_TOKEN_MSG}
${C}│${NC} Disable:  JUPYTER_DISABLE_TOKEN=${DISABLE_TOKEN}
${C}│${NC} LSP off:  JUPYTER_DISABLE_LSP=${DISABLE_LSP}
${C}│${NC} Python:   $(python3 --version 2>/dev/null | awk '{print $2}')
${C}│${NC} R:        $(R --version 2>/dev/null | awk 'NR==1{print $3}')
${C}│${NC} Persist:  ${PERSIST_NOTE}
${C}└──────────────────────────────────────────────────────────────┘${NC}
${B}OPEN THIS URL:${NC} ${LOCAL_URL}
${B}Alt (container):${NC} ${CONTAINER_URL}
EOB
fi

if [[ "$DISABLE_TOKEN" == "1" && "$LINK_ONLY" != "1" && "$PLAIN_URL" != "1" ]]; then
  echo -e "${Y}WARNING:${NC} Token disabled – anyone with port access can use the notebook." >&2
fi

if [[ "$JUPYTER_UI" == "notebook" ]]; then
  SUBCOMMAND="notebook"
else
  SUBCOMMAND="lab"
fi

# Launch Jupyter with the correct environment variables for the kernels
if [[ "$LINK_ONLY" == "1" ]]; then
  exec env "PYTHONPATH=${PYTHONPATH}" "PIP_TARGET=${PIP_TARGET}" "R_LIBS_USER=${R_LIBS_USER}" jupyter "$SUBCOMMAND" "${ARGS[@]}" >/dev/null 2>&1
else
  exec env "PYTHONPATH=${PYTHONPATH}" "PIP_TARGET=${PIP_TARGET}" "R_LIBS_USER=${R_LIBS_USER}" jupyter "$SUBCOMMAND" "${ARGS[@]}"
fi
