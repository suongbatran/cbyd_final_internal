#!/usr/bin/env bash
# Launch the reaction-route dashboard — VIEW (read-only) build.
#
#   ./scripts/dashboard/to_view/run.sh              # http://<this-host>:8001  (shared: binds 0.0.0.0)
#   ./scripts/dashboard/to_view/run.sh 8080         # pick a different port
#   HOST=127.0.0.1 ./scripts/dashboard/to_view/run.sh   # solo / localhost only
#
# Read-only: browse routes and inspect reactions/molecules. No editing, review, or
# staging — the original folder is only ever read. Safe to expose more broadly.
#
# Choose which folder to browse (default: data/routes):
#   ./scripts/dashboard/to_view/run.sh --routes data/cbyd_extra_batch
#   ./scripts/dashboard/to_view/run.sh 8080 --routes data/cbyd_extra_batch
# Path may be relative to the repo root or absolute. (Env var: DASH_ROUTES_DIR.)
#
# Keep it alive across logout with tmux: `tmux new -s dashview`. Ctrl+C to stop.
set -euo pipefail

PORT="8001"
HOST="${HOST:-0.0.0.0}"          # 0.0.0.0 = reachable on the network; 127.0.0.1 = local only

# Parse args: first bare number = PORT; --routes sets the folder to browse.
while [ $# -gt 0 ]; do
  case "$1" in
    --routes) export DASH_ROUTES_DIR="$2"; shift 2;;
    *) PORT="$1"; shift;;
  esac
done

# Run from the repo root so the `scripts` package and data/ paths resolve.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT"

# Activate the conda env if it isn't already active.
if [ "${CONDA_DEFAULT_ENV:-}" != "cbyd_test" ]; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate cbyd_test
fi

HOSTNAME_SHORT="$(hostname -s 2>/dev/null || hostname)"
echo "Dashboard (VIEW · read-only) starting on ${HOST}:${PORT}"
echo "  original dir: ${DASH_ROUTES_DIR:-data/routes}"
if [ "$HOST" = "0.0.0.0" ]; then
  echo "  open from another machine →  http://${HOSTNAME_SHORT}:${PORT}"
  IP="$(hostname -I 2>/dev/null | awk '{print $1}')"
  [ -n "${IP:-}" ] && echo "  or by IP                  →  http://${IP}:${PORT}"
else
  echo "  open →  http://${HOST}:${PORT}"
fi
# Optional login gate (same as the edit build): DASH_USERS / DASH_PASSCODE.
if [ -n "${DASH_USERS:-}" ]; then
  echo "  login: named accounts ON"
elif [ -n "${DASH_PASSCODE:-}" ]; then
  echo "  login: shared passcode ON"
fi
echo "  (Ctrl+C to stop)"

exec uvicorn scripts.dashboard.to_view.app:app --host "$HOST" --port "$PORT" \
  --reload --reload-dir scripts/dashboard/to_view
