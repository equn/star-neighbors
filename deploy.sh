#!/usr/bin/env bash
# Deploy: star-neighbors/web/  ->  NY server  /www/wwwroot/equn.me/star-neighbors/
# Public URL: https://equn.me/star-neighbors/
#
# Usage:  ./deploy.sh [-y] [--no-purge] [-n|--dry-run] [--skip-checks]
#
# What is published
# -----------------
# Only `web/` is a website. The repository root is a build tree: scripts/,
# data/processed/ (~44 MB of regenerable intermediates) and reference/ have no
# business on a web server. So this stages a copy rather than pointing rsync at
# the project root, and adds two files that live outside web/:
#
#   data/star_encounters.csv  the catalog deliverable - the in-app Methods
#   data/COLUMNS.md           panel points at "data/star_encounters.csv", which
#                             would be a dangling reference otherwise
#   README.md                 the methods write-up the app defers to
#
# Preflight
# ---------
# The pipeline has hard gates (scripts/05_validate.py, 08_verify.py) precisely so
# a half-built catalog cannot be published. Skipping them at the last step
# would defeat the point, so this refuses to deploy unless:
#
#   * 08_verify passes - the CSV, web bundle, coefficients and metadata all
#     describe the same catalog;
#   * the embedded validation status is "pass";
#   * web/data/ is no older than data/processed/orbits.csv, i.e. the bundle was
#     exported from the science currently on disk, not from a previous run.
#
# --skip-checks exists for republishing an unchanged site when the build tree is
# absent. It prints a loud warning; do not use it to push past a failure.
#
# Cache busting
# -------------
# app.js is referenced with a content hash (app.js?v=...) rewritten into the
# staged index.html only, so the source tree stays clean. Without it a returning
# visitor could run an old app.js against a new data bundle - the data itself is
# fetched with `cache: 'no-cache'` and always revalidates, but that policy lives
# *inside* app.js and cannot protect the file that carries it. The Cloudflare
# purge below handles the CDN layer; this handles browsers.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/../vpshosts/deploy-lib.sh"
_require_secrets

APP_NAME="star-neighbors"
SSH_HOST="ny"
REMOTE_TARGET="/www/wwwroot/equn.me/star-neighbors/"
CF_ZONE_ID="${CF_ZONE_EQUN_ME:-}"

EXCLUDES=(
  "${DEFAULT_EXCLUDES[@]}"
)

SKIP_CHECKS=false
DRY_RUN=false
ARGS=()
for a in "$@"; do
  case "$a" in
    --skip-checks) SKIP_CHECKS=true ;;
    -n|--dry-run)  DRY_RUN=true; ARGS+=("$a") ;;
    *) ARGS+=("$a") ;;
  esac
done

# --------------------------------------------------------------------- checks
preflight() {
  local proc="$SCRIPT_DIR/data/processed"
  local bundle="$SCRIPT_DIR/web/data"

  for f in "$bundle/meta.json" "$bundle/stars.json" "$bundle/cheb.bin" \
           "$SCRIPT_DIR/web/index.html" "$SCRIPT_DIR/web/app.js"; do
    [ -f "$f" ] || { _err "missing from the bundle: ${f#$SCRIPT_DIR/}"; exit 1; }
  done

  if [ ! -f "$proc/orbits.csv" ]; then
    _warn "data/processed/orbits.csv absent - cannot re-verify the build."
    _warn "Publishing the existing bundle unchecked."
    return 0
  fi

  _log "Verifying the build before publishing..."
  if ! python3 "$SCRIPT_DIR/scripts/08_verify.py" >/dev/null 2>&1; then
    _err "08_verify.py failed. Run it directly to see which artifact disagrees."
    exit 1
  fi

  python3 - "$SCRIPT_DIR" <<'PY' || exit 1
import json, pathlib, sys
root = pathlib.Path(sys.argv[1])
meta = json.loads((root / "web/data/meta.json").read_text())
val = meta.get("validation", {})
if val.get("status") != "pass":
    print(f"ERROR: embedded validation status is {val.get('status', 'absent')!r} "
          f"({val.get('n_failed', '?')} of {val.get('n_checks', '?')} checks failed)",
          file=sys.stderr)
    sys.exit(1)

# The bundle must not predate the science it claims to represent.
orbits = (root / "data/processed/orbits.csv").stat().st_mtime
stale = [p.name for p in (root / "web/data").iterdir()
         if p.is_file() and p.stat().st_mtime < orbits - 1]
if stale:
    print(f"ERROR: {stale} are older than data/processed/orbits.csv.",
          file=sys.stderr)
    print("       Re-run scripts/06_export.py (and 07_csv.py) before deploying.",
          file=sys.stderr)
    sys.exit(1)

print(f"    validation pass ({val['n_checks']} checks) · "
      f"{meta['n_stars_shipped']} stars shipped of {meta['n_stars_catalog']} · "
      f"retrieved {meta.get('retrieved', '?')}")
PY
}

# -------------------------------------------------------------------- staging
STAGE="$(mktemp -d "${TMPDIR:-/tmp}/star-neighbors-deploy.XXXXXX")"
trap 'rm -rf "$STAGE"' EXIT

build_stage() {
  _log "Staging the publishable tree..."
  rsync -a --exclude='.DS_Store' "$SCRIPT_DIR/web/" "$STAGE/"

  # The catalog deliverable, so the Methods panel's reference resolves.
  for f in "data/star_encounters.csv" "data/COLUMNS.md" "README.md"; do
    if [ -f "$SCRIPT_DIR/$f" ]; then
      cp "$SCRIPT_DIR/$f" "$STAGE/$(basename "$f")"
      [ "$f" = "data/star_encounters.csv" ] && mv "$STAGE/star_encounters.csv" "$STAGE/data/"
      [ "$f" = "data/COLUMNS.md" ] && mv "$STAGE/COLUMNS.md" "$STAGE/data/"
    else
      _warn "not staged (absent): $f"
    fi
  done

  # Fingerprint app.js in the staged index.html only.
  local hash
  hash=$(shasum -a 256 "$STAGE/app.js" | cut -c1-12)
  if grep -q 'src="app\.js"' "$STAGE/index.html"; then
    sed -i '' "s|src=\"app\.js\"|src=\"app.js?v=$hash\"|" "$STAGE/index.html"
    _log "    app.js fingerprinted as app.js?v=$hash"
  else
    _warn "could not fingerprint app.js - the <script src=\"app.js\"> tag has changed"
  fi

  _log "    staged $(find "$STAGE" -type f | wc -l | tr -d ' ') files, $(du -sh "$STAGE" | cut -f1)"
}

if [ "$SKIP_CHECKS" = true ]; then
  _warn "--skip-checks: publishing without re-verifying the build."
else
  preflight
fi
build_stage

LOCAL_DIR="$STAGE/"
deploy_app "${ARGS[@]+"${ARGS[@]}"}"

if [ "$DRY_RUN" = false ]; then
  printf '\n'
  _log "Live at https://equn.me/star-neighbors/"
fi
