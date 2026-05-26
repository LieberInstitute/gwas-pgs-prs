#!/usr/bin/env bash
set -euo pipefail

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

msg() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}

require_file() {
  local path=$1
  if [ ! -s "$path" ]; then
    if [ "${DRY_RUN:-0}" = "1" ]; then
      msg "dry-run warning: missing required file: $path"
      return 0
    fi
    die "missing required file: $path"
  fi
}

require_dir() {
  local path=$1
  [ -d "$path" ] || die "missing required directory: $path"
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "missing command: $1"
}

source_config() {
  local path=$1
  require_file "$path"
  ## shellcheck source=/dev/null
  source "$path"
}

run_cmd() {
  local dry_run=${DRY_RUN:-0}
  printf '+'
  printf ' %q' "$@"
  printf '\n'
  if [ "$dry_run" = "1" ]; then
    return 0
  fi
  "$@"
}

run_bash() {
  local dry_run=${DRY_RUN:-0}
  local cmd=$1
  printf '+ %s\n' "$cmd"
  if [ "$dry_run" = "1" ]; then
    return 0
  fi
  bash -o pipefail -c "$cmd"
}

manifest_rows() {
  local manifest=$1
  require_file "$manifest"
  awk -F '\t' 'NR > 1 && $0 !~ /^[[:space:]]*($|#)/' "$manifest"
}

field_index() {
  local manifest=$1
  local field=$2
  awk -v f="$field" -F '\t' 'NR == 1 {
    for (i = 1; i <= NF; i++) if ($i == f) { print i; exit }
  }' "$manifest"
}

require_fields() {
  local manifest=$1
  shift
  local field
  for field in "$@"; do
    [ -n "$(field_index "$manifest" "$field")" ] || die "$manifest lacks required column: $field"
  done
}
