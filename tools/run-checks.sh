#!/usr/bin/env bash
# No dependency downloads. Supply a prepared, pinned R environment.
# Example: R_BIN=/path/to/R RSCRIPT_BIN=/path/to/Rscript \
#   tools/run-checks.sh --baseline-source /path/to/pinned-source \
#   --output-dir /path/to/new-receipts --mode full --timeout-seconds 900
set -u
set -o pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
candidate_source=$(cd -- "$script_dir/.." && pwd)
baseline_source=""
output_dir=""
mode="full"
timeout_seconds=900
while (($#)); do
  if (($# < 2)); then printf 'Each option needs a value\n' >&2; exit 2; fi
  case "$1" in
    --baseline-source) baseline_source=$2 ;;
    --output-dir) output_dir=$2 ;;
    --mode) mode=$2 ;;
    --timeout-seconds) timeout_seconds=$2 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
  esac
  shift 2
done
if [[ -z "$baseline_source" || -z "$output_dir" || ! -f "$baseline_source/DESCRIPTION" ]]; then
  printf 'Required: --baseline-source PATH --output-dir NEW_PATH\n' >&2
  exit 2
fi
if [[ "$mode" != full && "$mode" != quick ]] || [[ ! "$timeout_seconds" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Mode must be full/quick and timeout a positive integer in seconds\n' >&2
  exit 2
fi
if [[ -e "$output_dir" ]]; then
  printf 'Output path already exists; use a new path to preserve prior receipts: %s\n' "$output_dir" >&2
  exit 2
fi
baseline_source=$(cd -- "$baseline_source" && pwd)
mkdir -p -- "$output_dir" || exit 2
output_dir=$(cd -- "$output_dir" && pwd)
r_bin=${R_BIN:-R}
rscript_bin=${RSCRIPT_BIN:-Rscript}
timeout_bin=${TIMEOUT_BIN:-timeout}
for program in "$r_bin" "$rscript_bin" "$timeout_bin"; do
  if ! command -v -- "$program" >/dev/null 2>&1; then
    printf 'NOT RUN: required executable unavailable: %s\n' "$program" | tee -a "$output_dir/setup.log"
    exit 2
  fi
done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
printf 'stage\tstatus\texit_code\telapsed_seconds\n' > "$output_dir/status.tsv"
overall=0
run_logged() {
  local label=$1
  shift
  local start=$SECONDS rc
  printf 'START %s\n' "$label"
  if [[ -x /usr/bin/time ]]; then
    /usr/bin/time -v -o "$output_dir/$label.resources.txt" \
      "$timeout_bin" --signal=TERM --kill-after=15 "$timeout_seconds" "$@" > "$output_dir/$label.log" 2>&1
    rc=$?
  else
    "$timeout_bin" --signal=TERM --kill-after=15 "$timeout_seconds" "$@" > "$output_dir/$label.log" 2>&1
    rc=$?
    printf 'NOT RUN: process peak RSS, /usr/bin/time unavailable\n' > "$output_dir/$label.resources.txt"
  fi
  local status=PASS
  if ((rc != 0)); then status=FAIL; overall=1; fi
  if ((rc == 124 || rc == 137)); then status=TIMEOUT; fi
  if [[ "$label" == *-check && "$rc" == 0 ]]; then
    while IFS= read -r line; do
      if [[ "$line" == "Status: "* && "$line" != "Status: OK" ]]; then
        status=COMPLETED_WITH_FINDINGS
        printf '%s\n' "$line" > "$output_dir/$label.findings.txt"
      fi
    done < "$output_dir/$label.log"
  fi
  printf '%s\t%s\t%s\t%s\n' "$label" "$status" "$rc" "$((SECONDS-start))" >> "$output_dir/status.tsv"
  printf '%s %s (exit %s)\n' "$status" "$label" "$rc"
  return "$rc"
}
skip_stage() { printf '%s\tNOT RUN\tNA\tNA\n' "$1" >> "$output_dir/status.tsv"; }
"$r_bin" --version > "$output_dir/R-version.txt" 2>&1
"$rscript_bin" --vanilla -e 'writeLines(capture.output(sessionInfo())); print(.libPaths()); ip <- installed.packages(); write.table(ip[,c("Package","Version","LibPath")],row.names=FALSE,sep="\t")' > "$output_dir/environment.txt" 2>&1
for label in baseline candidate; do
  if [[ "$label" == baseline ]]; then src=$baseline_source; else src=$candidate_source; fi
  lib="$output_dir/$label-library"
  mkdir -p -- "$lib" "$output_dir/$label-build"
  git -C "$src" rev-parse HEAD > "$output_dir/$label-commit.txt" 2>/dev/null || printf 'UNAVAILABLE\n' > "$output_dir/$label-commit.txt"
  git -C "$src" diff --binary > "$output_dir/$label-working-tree.patch" 2>/dev/null || true
  export DRFARM_SOURCE_COMMIT=$(cat "$output_dir/$label-commit.txt")
  if run_logged "$label-install" "$r_bin" CMD INSTALL --preclean --no-multiarch --library="$lib" "$src"; then
    compare_args=()
    if [[ "$label" == candidate && -f "$output_dir/baseline-example/results.rds" ]]; then
      compare_args=(--compare "$output_dir/baseline-example/results.rds")
    fi
    run_logged "$label-example" "$rscript_bin" --vanilla "$script_dir/reproduce.R" \
      --output "$output_dir/$label-example" --library "$lib" --mode "$mode" --label "$label" "${compare_args[@]}" || true
    if [[ "$label" == candidate ]]; then
      run_logged candidate-profile "$rscript_bin" --vanilla "$script_dir/profile.R" \
        --output "$output_dir/candidate-profile" --library "$lib" --kinship yes || true
    fi
  else
    skip_stage "$label-example"
    if [[ "$label" == candidate ]]; then skip_stage candidate-profile; fi
  fi
  # R CMD build writes its tarball to the working directory, never the source.
  previous_dir=$PWD
  cd -- "$output_dir/$label-build" || exit 2
  if run_logged "$label-build" "$r_bin" CMD build --no-manual "$src"; then
    tarballs=(./drfarm_*.tar.gz)
    if ((${#tarballs[@]} == 1)) && [[ -f "${tarballs[0]}" ]]; then
      run_logged "$label-check" env R_LIBS="$lib${R_LIBS:+:$R_LIBS}" \
        "$r_bin" CMD check --no-manual "${tarballs[0]}" || true
    else
      printf 'Expected exactly one built drfarm tarball\n' >> "$output_dir/$label-build.log"
      skip_stage "$label-check"
      overall=1
    fi
  else
    skip_stage "$label-check"
  fi
  cd -- "$previous_dir" || exit 2
done
cat "$output_dir/status.tsv"
exit "$overall"
