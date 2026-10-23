#!/usr/bin/env bash
# =============================================================================
# Download raw FASTQs + ENA metadata records for prenatal human forebrain
# RNA-seq from E-MTAB-6814 (Cardoso-Moreira et al. 2019, Nature 571:505).
#
# Selection: Homo sapiens, organism part = forebrain (prenatally this is
#            forebrain + midbrain), stage = "N week post conception", both sexes.
#            32 runs, single-end 100 nt, TruSeq Stranded mRNA (reverse strand).
#
# Usage:
#   bash download_prenatal_forebrain.sh [OUTDIR]
#
# Output layout:
#   OUTDIR/samples.tsv               sample sheet (run, name, biosample, wpc, sex)
#   OUTDIR/ena/read_run.tsv          ENA run report: accessions, MD5s, sizes,
#                                    read/base counts, submitted file names, ...
#   OUTDIR/ena/runs.xml              ENA XML records for the runs (ERR)
#   OUTDIR/ena/experiments.xml       ENA XML records for the experiments (ERX)
#   OUTDIR/ena/samples.xml           ENA XML records for the samples (ERS)
#   OUTDIR/fastq/<run>_<name>.fastq.gz
#   OUTDIR/download-fastq.log       script output and errors
#
# Environment options:
#   JOBS=2        parallel FASTQ downloads (default 2; EBI throttles heavy use)
#   PROTO=https   https (default) or ftp
#   LIMIT_RATE=   optional per-download bandwidth cap, e.g. 20M
#   DRY_RUN=1     print what would be fetched, download nothing
#   SKIP_MD5=1    skip MD5 verification
#
# Rate-limit handling:
#   * All ENA metadata comes from 4 batched requests, made once and cached.
#     MD5s are read from that cache, not requested per file.
#   * Only JOBS downloads run at once, and each start is staggered.
#   * curl retries 429/5xx/timeouts with exponential backoff and honours
#     Retry-After; on top of that each file gets up to 3 resumed attempts
#     with a growing pause between them.
#   * Re-running is safe: verified files are skipped, partial files resume.
# =============================================================================
set -euo pipefail

OUTDIR="${1:-fastq_prenatal_forebrain}"
JOBS="${JOBS:-2}"
PROTO="${PROTO:-https}"
LIMIT_RATE="${LIMIT_RATE:-}"
DRY_RUN="${DRY_RUN:-0}"
SKIP_MD5="${SKIP_MD5:-0}"

BASE="${PROTO}://ftp.sra.ebi.ac.uk/vol1/fastq"
PORTAL="https://www.ebi.ac.uk/ena/portal/api"
BROWSER="https://www.ebi.ac.uk/ena/browser/api"
ENA_DIR="${OUTDIR}/ena"
FQ_DIR="${OUTDIR}/fastq"
META="${ENA_DIR}/read_run.tsv"
LOG_FILE="${OUTDIR}/download-fastq.log"

mkdir -p "$ENA_DIR" "$FQ_DIR"
exec >>"$LOG_FILE" 2>&1

# Shared curl options. --retry covers 408/429/500/502/503/504 and timeouts;
# without --retry-delay curl backs off exponentially and obeys Retry-After.
CURL_OPTS=(-fL --silent --show-error --retry 8 --retry-connrefused
           --connect-timeout 30)

# -----------------------------------------------------------------------------
# Sample sheet (from the E-MTAB-6814 SDRF).
# wpc = "developmental stage" characteristic. Two samples are named "20w" but
# annotated 19 wpc (5668sTS, 5770sTS); the annotation is used here.
# -----------------------------------------------------------------------------
SAMPLE_SHEET="${OUTDIR}/samples.tsv"
cat > "$SAMPLE_SHEET" <<'EOF'
run	sample_name	biosample	wpc	sex
ERR2598164	3587sTS.Human.Forebrain.CS13.Male	SAMEA4690086	4	male
ERR2598165	5624sTS.Human.Forebrain.4w.Male	SAMEA4690087	4	male
ERR2598166	6096sTS.Human.Forebrain.4w.Male	SAMEA4690088	4	male
ERR2598180	3669sTS.Human.Forebrain.CS16.Male	SAMEA4690102	5	male
ERR2598181	5628sTS.Human.Forebrain.5w.Female	SAMEA4690103	5	female
ERR2598206	2026sTSm.Human.Forebrain.CS21.Male	SAMEA4690128	7	male
ERR2598207	3571sTS.Human.Forebrain.CS20.Female	SAMEA4690129	7	female
ERR2598208	3636sTS.Human.Forebrain.CS21.Male	SAMEA4690130	7	male
ERR2598226	2034sTS.Human.Forebrain.CS22.Female	SAMEA4690148	8	female
ERR2598227	2080sTS.Human.Brain.CS23.Male	SAMEA4690149	8	male
ERR2598228	2086sTSm.Human.Brain.CS23.Female	SAMEA4690150	8	female
ERR2598229	2035sTSm.Human.Forebrain.8w.Male	SAMEA4690151	8	male
ERR2598250	2041sTSb.Human.Brain.9w.Male	SAMEA4690173	9	male
ERR2598251	2050sTS.Human.Forebrain.9w.Female	SAMEA4690174	9	female
ERR2598252	2068sTS.Human.Brain.9w.Male	SAMEA4690175	9	male
ERR2598051	5664sTS.Human.Forebrain.10w.Female	SAMEA4689973	10	female
ERR2598069	2047sTS.Human.Forebrain.11w.Male	SAMEA4689991	11	male
ERR2598070	2111sTS.Human.Forebrain.11w.Male	SAMEA4689992	11	male
ERR2598071	5952sTS.Human.Forebrain.11w.Male	SAMEA4689993	11	male
ERR2598072	5953sTS.Human.Forebrain.11w.Female	SAMEA4689994	11	female
ERR2598089	5721sTS.Human.Forebrain.12w.Male	SAMEA4690011	12	male
ERR2598102	5711sTS.Human.Forebrain.13w.Male	SAMEA4690024	13	male
ERR2598103	5939sTS.Human.Forebrain.13w.Male	SAMEA4690025	13	male
ERR2598104	5955sTS.Human.Forebrain.13w.Female	SAMEA4690026	13	female
ERR2598120	5665sTSm.Human.Forebrain.16w.Female	SAMEA4690042	16	female
ERR2598121	5951sTS.Human.Forebrain.16w.Male	SAMEA4690043	16	male
ERR2598137	5667sTS.Human.Forebrain.18w.Male	SAMEA4690059	18	male
ERR2598148	5768sTS.Human.Brain.19w.Male	SAMEA4690070	19	male
ERR2598149	5782sTS.Human.Brain.19w.Female	SAMEA4690071	19	female
ERR2598150	5960sTS.Human.Brain.19w.Male	SAMEA4690072	19	male
ERR2598151	5668sTS.Human.Brain.20w.Male	SAMEA4690073	19	male
ERR2598152	5770sTS.Human.Brain.20w.Female	SAMEA4690074	19	female
EOF

N_RUNS=$(( $(wc -l < "$SAMPLE_SHEET") - 1 ))
RUN_LIST=$(tail -n +2 "$SAMPLE_SHEET" | cut -f1 | paste -sd, -)

# -----------------------------------------------------------------------------
# 1. ENA metadata: 4 batched requests, cached. Skipped if already present.
# -----------------------------------------------------------------------------
FIELDS="run_accession,experiment_accession,sample_accession,secondary_sample_accession,study_accession,secondary_study_accession,sample_alias,experiment_title,library_name,library_layout,library_strategy,library_source,library_selection,instrument_model,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp,submitted_ftp,submitted_md5,first_public"

column_of() {  # print column $2 (by header name) of TSV $1, comma-joined
  awk -F'\t' -v col="$2" 'NR==1{for(i=1;i<=NF;i++) if($i==col) c=i; next}
                           c && $c!="" {print $c}' "$1" | sort -u | paste -sd, -
}

fetch_xml() {  # $1 = comma-separated accessions, $2 = output file
  local acc="$1" out="$2"
  if [[ -s "$out" ]]; then echo "[skip] $(basename "$out") already cached"; return 0; fi
  if [[ "$DRY_RUN" == "1" ]]; then echo "[dry-run] ${BROWSER}/xml/<${acc%%,*},...> -> ${out}"; return 0; fi
  if [[ -z "$acc" ]]; then echo "[warn] no accessions for $(basename "$out"); skipped" >&2; return 0; fi
  echo "[meta] $(basename "$out")"
  curl "${CURL_OPTS[@]}" -o "${out}.part" "${BROWSER}/xml/${acc}" \
    && mv "${out}.part" "$out" \
    || { rm -f "${out}.part"; echo "[warn] could not fetch $(basename "$out")" >&2; }
  sleep 2   # be gentle between metadata calls
}

meta_rows() { [[ -s "$META" ]] && echo $(( $(wc -l < "$META") - 1 )) || echo 0; }

if [[ $(meta_rows) -eq $N_RUNS ]]; then
  echo "[skip] read_run.tsv already cached"
elif [[ "$DRY_RUN" == "1" ]]; then
  echo "[dry-run] ${PORTAL}/search?result=read_run&includeAccessions=<${N_RUNS} runs> -> ${META}"
else
  echo "[meta] read_run.tsv (one request for all ${N_RUNS} runs)"
  if ! curl "${CURL_OPTS[@]}" -o "${META}.part" -X POST "${PORTAL}/search" \
         --data-urlencode "result=read_run" \
         --data-urlencode "includeAccessions=${RUN_LIST}" \
         --data-urlencode "fields=${FIELDS}" \
         --data-urlencode "format=tsv" \
         --data-urlencode "limit=0" \
     || [[ $(( $(wc -l < "${META}.part") - 1 )) -ne $N_RUNS ]]; then
    # Fallback: one request per run, paced at 1 per second.
    echo "[meta] batched query failed or incomplete; falling back to per-run requests" >&2
    : > "${META}.part"
    for run in ${RUN_LIST//,/ }; do
      curl "${CURL_OPTS[@]}" "${PORTAL}/filereport?accession=${run}&result=read_run&fields=${FIELDS}&format=tsv" \
        | { if [[ -s "${META}.part" ]]; then tail -n +2; else cat; fi; } >> "${META}.part" \
        || echo "[warn] metadata request failed for ${run}" >&2
      sleep 1
    done
  fi
  mv "${META}.part" "$META"
  got=$(( $(wc -l < "$META") - 1 ))
  [[ $got -eq $N_RUNS ]] || echo "[warn] ENA metadata has ${got}/${N_RUNS} runs" >&2
fi

if [[ -s "$META" ]]; then
  fetch_xml "$(column_of "$META" run_accession)"              "${ENA_DIR}/runs.xml"
  fetch_xml "$(column_of "$META" experiment_accession)"       "${ENA_DIR}/experiments.xml"
  fetch_xml "$(column_of "$META" secondary_sample_accession)" "${ENA_DIR}/samples.xml"
elif [[ "$DRY_RUN" == "1" ]]; then
  for f in runs experiments samples; do fetch_xml "ERR..." "${ENA_DIR}/${f}.xml"; done
fi

# -----------------------------------------------------------------------------
# 2. FASTQs: throttled, staggered, resumable, MD5-verified from the cache.
# -----------------------------------------------------------------------------
md5_of() {
  if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'
  else md5 -q "$1"; fi   # macOS
}

meta_field() {  # look up column $2 for run $1 in the cached report
  [[ -s "$META" ]] || return 0
  awk -F'\t' -v r="$1" -v col="$2" 'NR==1{for(i=1;i<=NF;i++){if($i=="run_accession")a=i; if($i==col)m=i}; next}
                                    a && m && $a==r {print $m}' "$META"
}

file_size() { wc -c < "$1" | tr -d ' '; }

download_one() {
  local run="$1" name="$2" slot="$3"
  # ENA layout for 10-char accessions: vol1/fastq/<first 6>/00<last digit>/<run>/
  local url="${BASE}/${run:0:6}/00${run: -1}/${run}/${run}.fastq.gz"
  local out="${FQ_DIR}/${run}_${name}.fastq.gz"
  local done_flag="${out}.done"

  if [[ -f "$done_flag" ]]; then echo "[skip] ${run} already downloaded and verified"; return 0; fi
  if [[ "$DRY_RUN" == "1" ]]; then echo "[dry-run] ${url} -> ${out}"; return 0; fi

  # Stagger starts so parallel jobs don't open connections at the same moment.
  sleep $(( (slot % JOBS) * 5 + RANDOM % 5 ))

  # Use the expected size from ENA to decide what's on disk: a file that is
  # already complete (e.g. the script was stopped during the MD5 check) skips
  # the download; a file larger than expected is corrupt and starts over.
  local want_bytes have_bytes=0
  want_bytes=$(meta_field "$run" fastq_bytes)
  [[ -f "$out" ]] && have_bytes=$(file_size "$out")
  if [[ -n "$want_bytes" && "$have_bytes" -gt "$want_bytes" ]]; then
    echo "[warn] ${run}: local file larger than expected; starting over" >&2
    rm -f "$out"; have_bytes=0
  fi

  if [[ -n "$want_bytes" && "$have_bytes" -eq "$want_bytes" ]]; then
    echo "[have] ${run}: already complete, verifying"
  else
    local extra=() attempt
    [[ -n "$LIMIT_RATE" ]] && extra+=(--limit-rate "$LIMIT_RATE")
    for attempt in 1 2 3; do
      echo "[get ] ${run} (${name}) attempt ${attempt}"
      if curl "${CURL_OPTS[@]}" "${extra[@]}" -C - -o "$out" "$url"; then break; fi
      if [[ $attempt -eq 3 ]]; then echo "[FAIL] ${run}: download failed 3 times" >&2; return 1; fi
      echo "[wait] ${run}: pausing $(( attempt * 120 ))s before resuming" >&2
      sleep $(( attempt * 120 ))
    done
  fi

  echo "[done] ${run}: download complete; checksum validation follows"
}

validate_fastqs() {
  local failed=0 run name out done_flag expected actual

  while IFS=$'\t' read -r run name _; do
    out="${FQ_DIR}/${run}_${name}.fastq.gz"
    done_flag="${out}.done"

    if [[ ! -s "$out" ]]; then
      echo "[FAIL] ${run}: FASTQ is missing or empty"
      rm -f "$done_flag"
      failed=1
      continue
    fi

    if [[ "$SKIP_MD5" == "1" ]]; then
      echo "[skip] ${run}: checksum validation disabled"
      touch "$done_flag"
      continue
    fi

    expected=$(meta_field "$run" fastq_md5)
    if [[ -z "$expected" ]]; then
      echo "[FAIL] ${run}: no MD5 checksum in ENA metadata"
      rm -f "$done_flag"
      failed=1
      continue
    fi

    actual=$(md5_of "$out")
    if [[ "$expected" != "$actual" ]]; then
      echo "[FAIL] ${run}: MD5 mismatch (expected ${expected}, got ${actual}); deleting"
      rm -f "$out" "$done_flag"
      failed=1
      continue
    fi

    touch "$done_flag"
    echo "[ ok ] ${run}: MD5 verified"
  done < <(tail -n +2 "$SAMPLE_SHEET")

  return "$failed"
}

export -f download_one md5_of meta_field file_size
export FQ_DIR META BASE DRY_RUN SKIP_MD5 JOBS LIMIT_RATE
export CURL_OPTS_STR="${CURL_OPTS[*]}"

echo "Downloading ${N_RUNS} FASTQs to ${FQ_DIR}/ with ${JOBS} parallel jobs"
download_status=0
tail -n +2 "$SAMPLE_SHEET" | cut -f1,2 | awk -F'\t' '{print $1"\t"$2"\t"NR-1}' \
  | xargs -P "$JOBS" -L 1 bash -c 'read -ra CURL_OPTS <<< "$CURL_OPTS_STR"; download_one "$0" "$1" "$2"' \
  || download_status=$?

validation_status=0
if [[ "$DRY_RUN" != "1" ]] && ! validate_fastqs; then
  validation_status=1
fi

if [[ "$download_status" -ne 0 || "$validation_status" -ne 0 ]]; then
  echo "[FAIL] Download or checksum validation failed; re-run the script to resume."
  exit 1
fi

echo "All done. Sample sheet: ${SAMPLE_SHEET}  ENA metadata: ${ENA_DIR}/  Log: ${LOG_FILE}"
