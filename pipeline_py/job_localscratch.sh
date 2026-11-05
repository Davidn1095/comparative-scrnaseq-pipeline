# Per-job temporary folder, sourced by pipeline_py/run.sh (from its snapshot), pipeline_py/genotype.sh and
# pipeline_py/cellranger.sh.
#   - JOBTMP = this job's folder under $PIPELINE_LOCAL_TMP (default: $TMPDIR, else /tmp; point it at node-local disk on
#     a cluster); TMPDIR points to it, so R and Python temporary files land inside it, and it is deleted on exit
#   - ls_run <cmd...>: runs a long command in its own process group and waits for it; bash defers a trapped signal
#     until the foreground command ends, but 'wait' is interrupted at once
#   - USR1 (e.g. sent before the time limit with sbatch --signal=B:USR1@300), TERM or INT: stops the running command's
#     process group, deletes JOBTMP and the folders in LS_EXTRA, and exits 1; EXIT deletes them too
#   - LS_ON_FAIL (empty by default): a command run once if the job stops on a signal or exits with an error
LS_TMP=${PIPELINE_LOCAL_TMP:-${TMPDIR:-/tmp}}
JOBTMP=$LS_TMP/$(id -un)_job_${SLURM_JOB_ID:-local$$}
LS_EXTRA=()
LS_KILL=()                           # background helpers (e.g. a disk monitor) stopped on exit
LS_ON_FAIL=""                        # optional command run once when the job ends with an error (argument: reason)
_ls_pid=""

ls_clean() { local d; for d in "${LS_KILL[@]}"; do kill "$d" 2>/dev/null || true; done; for d in "$JOBTMP" "${LS_EXTRA[@]}"; do [ -n "$d" ] && rm -rf --one-file-system "$d" || true; done; }

ls_stop() {
  echo "$(date '+%F %T') received $1: stopping the running command and deleting $JOBTMP ${LS_EXTRA[*]}"
  trap '' USR1 TERM INT
  if [ -n "$_ls_pid" ]; then
    kill -TERM -- -"$_ls_pid" 2>/dev/null || true
    local i; for i in $(seq 20); do kill -0 -- -"$_ls_pid" 2>/dev/null || break; sleep 1; done
    kill -KILL -- -"$_ls_pid" 2>/dev/null || true
  fi
  [ -n "$LS_ON_FAIL" ] && { $LS_ON_FAIL "stopped by $1" || true; }
  ls_clean; trap - EXIT; exit 1
}

ls_run() {
  setsid "$@" & _ls_pid=$!
  local rc=0; wait "$_ls_pid" || rc=$?
  _ls_pid=""; return $rc
}

mkdir -p "$JOBTMP"; export TMPDIR=$JOBTMP APPTAINERENV_TMPDIR=$JOBTMP
trap '_rc=$?; [ $_rc -ne 0 ] && [ -n "$LS_ON_FAIL" ] && { $LS_ON_FAIL "exit status $_rc" || true; }; ls_clean' EXIT
trap 'ls_stop USR1' USR1; trap 'ls_stop TERM' TERM; trap 'ls_stop INT' INT
