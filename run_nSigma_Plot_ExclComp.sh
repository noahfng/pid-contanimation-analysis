#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${DEBUG:-0}" -eq 1 ]]; then
  PS4='+ ${BASH_SOURCE##*/}:${LINENO}:${FUNCNAME[0]:-MAIN}() '
  set -x
fi

ulimit -c unlimited || true

on_error() {
  local ec=$?
  echo "ERROR: Exit ${ec} at ${BASH_SOURCE[1]}:${BASH_LINENO[0]}" >&2
  echo "Last command: ${BASH_COMMAND}" >&2
  echo "Context: run=${CURRENT_RUN:-?} tag=${RUN_TAG:-?} P=${PSTART:-?}:${PEND:-?} X=${XMIN:-?}:${XMAX:-?} NGAUSS=${MANUAL_NGAUSS:-?}" >&2
}
trap on_error ERR

PRANGES=("0.35:0.45" "0.45:0.55" "0.55:0.65" "0.65:0.75" "0.75:0.85" "0.85:0.95" "0.95:1.05" "1.05:1.15"    )
XRANGES=("-12:16" "-12:8" "-12:5" "-10:15" "-10:6" "-10:7" "-10:5" "-10:4")
MANUAL_SETS=(
  "4|-6,-5,0,4|1,1,1,4|1e6,1e6,1e5,1e4|mu,pi,e,K"
  "4|-6,-5,0,2|1,1,1,2|1e6,1e6,1e5,1e5|mu,pi,e,K"
  "4|-7,-5,-2,0|1,1,1,1|1e6,1e6,1e5,1e4|mu,pi,e,K"
  "5|-6,-5,-3,0,8|1,1,1,1,4|1e6,1e6,1e5,1e5,1e4|mu,pi,e,K,p"
  "5|-6,-5,-2,0,4|1,1,1,1,2|1e6,1e6,1e4,1e4,1e4|mu,pi,e,K,p"
  "4|-6,-5,0,2|1,1,1,1|1e6,1e6,5e4,5e4|mu,pi,e,K+p"
  "4|-6,-5,-1,1|1,1,1,1|3e5,3e5,1e4,1e4|mu,pi,e,K+p"
  "4|-6,-5,-1,0|1,1,1,1|3e5,3e5,1e4,1e4|mu,pi,e,K+p"
)

FITKAON=("on" "on" "on" "on" "on" "on" "on" "on")  
FITPROTON=("off" "off" "off" "on" "on" "on" "on" "on")

SRC="nSigma_Plot_ExclComp.Cpp"
FUNC="nSigma_Plot_ExclComp"

BASE_OUTDIR="overnight_runs_ExlComp"
MACRO_DIR=$(cd "$(dirname "$SRC")" && pwd)

SESSION_TS=$(date +"%Y-%m-%dT%H:%M:%S%z")   
LOG_DIR="${BASE_OUTDIR}/logs"
mkdir -p "$LOG_DIR"
MASTER_LOG="${LOG_DIR}/session-${SESSION_TS//[:+]/_}.log"
exec > >(tee -a "$MASTER_LOG") 2>&1


command -v root >/dev/null 2>&1 || { echo "ERROR: 'root' not in PATH"; exit 1; }
[[ -f "$SRC" ]] || { echo "ERROR: Macro not found: $SRC"; exit 1; }

LEN=${#PRANGES[@]}
[[ $LEN -eq ${#XRANGES[@]} && $LEN -eq ${#MANUAL_SETS[@]} ]] \
  || { echo "ERROR: PRANGES/XRANGES/MANUAL_SETS have to be same length."; exit 1; }

broadcast_or_check() {
  local name="$1"; shift
  local -n arr="$name"
  if [[ ${#arr[@]} -eq 1 ]]; then
    local tmp=()
    for ((i=0;i<LEN;i++)); do tmp+=("${arr[0]}"); done
    arr=("${tmp[@]}")
  elif [[ ${#arr[@]} -ne $LEN ]]; then
    echo "ERROR: ${name} must be length 1 (Broadcast) or ${LEN}."
    exit 1
  fi
}
broadcast_or_check FITKAON
broadcast_or_check FITPROTON

mkdir -p "$BASE_OUTDIR"

export LC_ALL=C
export LC_NUMERIC=C
export ROOT_HIST=0
export INTERACTIVE_MANUAL_PEAKS=0

for ((i=0; i<LEN; i++)); do
  CURRENT_RUN="$((i+1))"
  IFS=: read -r PSTART PEND <<<"${PRANGES[i]}"
  IFS=: read -r XMIN XMAX   <<<"${XRANGES[i]}"
  IFS='|' read -r NG MEANS SIGMAS AMPS LABELS <<<"${MANUAL_SETS[i]}"

  RUN_TAG=$(printf "run%02d-%s<p<%s-%s<x<%s-KaExcl%s-PrExcl%s" \
    "$((i+1))" "$PSTART" "$PEND" "$XMIN" "$XMAX" "${FITKAON[i]}" "${FITPROTON[i]}" \
    | tr ':' '_' | tr -d '+')
  OUTDIR="${BASE_OUTDIR}/${RUN_TAG}"
  mkdir -p "$OUTDIR"

  export XMIN XMAX PSTART PEND
  export FITKAONEXCLCOMP="${FITKAON[i]}"       
  export FITPROTONEXCLCOMP="${FITPROTON[i]}"
  export MANUAL_NGAUSS="$NG"
  export MANUAL_MEANS="$MEANS"
  export MANUAL_SIGMAS="$SIGMAS"
  export MANUAL_AMPS="$AMPS"
  export MANUAL_LABELS="$LABELS"

echo "=== [$((i+1))/$LEN] ${RUN_TAG} ==="
(
  root nSigma_Plot_ExclComp.Cpp -l -q -b
)

find "$MACRO_DIR" -maxdepth 1 -type f \
\( -name '*.pdf' -o -name 'nSigmaSummary_*.txt' \) \
  -exec mv -f {} "${OUTDIR}/" \;
done
echo "Finished, results in: ${BASE_OUTDIR}"