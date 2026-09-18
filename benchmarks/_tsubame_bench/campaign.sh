#!/bin/sh
#=======================================================================
# TSUBAME4.0 scaling campaign driver, 202607.
#
# Counterpart of bench-202607/scripts/campaign.sh on Miyabi.  It does nothing
# clever: it walks the configuration table and calls each tree's own
# job.sh-tsubame, which is what knows how to re-submit itself under AGE.
#
#   ./campaign.sh list                 # print every qsub that would happen
#   ./campaign.sh fig1                 # slab + pencil + P3DFFT, strong scaling
#   ./campaign.sh fig2                 # gallope sep/uni + calliope
#   ./campaign.sh fig3                 # slab, weak scaling
#   ./campaign.sh placement            # the placement control experiment
#   ./campaign.sh all
#   TRIALS="1" ./campaign.sh fig2      # one trial first, then decide
#   GPU_SPAN=intra ./campaign.sh fig1  # only what the inter-node gate does not block
#   GPU_SPAN=inter ./campaign.sh fig1  # the rest, once the gate has passed
#
# TRIALS defaults to "1 2 3": three INDEPENDENT submissions per configuration,
# which is what gives the median and the min-max whiskers.  Trials are never
# looped over inside a job.
#
# Resource choice per GPU count.  A node holds 4 H100, so 1 and 2 GPUs are
# taken from node_q / node_h rather than occupying a whole node_f and paying
# for three idle GPUs:
#
#    GPUs :  1       2       3   4  |  8      12      16      32
#    alloc: node_q  node_h   node_f |  node_f x 2, 3, 4, 8
#    nodes:  1       1        1     |  2       3       4       8
#           +---- NVLink only ----+ +----- + InfiniBand -----+
#
# The 3- and 12-GPU points are why Fig. 1 uses 768^3 / 1536^3 here instead of
# the 1024^3 / 2048^3 of Miyabi: cufftmp_r2c.f90:142 refuses mod(ny,nproc) > 0,
# and 768 = 2^8*3, 1536 = 2^9*3 divide by every count in the sweep.
#=======================================================================
set -u

WORK=${WORK:-${HOME}/work}
SLAB=${WORK}/p3dfft_vs_cufftmp/cufftmp-slab
PENCIL=${WORK}/p3dfft_vs_cufftmp/cufftmp-pencil
P3DFFT=${WORK}/p3dfft_vs_cufftmp/p3dfft
# The application trees are named *_performance on TSUBAME, not *_miyabi; the
# unified-memory build of Gallope has no tree of its own here, so GALLOPE_DIRS
# holds one entry and fig2 loops over whatever is present rather than assuming
# both.  Point GALLOPE_DIRS at two directories the day a unified build exists.
GALLOPE=${WORK}/gallope_performance
# "-" and not ":-": GALLOPE_DIRS="" must mean "no Gallope runs at all" (the way
# to submit the CPU half of fig2 on its own), and ":-" would silently treat the
# empty value as unset and submit the whole GPU side anyway.  It did, once.
GALLOPE_DIRS=${GALLOPE_DIRS-"${GALLOPE}"}
CALLIOPE=${WORK}/calliope_performance
# Same knob on the CPU side, and for the mirror-image reason: when one half of
# fig2 has to be resubmitted (a stale binary, a bad build) the other half is
# usually still queued, and re-running the whole target would submit it twice
# under the same trial index.  CALLIOPE_DIRS="" submits the GPU half alone.
CALLIOPE_DIRS=${CALLIOPE_DIRS-"${CALLIOPE}"}

TRIALS=${TRIALS:-"1 2 3"}
RESULTS=${RESULTS:-results-202607/tsubame}
DRY=0

# Grids.  Fig. 1 and Fig. 3 keep the Miyabi protocol; only Fig. 1's box changed,
# for the divisibility reason in the header.  Fig. 2 uses 384^3 = 2^7*3, which
# divides by every GPU count in the sweep AND fits one H100 (~23 GB), so the
# 1-GPU end of the application curve exists here -- at 1024^3 it would not.
FIG1_GRIDS=${FIG1_GRIDS:-768:1536}
FIG3_GRIDS=${FIG3_GRIDS:-512:1024}
FIG2_GRIDS=${FIG2_GRIDS:-384}

# GPU_SPAN selects which GPU submissions go out.  cuFFTMp reaches off-node
# through NVSHMEM, and until a 2-node job is seen to complete, those runs would
# either fail or -- worse -- report numbers from a transport that quietly fell
# back.  The single-node GPU points and the whole CPU side never touch NVSHMEM,
# so they need not wait behind a gate that does not apply to them.
#
#   intra  single-node GPU + all CPU   -- safe to run before the gate
#   inter  multi-node GPU only         -- exactly the complement, so running
#                                         "intra" then "inter" covers the set
#                                         once each and never re-runs a point
#   all    everything                  [default]
#
# The two halves are disjoint on purpose: re-running a whole figure after the
# gate would resubmit the single-node points under the same trial index and
# overwrite results that are already good.
GPU_SPAN=${GPU_SPAN:-all}

sub () {   # $1=dir  $2=script  rest=VAR=VAL ...
    dir=$1; script=$2; shift 2
    if [ "${GPU_SPAN}" != all ]; then
        is_gpu=0; nodes=1
        for kv in "$@"; do
            case ${kv} in
                GPUS=*)  is_gpu=1 ;;
                NODES=*) nodes=${kv#NODES=} ;;
            esac
        done
        if [ ${is_gpu} -eq 1 ] && [ "${nodes}" -gt 1 ] && [ "${GPU_SPAN}" = intra ]; then
            echo "  (held for the inter-node gate) ${dir} $*"
            return 0
        fi
        # "inter" is GPU-only: the CPU curves come out with the intra pass.
        if [ "${GPU_SPAN}" = inter ] && { [ ${is_gpu} -eq 0 ] || [ "${nodes}" -le 1 ]; }; then
            return 0
        fi
    fi
    if [ ${DRY} -eq 1 ]; then
        echo "  (cd ${dir} && $* RESULTS=${RESULTS} ${script})"
    else
        ( cd "${dir}" && env "$@" RESULTS="${RESULTS}" "${script}" )
    fi
}

# ---------------------------------------------------------------- Fig. 1
fig1 () {
for t in ${TRIALS}; do
    # slab: every GPU count
    sub "${SLAB}" ./job.sh-tsubame RESOURCE=node_q NODES=1 GPUS=1      GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=1:00:00
    sub "${SLAB}" ./job.sh-tsubame RESOURCE=node_h NODES=1 GPUS=2      GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=1:00:00
    sub "${SLAB}" ./job.sh-tsubame RESOURCE=node_f NODES=1 GPUS=3:4    GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=1:00:00
    sub "${SLAB}" ./job.sh-tsubame RESOURCE=node_f NODES=2 GPUS=8      GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=0:30:00
    sub "${SLAB}" ./job.sh-tsubame RESOURCE=node_f NODES=3 GPUS=12     GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=0:30:00
    sub "${SLAB}" ./job.sh-tsubame RESOURCE=node_f NODES=4 GPUS=16     GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=0:30:00
    sub "${SLAB}" ./job.sh-tsubame RESOURCE=node_f NODES=8 GPUS=32     GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=0:30:00
    # pencil: perfect squares only (nranks1d = sqrt(nproc) must be exact)
    sub "${PENCIL}" ./job.sh-tsubame RESOURCE=node_q NODES=1 GPUS=1    GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=1:00:00
    sub "${PENCIL}" ./job.sh-tsubame RESOURCE=node_f NODES=1 GPUS=4    GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=1:00:00
    sub "${PENCIL}" ./job.sh-tsubame RESOURCE=node_f NODES=4 GPUS=16   GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=0:30:00
    # P3DFFT: the CPU reference, on cpu_160 so no H100 is charged for idling
    for n in 1 2 3 4 8; do
        sub "${P3DFFT}" ./job.sh-tsubame RESOURCE=cpu_160 NODES=$n     GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=1:00:00
    done
done
}

# ---------------------------------------------------------------- Fig. 3
# Weak scaling holds the work per GPU fixed and doubles the box x -> y -> z,
# so only powers of two are meaningful; 3 and 12 are absent by construction,
# not by omission.
fig3 () {
for t in ${TRIALS}; do
    sub "${SLAB}" ./job.sh-tsubame MODE=weak RESOURCE=node_q NODES=1 GPUS=1  GRIDS=${FIG3_GRIDS} TRIAL=$t H_RT=1:00:00
    sub "${SLAB}" ./job.sh-tsubame MODE=weak RESOURCE=node_h NODES=1 GPUS=2  GRIDS=${FIG3_GRIDS} TRIAL=$t H_RT=1:00:00
    sub "${SLAB}" ./job.sh-tsubame MODE=weak RESOURCE=node_f NODES=1 GPUS=4  GRIDS=${FIG3_GRIDS} TRIAL=$t H_RT=1:00:00
    sub "${SLAB}" ./job.sh-tsubame MODE=weak RESOURCE=node_f NODES=2 GPUS=8  GRIDS=${FIG3_GRIDS} TRIAL=$t H_RT=0:30:00
    sub "${SLAB}" ./job.sh-tsubame MODE=weak RESOURCE=node_f NODES=4 GPUS=16 GRIDS=${FIG3_GRIDS} TRIAL=$t H_RT=0:30:00
    sub "${SLAB}" ./job.sh-tsubame MODE=weak RESOURCE=node_f NODES=8 GPUS=32 GRIDS=${FIG3_GRIDS} TRIAL=$t H_RT=0:30:00
done
}

# ---------------------------------------------------------------- Fig. 2
# One directory per Gallope build.  On Miyabi Fig. 2 carries a unified-memory
# curve as well, but -gpu=mem:unified on x86+H100 migrates over PCIe and is not
# the hardware-coherent GH200 unified memory that curve is about, so no such
# build was made here.  GALLOPE_DIRS is the knob if one ever is.
fig2 () {
for t in ${TRIALS}; do
    for d in ${GALLOPE_DIRS}; do
        sub "${d}" ./batch/job.sh-tsubame RESOURCE=node_q NODES=1 GPUS=1   GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=2:00:00
        sub "${d}" ./batch/job.sh-tsubame RESOURCE=node_h NODES=1 GPUS=2   GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=1:00:00
        sub "${d}" ./batch/job.sh-tsubame RESOURCE=node_f NODES=1 GPUS=3:4 GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=1:30:00
        sub "${d}" ./batch/job.sh-tsubame RESOURCE=node_f NODES=2 GPUS=8   GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=0:45:00
        sub "${d}" ./batch/job.sh-tsubame RESOURCE=node_f NODES=3 GPUS=12  GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=0:45:00
        sub "${d}" ./batch/job.sh-tsubame RESOURCE=node_f NODES=4 GPUS=16  GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=0:45:00
        sub "${d}" ./batch/job.sh-tsubame RESOURCE=node_f NODES=8 GPUS=32  GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=0:45:00
    done
    for d in ${CALLIOPE_DIRS}; do
        for n in 1 2 3 4 8; do
            sub "${d}" ./batch/job.sh-tsubame RESOURCE=cpu_160 NODES=$n GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=2:00:00
        done
    done
done
}

# ---------------------------------------------------------------- placement
# The decisive control.  Total GPU count and work per GPU are held fixed and
# only the number of nodes they are spread over changes, so nothing varies
# except whether the all-to-all crosses InfiniBand.  If the knee in the strong
# scaling curve is communication, these points reproduce it at constant GPU
# count; if it is anything else, they will not.
#
# The 4-GPU/1-node and 8-GPU/2-node cases are already in fig1/fig2, so they are
# not repeated here.
placement () {
for t in ${TRIALS}; do
    sub "${SLAB}" ./job.sh-tsubame RESOURCE=node_f NODES=2 GPUS=4 PPN=2 GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=0:30:00
    sub "${SLAB}" ./job.sh-tsubame RESOURCE=node_f NODES=4 GPUS=4 PPN=1 GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=0:30:00
    sub "${SLAB}" ./job.sh-tsubame RESOURCE=node_f NODES=4 GPUS=8 PPN=2 GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=0:30:00
    sub "${SLAB}" ./job.sh-tsubame RESOURCE=node_f NODES=8 GPUS=8 PPN=1 GRIDS=${FIG1_GRIDS} TRIAL=$t H_RT=0:30:00
    sub "${GALLOPE}" ./batch/job.sh-tsubame RESOURCE=node_f NODES=2 GPUS=4 PPN=2 GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=1:00:00
    sub "${GALLOPE}" ./batch/job.sh-tsubame RESOURCE=node_f NODES=4 GPUS=4 PPN=1 GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=1:00:00
    sub "${GALLOPE}" ./batch/job.sh-tsubame RESOURCE=node_f NODES=4 GPUS=8 PPN=2 GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=0:45:00
    sub "${GALLOPE}" ./batch/job.sh-tsubame RESOURCE=node_f NODES=8 GPUS=8 PPN=1 GRIDS=${FIG2_GRIDS} TRIAL=$t H_RT=0:45:00
done
}

case "${1:-list}" in
    list)      DRY=1; echo "--- fig1"; fig1; echo "--- fig3"; fig3
               echo "--- fig2"; fig2; echo "--- placement"; placement ;;
    fig1)      fig1 ;;
    fig2)      fig2 ;;
    fig3)      fig3 ;;
    placement) placement ;;
    all)       fig1; fig3; fig2; placement ;;
    *) echo "usage: $0 {list|fig1|fig2|fig3|placement|all}" >&2; exit 1 ;;
esac
