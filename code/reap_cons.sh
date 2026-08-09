#!/usr/bin/env bash
# Background janitor: delete cons.* files as MFC writes them, so they never
# accumulate on a space-constrained disk. The analysis needs only prim.2, so
# conservative variables are pure waste here. Run this in a second terminal
# BEFORE launching the simulation; it polls until the run finishes.
#
#   bash reap_cons.sh /data/GitHub/MFC/3D_sphbubcollapse_choke/D
#
# Stop it with Ctrl-C after the run completes (or it auto-exits after 20 idle
# polls with no new files).
set -uo pipefail
DIR="${1:?usage: reap_cons.sh <case>/D}"
echo "reaping cons.* in $DIR every 15s (keeps prim.* only)"
idle=0
while true; do
    n=$(ls "$DIR"/cons.*.dat 2>/dev/null | wc -l)
    if [ "$n" -gt 0 ]; then
        rm -f "$DIR"/cons.*.dat
        echo "$(date +%T)  reaped $n cons files; $(df -h "$DIR" | awk 'NR==2{print $4" free"}')"
        idle=0
    else
        idle=$((idle+1))
        [ $((idle % 8)) -eq 0 ] && echo "$(date +%T)  idle ($idle) $(df -h "$DIR" | awk 'NR==2{print $4" free"}')"
    fi
    # exit if run appears done: no cons for many polls AND prim files exist
    if [ "$idle" -ge 40 ] && [ "$(ls "$DIR"/prim.*.dat 2>/dev/null | wc -l)" -gt 0 ]; then
        echo "no new cons for 10 min; assuming run finished. exiting."
        break
    fi
    sleep 15
done
