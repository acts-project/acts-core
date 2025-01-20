#!/usr/bin/env bash
set -euo pipefail

N=$1
shift
declare -a pids=()

# Start all processes in the background
for i in $(seq 1 "$N"); do
    # Replace `sleep 10` with the actual command you want to run.
    # For demonstration, we are using a command that sleeps for 10 seconds.
    # Make sure it runs in the background with '&'.
    "$@" &
    pids+=($!)
done

# Wait for all processes to finish, if any fails, kill them all
for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
        echo "Process $pid failed. Terminating all remaining processes..."
        # Kill all started processes
        kill "${pids[@]}" 2>/dev/null || true
        exit 1
    fi
done

echo "All processes completed successfully."
exit 0
