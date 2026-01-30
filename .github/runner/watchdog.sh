#!/bin/bash
# GitHub Actions Runner Watchdog
# Monitors runner health and restarts on socket timeout issues
# Usage: ./watchdog.sh [--daemon]

# Configuration - set RUNNER_DIR env var or edit this default
RUNNER_DIR="${RUNNER_DIR:-$HOME/actions-runner}"
LOG_FILE="$RUNNER_DIR/_diag/watchdog.log"
PID_FILE="$RUNNER_DIR/.watchdog.pid"
CHECK_INTERVAL=60  # Check every 60 seconds
STALL_THRESHOLD=300  # Consider stalled if no log activity for 5 minutes during a job

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

is_runner_running() {
    pgrep -f "Runner.Listener" > /dev/null 2>&1
}

get_runner_pid() {
    pgrep -f "Runner.Listener" | head -1
}

get_latest_log() {
    ls -t "$RUNNER_DIR/_diag/Runner_"*.log 2>/dev/null | head -1
}

check_for_socket_errors() {
    local log_file=$(get_latest_log)
    if [[ -f "$log_file" ]]; then
        # Check last 100 lines for socket timeout errors in the last 5 minutes
        local recent_errors=$(tail -100 "$log_file" | grep -c "Socket Error: TimedOut\|maximum number of attempts has been reached" 2>/dev/null || echo "0")
        echo "$recent_errors"
    else
        echo "0"
    fi
}

check_if_stalled() {
    local log_file=$(get_latest_log)
    if [[ -f "$log_file" ]]; then
        local last_modified=$(stat -f %m "$log_file" 2>/dev/null || stat -c %Y "$log_file" 2>/dev/null)
        local now=$(date +%s)
        local age=$((now - last_modified))

        # If log hasn't been updated in STALL_THRESHOLD seconds, might be stalled
        if [[ $age -gt $STALL_THRESHOLD ]]; then
            echo "1"
        else
            echo "0"
        fi
    else
        echo "0"
    fi
}

restart_runner() {
    log "Restarting runner..."

    # Kill existing processes
    pkill -f "Runner.Listener" 2>/dev/null
    pkill -f "Runner.Worker" 2>/dev/null
    sleep 2

    # Double-check they're dead
    pkill -9 -f "Runner.Listener" 2>/dev/null
    pkill -9 -f "Runner.Worker" 2>/dev/null
    sleep 1

    # Start runner
    cd "$RUNNER_DIR"
    nohup ./run.sh >> "$RUNNER_DIR/_diag/launchd-stdout.log" 2>&1 &

    log "Runner restarted with PID $(get_runner_pid)"
}

check_github_connectivity() {
    # Quick connectivity test to GitHub
    if curl -s --max-time 10 "https://api.github.com/zen" > /dev/null 2>&1; then
        echo "1"
    else
        echo "0"
    fi
}

watchdog_loop() {
    log "Watchdog started (PID $$)"
    echo $$ > "$PID_FILE"

    local consecutive_errors=0
    local max_consecutive_errors=3

    while true; do
        sleep $CHECK_INTERVAL

        # Check if runner is running at all
        if ! is_runner_running; then
            log "Runner not running, starting it..."
            restart_runner
            consecutive_errors=0
            continue
        fi

        # Check for socket timeout errors
        local socket_errors=$(check_for_socket_errors)
        if [[ "$socket_errors" -gt 5 ]]; then
            log "Detected $socket_errors socket timeout errors in recent logs"
            ((consecutive_errors++))
        else
            consecutive_errors=0
        fi

        # Check if stalled
        local is_stalled=$(check_if_stalled)
        if [[ "$is_stalled" == "1" ]]; then
            log "Runner appears stalled (no log activity for ${STALL_THRESHOLD}s)"
            ((consecutive_errors++))
        fi

        # Check GitHub connectivity
        local github_ok=$(check_github_connectivity)
        if [[ "$github_ok" != "1" ]]; then
            log "GitHub connectivity issue detected"
            # Don't increment errors here - just log it
        fi

        # Restart if too many consecutive errors
        if [[ $consecutive_errors -ge $max_consecutive_errors ]]; then
            log "Too many consecutive errors ($consecutive_errors), forcing restart"
            restart_runner
            consecutive_errors=0
        fi

        # Periodic health log (every 10 checks = ~10 minutes)
        if [[ $((RANDOM % 10)) -eq 0 ]]; then
            local pid=$(get_runner_pid)
            log "Health check: Runner PID=$pid, Errors=$consecutive_errors, GitHub=${github_ok}"
        fi
    done
}

stop_watchdog() {
    if [[ -f "$PID_FILE" ]]; then
        local wpid=$(cat "$PID_FILE")
        if kill -0 "$wpid" 2>/dev/null; then
            log "Stopping watchdog (PID $wpid)"
            kill "$wpid"
            rm -f "$PID_FILE"
        fi
    fi
}

case "${1:-}" in
    --daemon|-d)
        # Run in background
        nohup "$0" >> "$LOG_FILE" 2>&1 &
        echo "Watchdog started in background (PID $!)"
        ;;
    --stop)
        stop_watchdog
        ;;
    --status)
        if [[ -f "$PID_FILE" ]] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
            echo "Watchdog running (PID $(cat "$PID_FILE"))"
            echo "Runner: $(is_runner_running && echo "running (PID $(get_runner_pid))" || echo "not running")"
        else
            echo "Watchdog not running"
        fi
        ;;
    --restart-now)
        log "Manual restart requested"
        restart_runner
        ;;
    *)
        # Run in foreground
        watchdog_loop
        ;;
esac
