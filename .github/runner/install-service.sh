#!/bin/bash
# Install GitHub Actions Runner as a managed macOS service
# This script sets up:
# 1. launchd service for auto-restart
# 2. Watchdog for health monitoring
# 3. Network keepalive settings

set -e

# Configuration - set these env vars or edit defaults
RUNNER_DIR="${RUNNER_DIR:-$HOME/actions-runner}"
RUNNER_NAME="${RUNNER_NAME:-wasp2}"
PLIST_NAME="com.github.actions.runner.${RUNNER_NAME}"
PLIST_SRC="$(dirname "$0")/com.github.actions.runner.plist"
PLIST_DST="$HOME/Library/LaunchAgents/$PLIST_NAME.plist"
WATCHDOG_PLIST_NAME="com.github.actions.runner.${RUNNER_NAME}.watchdog"

echo "=== GitHub Actions Runner Service Installer ==="
echo ""

# Step 1: Stop existing runner processes
echo "[1/6] Stopping existing runner processes..."
pkill -f "Runner.Listener" 2>/dev/null || true
pkill -f "Runner.Worker" 2>/dev/null || true
launchctl unload "$PLIST_DST" 2>/dev/null || true
sleep 2

# Step 2: Make scripts executable
echo "[2/6] Setting permissions..."
chmod +x "$RUNNER_DIR/run.sh"
chmod +x "$RUNNER_DIR/watchdog.sh"
chmod +x "$RUNNER_DIR/config.sh"

# Step 3: Create LaunchAgents directory if needed
echo "[3/6] Setting up LaunchAgents..."
mkdir -p "$HOME/Library/LaunchAgents"

# Step 4: Generate and install the runner plist from template
echo "[4/6] Installing launchd service..."
sed -e "s|RUNNER_NAME|${RUNNER_NAME}|g" \
    -e "s|RUNNER_DIR|${RUNNER_DIR}|g" \
    -e "s|HOME_DIR|${HOME}|g" \
    "$PLIST_SRC" > "$PLIST_DST"

# Copy watchdog script to runner directory if not already there
SCRIPT_DIR="$(dirname "$0")"
if [[ -f "$SCRIPT_DIR/watchdog.sh" ]] && [[ ! -f "$RUNNER_DIR/watchdog.sh" ]]; then
    cp "$SCRIPT_DIR/watchdog.sh" "$RUNNER_DIR/watchdog.sh"
    chmod +x "$RUNNER_DIR/watchdog.sh"
fi

# Step 5: Create watchdog launchd plist
echo "[5/6] Installing watchdog service..."
cat > "$HOME/Library/LaunchAgents/$WATCHDOG_PLIST_NAME.plist" << WATCHDOG_PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>Label</key>
    <string>$WATCHDOG_PLIST_NAME</string>
    <key>ProgramArguments</key>
    <array>
        <string>$RUNNER_DIR/watchdog.sh</string>
    </array>
    <key>WorkingDirectory</key>
    <string>$RUNNER_DIR</string>
    <key>EnvironmentVariables</key>
    <dict>
        <key>RUNNER_DIR</key>
        <string>$RUNNER_DIR</string>
    </dict>
    <key>RunAtLoad</key>
    <true/>
    <key>KeepAlive</key>
    <true/>
    <key>ThrottleInterval</key>
    <integer>30</integer>
    <key>StandardOutPath</key>
    <string>$RUNNER_DIR/_diag/watchdog-stdout.log</string>
    <key>StandardErrorPath</key>
    <string>$RUNNER_DIR/_diag/watchdog-stderr.log</string>
</dict>
</plist>
WATCHDOG_PLIST

# Step 6: Load and start services
echo "[6/6] Starting services..."
launchctl load "$PLIST_DST"
launchctl load "$HOME/Library/LaunchAgents/$WATCHDOG_PLIST_NAME.plist"

sleep 3

echo ""
echo "=== Installation Complete ==="
echo ""
echo "Services installed:"
echo "  - Runner:   $PLIST_NAME"
echo "  - Watchdog: $WATCHDOG_PLIST_NAME"
echo ""
echo "Useful commands:"
echo "  Check status:     launchctl list | grep actions"
echo "  View runner logs: tail -f $RUNNER_DIR/_diag/Runner_*.log"
echo "  View watchdog:    tail -f $RUNNER_DIR/_diag/watchdog.log"
echo "  Stop runner:      launchctl unload ~/Library/LaunchAgents/$PLIST_NAME.plist"
echo "  Start runner:     launchctl load ~/Library/LaunchAgents/$PLIST_NAME.plist"
echo "  Restart runner:   $RUNNER_DIR/watchdog.sh --restart-now"
echo ""

# Verify
echo "Current status:"
launchctl list | grep -E "actions|runner" || echo "  (services starting...)"
echo ""
pgrep -fl "Runner.Listener" || echo "Runner process starting..."
