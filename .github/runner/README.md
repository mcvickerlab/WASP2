# Self-Hosted GitHub Actions Runner Setup

Scripts to run a self-hosted GitHub Actions runner as a managed macOS service with automatic restart and health monitoring.

## Problem Solved

Self-hosted runners often get "stuck" in a busy state due to:
- Socket timeout errors to GitHub's broker service
- Network connectivity issues causing silent connection drops
- macOS power management interrupting long-poll connections

This setup provides:
1. **launchd service** - Auto-restarts runner on crash/exit
2. **Watchdog** - Monitors for stuck states and socket errors, force-restarts if needed

## Quick Install

```bash
# Set your runner directory and name
export RUNNER_DIR="$HOME/actions-runner"
export RUNNER_NAME="wasp2"

# Run installer
.github/runner/install-service.sh
```

## Manual Installation

1. Download and configure the GitHub Actions runner:
   ```bash
   mkdir ~/actions-runner && cd ~/actions-runner
   # Download from https://github.com/actions/runner/releases
   ./config.sh --url https://github.com/YOUR/REPO --token YOUR_TOKEN
   ```

2. Copy scripts to runner directory:
   ```bash
   cp .github/runner/watchdog.sh ~/actions-runner/
   chmod +x ~/actions-runner/watchdog.sh
   ```

3. Install services:
   ```bash
   RUNNER_DIR=~/actions-runner RUNNER_NAME=myrunner .github/runner/install-service.sh
   ```

## Commands

```bash
# Check service status
launchctl list | grep actions

# View runner logs
tail -f ~/actions-runner/_diag/Runner_*.log

# View watchdog logs
tail -f ~/actions-runner/_diag/watchdog.log

# Manual restart
~/actions-runner/watchdog.sh --restart-now

# Watchdog status
~/actions-runner/watchdog.sh --status

# Stop services
launchctl unload ~/Library/LaunchAgents/com.github.actions.runner.*.plist

# Start services
launchctl load ~/Library/LaunchAgents/com.github.actions.runner.wasp2.plist
launchctl load ~/Library/LaunchAgents/com.github.actions.runner.wasp2.watchdog.plist
```

## How It Works

### Watchdog Detection

The watchdog checks every 60 seconds for:
- Runner process not running → immediate restart
- 5+ socket timeout errors in logs → increment error counter
- No log activity for 5 minutes → increment error counter
- 3 consecutive error cycles → force restart

### launchd Service

The launchd plist provides:
- Auto-start on boot
- Auto-restart on crash (10s throttle)
- Higher process priority to avoid being killed under memory pressure
- Proper environment variables for build tools

## Troubleshooting

**Runner shows "busy" but nothing running:**
```bash
# Force restart
~/actions-runner/watchdog.sh --restart-now
```

**Services not starting:**
```bash
# Check for errors
launchctl list | grep actions
cat ~/actions-runner/_diag/launchd-stderr.log
```

**Reinstall services:**
```bash
launchctl unload ~/Library/LaunchAgents/com.github.actions.runner.*.plist
rm ~/Library/LaunchAgents/com.github.actions.runner.*.plist
RUNNER_DIR=~/actions-runner RUNNER_NAME=wasp2 .github/runner/install-service.sh
```
