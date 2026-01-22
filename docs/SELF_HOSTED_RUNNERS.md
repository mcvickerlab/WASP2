# Self-Hosted Runner Configuration for WASP2

This guide explains how to configure self-hosted GitHub Actions runners for WASP2, particularly useful for HPC environments and GPU-enabled processing of large BAM files.

## Why Self-Hosted Runners?

WASP2 processes genomics data that can be extremely large (BAM files often exceed 100GB). Self-hosted runners provide:

- **Access to large storage**: HPC clusters typically have terabytes of scratch space
- **GPU acceleration**: For compute-intensive variant calling operations
- **Network locality**: Direct access to shared genomics data without transfer
- **Custom environments**: Pre-configured bioinformatics tools and reference genomes

## Setting Up a Self-Hosted Runner

### Prerequisites

1. A Linux machine with:
   - At least 32GB RAM (64GB+ recommended for large BAM files)
   - 500GB+ scratch storage
   - Network access to GitHub
   - (Optional) NVIDIA GPU with CUDA 11.0+

2. Required software:
   - Python 3.10+
   - Rust toolchain
   - bcftools, samtools, bedtools
   - (Optional) CUDA toolkit

### Installation Steps

#### 1. Create a dedicated user

```bash
sudo useradd -m -s /bin/bash github-runner
sudo usermod -aG docker github-runner  # If using Docker
```

#### 2. Download and configure the runner

```bash
# As the github-runner user
cd ~
mkdir actions-runner && cd actions-runner

# Download latest runner (check releases for current version)
curl -o actions-runner-linux-x64-2.311.0.tar.gz -L \
  https://github.com/actions/runner/releases/download/v2.311.0/actions-runner-linux-x64-2.311.0.tar.gz

tar xzf ./actions-runner-linux-x64-2.311.0.tar.gz
```

#### 3. Register the runner

Navigate to your repository settings:
`Settings` → `Actions` → `Runners` → `New self-hosted runner`

```bash
./config.sh --url https://github.com/mcvickerlab/WASP2 \
  --token YOUR_TOKEN \
  --name "hpc-runner-01" \
  --labels "self-hosted,linux,x64,hpc,gpu" \
  --work "_work"
```

#### 4. Install as a service

```bash
sudo ./svc.sh install
sudo ./svc.sh start
```

## Runner Labels

Configure runners with appropriate labels for job targeting:

| Label | Description |
|-------|-------------|
| `hpc` | High-memory HPC nodes |
| `gpu` | GPU-enabled nodes |
| `large-storage` | Nodes with >1TB scratch |
| `bioinformatics` | Pre-configured with bio tools |

## Workflow Configuration for Self-Hosted Runners

### Example: Large BAM Processing Job

```yaml
jobs:
  process-large-bam:
    runs-on: [self-hosted, linux, hpc, large-storage]
    steps:
      - uses: actions/checkout@v4

      - name: Process BAM files
        run: |
          # Access shared data on HPC
          wasp2-map make-reads /shared/data/sample.bam /shared/ref/variants.vcf \
            --threads ${{ runner.cpus }} \
            --out_dir ${{ runner.temp }}/output
```

### Example: GPU-Accelerated Job

```yaml
jobs:
  gpu-analysis:
    runs-on: [self-hosted, linux, gpu]
    steps:
      - uses: actions/checkout@v4

      - name: Run GPU-accelerated analysis
        run: |
          nvidia-smi  # Verify GPU availability
          # Your GPU-enabled analysis here
```

## Test Data Caching

For CI/CD with test data, configure local caching:

### 1. Create a persistent cache directory

```bash
sudo mkdir -p /data/github-runner-cache
sudo chown github-runner:github-runner /data/github-runner-cache
```

### 2. Configure the runner environment

Add to `/etc/environment` or runner's `.bashrc`:

```bash
export WASP2_TEST_DATA_CACHE=/data/github-runner-cache/test-data
export WASP2_REFERENCE_CACHE=/data/github-runner-cache/references
```

### 3. Use caching in workflows

```yaml
- name: Cache test data
  uses: actions/cache@v4
  with:
    path: |
      ${{ env.WASP2_TEST_DATA_CACHE }}
    key: test-data-${{ hashFiles('tests/data/checksums.txt') }}
    restore-keys: |
      test-data-

- name: Download test data if not cached
  run: |
    if [ ! -f "$WASP2_TEST_DATA_CACHE/sample.bam" ]; then
      ./scripts/download_test_data.sh
    fi
```

## Security Considerations

### Repository Access

Self-hosted runners have full access to repository secrets. For public repositories:

1. **Never run untrusted code** on self-hosted runners
2. Configure runners to only accept jobs from protected branches
3. Use repository environments with required reviewers

### Network Security

```bash
# Restrict outbound connections (example with iptables)
sudo iptables -A OUTPUT -m owner --uid-owner github-runner -d github.com -j ACCEPT
sudo iptables -A OUTPUT -m owner --uid-owner github-runner -d api.github.com -j ACCEPT
sudo iptables -A OUTPUT -m owner --uid-owner github-runner -d pypi.org -j ACCEPT
# Add other required destinations...
```

## Monitoring and Maintenance

### Health Checks

Create a simple health check script:

```bash
#!/bin/bash
# /usr/local/bin/check-runner-health.sh

RUNNER_DIR=/home/github-runner/actions-runner

# Check if runner process is running
if ! pgrep -f "Runner.Listener" > /dev/null; then
    echo "Runner not running, restarting..."
    sudo systemctl restart actions.runner.mcvickerlab-WASP2.hpc-runner-01.service
fi

# Check disk space
DISK_USAGE=$(df /data | tail -1 | awk '{print $5}' | sed 's/%//')
if [ "$DISK_USAGE" -gt 90 ]; then
    echo "Warning: Disk usage at ${DISK_USAGE}%"
    # Clean old work directories
    find "${RUNNER_DIR}/_work" -type d -mtime +7 -exec rm -rf {} \;
fi
```

### Log Rotation

Add to `/etc/logrotate.d/github-runner`:

```
/home/github-runner/actions-runner/_diag/*.log {
    daily
    rotate 7
    compress
    missingok
    notifempty
}
```

## Troubleshooting

### Common Issues

**Runner not picking up jobs:**
```bash
# Check runner status
sudo systemctl status actions.runner.*.service

# Check connectivity
curl -s https://api.github.com/zen
```

**Build failures due to missing tools:**
```bash
# Verify environment
which python3 rustc cargo bcftools samtools bedtools

# Check versions
python3 --version
rustc --version
```

**Out of memory errors:**
```bash
# Monitor during job execution
watch -n 1 'free -h; echo; top -b -n 1 | head -20'
```

## Additional Resources

- [GitHub Actions Self-Hosted Runners Documentation](https://docs.github.com/en/actions/hosting-your-own-runners)
- [WASP2 Development Guide](./CONTRIBUTING.md)
- [HPC Best Practices for Genomics](https://hpc.nih.gov/apps/WASP.html)
