"""
WASP2 Configuration Management.

Handles loading and saving user configuration from ~/.wasp2/config.yaml.
Configuration values are merged with CLI arguments, with CLI taking precedence.
"""

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional
import logging
import os

# Optional YAML support - gracefully degrade if not available
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False


# Default configuration directory
CONFIG_DIR = Path.home() / ".wasp2"
CONFIG_FILE = CONFIG_DIR / "config.yaml"


@dataclass
class WASP2Config:
    """WASP2 configuration settings.

    Attributes:
        threads: Default number of threads for parallel operations
        use_rust: Whether to use Rust acceleration by default
        temp_dir: Default temporary directory (None = system default)
        output_format: Default output format ('tsv', 'json', 'yaml')
        log_level: Logging level ('DEBUG', 'INFO', 'WARNING', 'ERROR')
        log_file: Optional log file path
        include_indels: Whether to include indels by default
        pseudocount: Default pseudocount for imbalance analysis
        min_count: Default minimum count threshold
    """
    threads: int = 1
    use_rust: bool = True
    temp_dir: Optional[str] = None
    output_format: str = "tsv"
    log_level: str = "WARNING"
    log_file: Optional[str] = None
    include_indels: bool = False
    pseudocount: int = 1
    min_count: int = 10

    def to_dict(self) -> dict:
        """Convert config to dictionary, excluding None values."""
        return {k: v for k, v in asdict(self).items() if v is not None}


def get_config_dir() -> Path:
    """Get the WASP2 configuration directory, creating if needed."""
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    return CONFIG_DIR


def load_config() -> WASP2Config:
    """Load configuration from file, returning defaults if not found."""
    if not CONFIG_FILE.exists():
        return WASP2Config()

    if not YAML_AVAILABLE:
        logging.debug("PyYAML not installed, using default config")
        return WASP2Config()

    try:
        with open(CONFIG_FILE) as f:
            data = yaml.safe_load(f) or {}

        # Filter to only valid fields
        valid_fields = {f.name for f in WASP2Config.__dataclass_fields__.values()}
        filtered = {k: v for k, v in data.items() if k in valid_fields}

        return WASP2Config(**filtered)
    except Exception as e:
        logging.warning(f"Error loading config file: {e}. Using defaults.")
        return WASP2Config()


def save_config(config: WASP2Config) -> Path:
    """Save configuration to file."""
    if not YAML_AVAILABLE:
        raise RuntimeError("PyYAML is required to save config. Install with: pip install pyyaml")

    get_config_dir()  # Ensure directory exists

    with open(CONFIG_FILE, "w") as f:
        yaml.safe_dump(config.to_dict(), f, default_flow_style=False, sort_keys=False)

    return CONFIG_FILE


def get_config_path() -> Path:
    """Get the path to the config file."""
    return CONFIG_FILE


# Verbosity level mapping for CLI
VERBOSITY_LEVELS = {
    0: logging.WARNING,   # Default: warnings and errors only
    1: logging.INFO,      # -v: info messages
    2: logging.DEBUG,     # -vv: debug messages
    3: logging.DEBUG,     # -vvv: debug with extra detail (same level, different behavior)
}


def setup_logging(verbosity: int = 0, log_file: Optional[str] = None) -> None:
    """Configure logging based on verbosity level.

    Args:
        verbosity: 0 = WARNING, 1 = INFO, 2+ = DEBUG
        log_file: Optional file to write logs to
    """
    level = VERBOSITY_LEVELS.get(verbosity, logging.DEBUG)

    # Format varies by verbosity
    if verbosity >= 2:
        fmt = "%(asctime)s [%(levelname)s] %(name)s:%(lineno)d - %(message)s"
    elif verbosity >= 1:
        fmt = "%(levelname)s: %(message)s"
    else:
        fmt = "%(message)s"

    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=level,
        format=fmt,
        handlers=handlers,
        force=True,  # Override any existing configuration
    )


# Global config instance (lazy-loaded)
_config: Optional[WASP2Config] = None


def get_config() -> WASP2Config:
    """Get the global configuration, loading from file if needed."""
    global _config
    if _config is None:
        _config = load_config()
    return _config


def reset_config() -> None:
    """Reset the global config to force reload."""
    global _config
    _config = None
